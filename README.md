############################################################
# OECD Panel Time-Series (2002–2021)
# Author: Dr. Fariborz Aref
# Focus: Income (gini_income), Health (health_ineq),
#        and Labor (labor_ineq) inequality dynamics
# Methods: Panel unit-root tests, Arellano–Bond (GMM), Panel VAR + IRFs
############################################################

# ---- 0) Packages ----
required <- c("tidyverse","plm","sandwich","lmtest","panelvar","data.table")
to_install <- setdiff(required, rownames(installed.packages()))
if(length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

library(tidyverse)
library(plm)
library(sandwich)
library(lmtest)
library(panelvar)
library(data.table)

# ---- 1) Load & Inspect Data ----
# Place the CSV in repo root: oecd_inequality_2002_2021.csv
dat <- fread("oecd_inequality_2002_2021.csv") |>
  janitor::clean_names()

# Basic checks
stopifnot(all(c("country","year","gini_income","health_ineq","labor_ineq",
                "gdp_pc","unemployment","openness") %in% names(dat)))

# Keep 2002–2021 for 27 countries
dat <- dat |>
  filter(year >= 2002, year <= 2021) |>
  mutate(country = as.factor(country)) |>
  arrange(country, year)

# Balanced panel check
panel_idx <- dat |>
  count(country) |>
  pull(n)
if(length(unique(panel_idx)) != 1) {
  warning("Panel is not balanced. Proceeding with available balance. Consider balancing or using methods robust to unbalanced panels.")
}

# Convert to pdata.frame
pdat <- pdata.frame(dat, index = c("country","year"))

# ---- 2) Panel Unit-Root Diagnostics (Levin-Lin-Chu, IPS) ----
vars_to_test <- c("gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
unitroot_results <- lapply(vars_to_test, function(v) {
  x <- pdat[[v]]
  list(
    var = v,
    LLC = purtest(x, pdat, test = "levinlin", exo = "trend", lags = "AIC"),
    IPS = purtest(x, pdat, test = "ips", exo = "trend", lags = "AIC")
  )
})
cat("\n=== PANEL UNIT ROOT TESTS (LLC & IPS) ===\n")
for (res in unitroot_results) {
  cat("\n---", res$var, "---\n")
  print(summary(res$LLC))
  print(summary(res$IPS))
}

# (If nonstationary, consider first-differencing or PMG/ARDL; here we proceed with
# difference GMM in dynamic model and de-meaned PVAR that tolerate some I(1) behavior.)

# ---- 3) Dynamic Panel: Arellano–Bond Difference GMM (pgmm) ----
# Model: gini_income_t = α*gini_{t-1} + β1*health_ineq + β2*labor_ineq + β3*gdp_pc
#        + β4*unemployment + β5*openness + country FE + time FE
# Instruments: deeper lags of gini_income in differences
form_ab <- gini_income ~ lag(gini_income, 1) + health_ineq + labor_ineq + gdp_pc + unemployment + openness

ab_gmm <- pgmm(
  formula = form_ab,
  data = pdat,
  effect = "twoways",       # country & time effects
  model  = "twosteps",      # Windmeijer-corrected SE advisable (on twostep vcov)
  transformation = "d",     # first differences (Arellano–Bond)
  collapse = TRUE,
  # Instruments: gini_income in levels as IVs for differenced lag
  # plm::pgmm handles lag depth with 'lag.form' inside formula; we use default deep lags
  # If you want explicit: ~ lag(gini_income, 2:5) as instruments
)

cat("\n=== ARELLANO–BOND GMM: gini_income ===\n")
coeftest(ab_gmm, vcov = vcovHC(ab_gmm, method = "arellano", type = "HC0"))

# Diagnostics: AR(1), AR(2) in differences; Hansen/Sargan tests
cat("\n-- AB Tests --\n")
print(mtest(ab_gmm, order = 1))  # AR(1)
print(mtest(ab_gmm, order = 2))  # AR(2) should be > 0.05 ideally
cat("\n-- Sargan/Hansen --\n")
print(sargan(ab_gmm))

# ---- 4) Panel VAR (PVAR) with {panelvar} ----
# Variables: gini_income, health_ineq, labor_ineq
# We demean (fixed effects) and use 2 lags (adjust if T small/large)
PVAR_VARS <- c("gini_income","health_ineq","labor_ineq")

# Ensure complete cases for PVAR on selected vars
pvar_data <- pdat[, c("country","year", PVAR_VARS)]
pvar_data <- na.omit(pvar_data)

# Choose lags (can compare AIC/BIC manually; here set to 2 for medium T)
pvar_fit <- pvarfeols(
  dependent_vars = PVAR_VARS,
  lags = 2,
  transformation = "demean",
  data = as.data.frame(pvar_data),
  panel_identifier = c("country","year")
)

cat("\n=== PVAR (FE-OLS) Summary ===\n")
print(summary(pvar_fit))

# Robust covariance & Wald tests per equation
cat("\n=== Equation-wise Robust Inference ===\n")
for (dv in PVAR_VARS) {
  eq <- pvar_fit$coefficients[[dv]]
  cat("\n--- Dependent:", dv, "---\n")
  print(eq)
}

# ---- 5) Impulse Response Functions (IRFs) ----
# Bootstrap IRFs for shocks to each variable (H steps ahead)
set.seed(2025)
irf_res <- pvarirf(
  x = pvar_fit,
  impulse = PVAR_VARS, response = PVAR_VARS,
  n.ahead = 8, ci = 0.95, runs = 500, cumulative = FALSE,
  ortho = TRUE # Cholesky identification; order matters
)

cat("\n=== IRFs computed (store & plot as needed) ===\n")

# ---- 6) Plot IRFs (ggplot) ----
plot_irf <- function(irf_list, impulse, response) {
  df <- as.data.frame(irf_list$irf[[impulse]][[response]])
  names(df) <- "irf"
  df$h <- seq_len(nrow(df)) - 1
  # CIs
  lower <- as.data.frame(irf_list$Lower[[impulse]][[response]]); names(lower) <- "lower"
  upper <- as.data.frame(irf_list$Upper[[impulse]][[response]]); names(upper) <- "upper"
  df$lower <- lower$lower
  df$upper <- upper$upper

  ggplot(df, aes(x = h, y = irf)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
    geom_line(size = 0.9) +
    labs(
      title = paste0("IRF: Shock to ", impulse, " → ", response),
      x = "Horizon (years)", y = "Response"
    ) +
    theme_minimal(base_size = 12)
}

# Example: response of gini to shocks in each variable
# g1 <- plot_irf(irf_res, impulse = "gini_income", response = "gini_income")
# g2 <- plot_irf(irf_res, impulse = "health_ineq", response = "gini_income")
# g3 <- plot_irf(irf_res, impulse = "labor_ineq", response = "gini_income")
# print(g1); print(g2); print(g3)

# ---- 7) Policy-Style Interpretation Scaffold (print hooks) ----
cat("\n=== INTERPRETATION GUIDE ===\n")
cat("* Arellano–Bond: α on lag(gini) captures persistence of income inequality;\n")
cat("  significant β on health_ineq/labor_ineq implies cross-domain influence on income inequality.\n")
cat("* AR(2) p>0.05 + valid Hansen/Sargan → instrument set plausible.\n")
cat("* PVAR IRFs: trace how shocks to health or labor inequality propagate into gini over 2–8 years;\n")
cat("  look for sign, magnitude, and duration to inform OECD-level policy levers.\n")

# ---- 8) Save Model Artifacts (optional) ----
saveRDS(list(ab_gmm = ab_gmm, pvar = pvar_fit, irf = irf_res), file = "TSM/model_artifacts_oecd_2002_2021.rds")

############################################################
# End of script
############################################################
