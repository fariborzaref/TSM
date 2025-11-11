############################################################
# OECD PANEL TIME SERIES (2002 to 2021)
# Author: Dr. Fariborz Aref
# Focus: Income (gini_income), Health (health_ineq), Labor (labor_ineq)
# Methods: Panel unit root tests, Arellano–Bond GMM, Panel VAR with IRFs
############################################################

# 0) Packages
required <- c(
  "tidyverse","plm","sandwich","lmtest","panelvar","data.table","janitor","ggplot2","patchwork"
)
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))

set.seed(2025)
options(stringsAsFactors = FALSE)

# Compact academic theme
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 11, base_family = "serif") +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      axis.title   = ggplot2::element_text(size = 10),
      axis.text    = ggplot2::element_text(size = 9),
      strip.text   = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
)

# 1) Load data
csv_path <- "oecd_inequality_2002_2021.csv"   # place this next to tsm_run.R or set full path
if (!file.exists(csv_path)) stop("File not found: ", csv_path)

dat <- data.table::fread(csv_path) |>
  janitor::clean_names()

need <- c("country","year","gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
miss <- setdiff(need, names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

dat <- dat |>
  dplyr::filter(year >= 2002, year <= 2021) |>
  dplyr::mutate(country = as.factor(country)) |>
  dplyr::arrange(country, year)

# Panel balance check
panel_len <- dat |>
  dplyr::count(country) |>
  dplyr::pull(n)
if (length(unique(panel_len)) != 1) {
  warning("Panel is not balanced. Methods used tolerate some unbalanced panels.")
}

# Convert to pdata.frame
pdat <- plm::pdata.frame(dat, index = c("country","year"))

# 2) Panel unit root diagnostics
vars_to_test <- c("gini_income","health_ineq","labor_ineq","gdp_pc","unemployment","openness")
cat("\n=== PANEL UNIT ROOT TESTS (LLC and IPS) ===\n")
ur_list <- lapply(vars_to_test, function(v) {
  x <- pdat[[v]]
  list(
    var = v,
    LLC = plm::purtest(x, pdat, test = "levinlin", exo = "trend", lags = "AIC"),
    IPS = plm::purtest(x, pdat, test = "ips",      exo = "trend", lags = "AIC")
  )
})
for (res in ur_list) {
  cat("\n---", res$var, "---\n"); print(summary(res$LLC)); print(summary(res$IPS))
}

# Note: if clear nonstationarity is detected, consider differencing, PMG, or ARDL.
# Here we continue with difference GMM and de-meaned PVAR.

# 3) Arellano–Bond difference GMM for dynamic income inequality
form_ab <- gini_income ~ lag(gini_income, 1) + health_ineq + labor_ineq +
  gdp_pc + unemployment + openness

ab_gmm <- plm::pgmm(
  formula = form_ab,
  data = pdat,
  effect = "twoways",       # country and time effects
  model  = "twosteps",      # two step GMM
  transformation = "d",     # first differences
  collapse = TRUE
)

cat("\n=== ARELLANO–BOND GMM: gini_income ===\n")
print(lmtest::coeftest(ab_gmm, vcov = plm::vcovHC(ab_gmm, method = "arellano", type = "HC0")))

cat("\n-- AB tests --\n")
print(plm::mtest(ab_gmm, order = 1))   # AR(1)
print(plm::mtest(ab_gmm, order = 2))   # AR(2) ideally p > 0.05

cat("\n-- Sargan/Hansen --\n")
print(plm::sargan(ab_gmm))

# 4) Panel VAR with fixed effects
pvar_vars <- c("gini_income","health_ineq","labor_ineq")
pvar_data <- pdat[, c("country","year", pvar_vars)]
pvar_data <- stats::na.omit(pvar_data)

pvar_fit <- panelvar::pvarfeols(
  dependent_vars = pvar_vars,
  lags = 2,
  transformation = "demean",
  data = as.data.frame(pvar_data),
  panel_identifier = c("country","year")
)

cat("\n=== PVAR summary ===\n")
print(summary(pvar_fit))

# 5) IRFs with bootstrap
set.seed(2025)
irf_res <- panelvar::pvarirf(
  x = pvar_fit,
  impulse = pvar_vars,
  response = pvar_vars,
  n.ahead = 8,
  ci = 0.95,
  runs = 500,
  cumulative = FALSE,
  ortho = TRUE
)

cat("\nIRFs computed. Use plot_irf to visualize and save.\n")

# 6) Plot IRFs and save a grid for gini responses
plot_irf <- function(irf_list, impulse, response) {
  df  <- as.data.frame(irf_list$irf[[impulse]][[response]]); names(df) <- "irf"
  df$h <- seq_len(nrow(df)) - 1
  lo  <- as.data.frame(irf_list$Lower[[impulse]][[response]]); names(lo) <- "lower"
  up  <- as.data.frame(irf_list$Upper[[impulse]][[response]]); names(up) <- "upper"
  df$lower <- lo$lower; df$upper <- up$upper

  ggplot2::ggplot(df, ggplot2::aes(x = h, y = irf)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.15) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::labs(title = paste0("IRF: shock to ", impulse, " to ", response),
                  x = "Horizon (years)", y = "Response")
}

g1 <- plot_irf(irf_res, "gini_income",  "gini_income")
g2 <- plot_irf(irf_res, "health_ineq",  "gini_income")
g3 <- plot_irf(irf_res, "labor_ineq",   "gini_income")

irf_grid <- g1 + g2 + g3 + patchwork::plot_layout(ncol = 3)

if (!dir.exists("TSM/figs")) dir.create("TSM/figs", recursive = TRUE)
ggplot2::ggsave("TSM/figs/irf_gini_grid.png", irf_grid, width = 9.5, height = 3.6, dpi = 300)

# 7) Interpretation hooks
cat("\n=== INTERPRETATION GUIDE ===\n")
cat("* Arellano–Bond: lag coefficient on gini reflects persistence. Significant health_ineq or labor_ineq implies cross domain effects.\n")
cat("* Diagnostics: AR(2) near null and acceptable Sargan or Hansen suggest instrument validity.\n")
cat("* PVAR IRFs: check sign, size, and duration of responses of gini to health and labor shocks for policy relevance.\n")

# 8) Save artifacts
if (!dir.exists("TSM/out"))  dir.create("TSM/out",  recursive = TRUE)
saveRDS(list(
  ab_gmm = ab_gmm,
  pvar   = pvar_fit,
  irf    = irf_res
), file = "TSM/out/tsm_oecd_models.rds")

cat("\nDone. Artifacts saved in TSM/out and figures in TSM/figs.\n")
