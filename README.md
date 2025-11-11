# OECD Panel Time Series, 2002 to 2021

**Author:** Dr. Fariborz Aref  
**Discipline:** Quantitative Sociology and Inequality Research  
**License:** MIT  

### Purpose
Dynamic panel analysis of inequality across OECD members. Implements panel unit root tests, Arellano–Bond difference GMM for income inequality persistence, and Panel VAR with impulse responses for cross domain effects among income, health, and labor inequality.

### Structure
TSM/
├── tsm_run.R
├── oecd_inequality_2002_2021.csv
├── out/
│ └── tsm_oecd_models.rds
└── figs/
└── irf_gini_grid.png

### Methods
- Panel unit root tests Levin Lin Chu and IPS  
- Arellano–Bond difference GMM with two way effects  
- Panel VAR with fixed effects and two lags  
- Bootstrap impulse responses with orthogonal identification

### Quick example
```r
library(plm)
ab <- pgmm(gini_income ~ lag(gini_income, 1) + health_ineq + labor_ineq + gdp_pc + unemployment + openness,
           data = pdata.frame(dat, index = c("country","year")),
           effect = "twoways", model = "twosteps", transformation = "d", collapse = TRUE)
summary(ab)

