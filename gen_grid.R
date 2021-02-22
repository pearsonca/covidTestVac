suppressPackageStartupMessages({
  require(data.table)
})

.args <- if (interactive()) c(
  "grid.rds"
) else commandArgs(trailingOnly = TRUE)

ret <- as.data.table(expand.grid(
  sens = seq(0.7, 1, by=0.01),
  spec = seq(0.7, 1, by=0.01),
  seropos = seq(0.2, 0.8, by=0.2),
  test_cost_frac = seq(0.2, 0.8, by=0.1)
))

expanded_protection_factor <- function(
  TPR, TNR, seropositivity
) TNR + seropositivity*(1 + TPR - TNR)

cost_per_dose <- function(
  TPR, TNR, seropositivity, test_cost_fraction
) 1 + test_cost_fraction/(1+(1-seropositivity)*TNR + seropositivity*(1-TPR))

ret[, protect_mul := expanded_protection_factor(sens, spec, seropos) ]

ret[, cost_mul := cost_per_dose(sens, spec, seropos, test_cost_frac)/protect_mul - 1 ]

saveRDS(ret, tail(.args, 1))