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

saveRDS(ret, tail(.args, 1))