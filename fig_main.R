
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

#' panel A: with increasing seroprevalence, increasing
#' coverage enabled by testing

PPDmul <- function(TNR, TPR, seropos) TNR + seropos*(1+TPR-TNR)

ref <- data.table(seropos = seq(0,1,by=0.01))

covaxcoverage <- .2

ref[, notest_sp := covaxcoverage*seropos ]
ref[, notest_sn := covaxcoverage*(1-seropos) ]
ref[, test70_sp := notest_sp*(PPDmul(1, .7, seropos) - 1) ]
ref[, test70_sn := notest_sn*(PPDmul(1, .7, seropos) - 1) ]
ref[, test90_sp := notest_sp*(PPDmul(1, .9, seropos)-1) - test70_sp ]
ref[, test90_sn := notest_sn*(PPDmul(1, .9, seropos)-1) - test70_sn ]

ref.mlt <- melt.data.table(
  ref,
  id.vars = "seropos"
)[, c("scenario", "serostatus") := tstrsplit(variable, "_")]

ref.mlt[,
  ord := 2*c("notest"=2,"test70"=1,"test90"=0)[scenario]+c(sp=0, sn=1)[serostatus]
]

ref[, costlim70 := PPDmul(1, .7, seropos) - 1 ]
ref[, costlim90 := PPDmul(1, .9, seropos) - 1 ]

p.coverage <- ggplot(ref.mlt) +
  aes(seropos, value, group = ord, fill = serostatus, alpha = scenario) +
  geom_area() +
  coord_cartesian(xlim = c(0.2, 0.8), ylim = c(0, 0.5), expand = FALSE) +
  scale_alpha_manual(
    "Test Performance",
    labels = c(
      "notest"="(No Test)",
      "test70"="70% Sensitivity",
      "test90"="90% Sensitivity"
    ),
    values = c("notest"=.25, "test70"=.5, "test90"=.9),
    guide = guide_legend(override.aes = list(fill="grey5"))
  ) +
  scale_fill_manual(
    "Vaccinee Serostatus",
    breaks = c("sp","sn"),
    labels = c(sn="Negative", sp="Positive"),
    values = c(sn="dodgerblue", sp="firebrick")
  ) +
  geom_line(aes(seropos, y=covaxcoverage, linetype = "notest"), ref, inherit.aes = FALSE, show.legend = F) +
  geom_line(aes(seropos, y=(costlim70+1)*covaxcoverage, linetype = "test70"), ref, inherit.aes = FALSE, show.legend = F) +
  geom_line(aes(seropos, y=(costlim90+1)*covaxcoverage, linetype = "test90"), ref, inherit.aes = FALSE, show.legend = F) +
  scale_y_continuous("Total Immunized Fraction", breaks = seq(0,0.5,0.1), minor_breaks = NULL) +
  scale_x_continuous("Seroprevalence") +
  scale_linetype_manual(
    values = c(notest="solid", test70="longdash", test90="dotted")
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1)
  )

p.cost <- ggplot(ref) +
  aes(seropos) +
  geom_blank(aes(linetype = "notest")) +
  geom_ribbon(aes(ymin=0, ymax=costlim70, alpha = "test70"), fill = "grey5", show.legend = F) +
  geom_ribbon(aes(ymin=costlim70, ymax=costlim90, alpha = "test90"), fill = "grey5", show.legend = F) +
  geom_line(aes(y=costlim70, linetype = "test70")) +
  geom_line(aes(y=costlim90, linetype = "test90")) +
  coord_cartesian(xlim = c(0.2, 0.8), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous("Seroprevalence") +
  scale_y_continuous("Cost-Limit for Testing,\nas a Fraction of Dose Cost") +
  scale_linetype_manual(
    "Test Performance",
    breaks = c("notest", "test70", "test90"),
    labels = c(
      "notest"="(No Test)",
      "test70"="70% Sensitivity",
      "test90"="90% Sensitivity"
    ),
    values = c("solid", "longdash", "dotted")
  ) +
  scale_alpha_manual(values = c("test90"=.9, "test70"=.5)) +
  theme_minimal() +
  theme(
    legend.position = c(0, 1), legend.justification = c(0, 1)
  )

res.p <- p.coverage / p.cost +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

# ref[, test70_sp := notest_sp*PPDmul(1, .7, seropos) ]
# ref[, test70_sn := notest_sn*PPDmul(1, .7, seropos) ]
# ref[, test90_sp := notest_sp*PPDmul(1, .9, seropos) ]
# ref[, test90_sn := notest_sn*PPDmul(1, .9, seropos) ]

#' want rectangle for the no-test coverage, sn expanding below 0,
#' sp expanding above rectangle

ref[, c(
  "allsn", "sn90", "sn70", "sn",
  "sp", "sp70", "sp90", "allsp"
) := .(
  notest_sn - (1-seropos),
  notest_sn - test90_sn,
  notest_sn - test70_sn,
  notest_sn,
  notest_sn + notest_sp,
  notest_sn + test70_sp,
  notest_sn + test90_sp,
  seropos
)]

ref.mlt <- melt.data.table(
  ref, id.vars = "seropos"
)

ggplot(ref) + aes(seropos) +
  geom_ribbon(aes(
    ymin = allsn, ymax = sn90,
    fill = "seronegative", alpha = "all")
  ) +
  geom_ribbon(aes(
    ymin=sn90, ymax=sn70,
    fill = "seronegative", alpha = "90"
  )) +
  geom_ribbon(aes(
    ymin=sn70, ymax = 0,
    fill = "seronegative", alpha = "70"
  )) +
  geom_ribbon(aes(
    ymax=allsp, ymin = sp90,
    fill = "seropositive", alpha = "all"
  )) +
  geom_ribbon(aes(
    ymax = sp90, ymin = sp70,
    fill = "seropositive", alpha = "90"
  )) +
  geom_ribbon(aes(
    ymax = sp70, ymin = sp,
    fill = "seropositive", alpha = "70"
  )) +
  geom_ribbon(aes(
    ymax = sp, ymin = sn,
    fill = "seropositive", alpha = "notest"
  )) +
  geom_ribbon(aes(
    ymax = sn, ymin = 0,
    fill = "seronegative", alpha = "notest"
  )) +
  coord_cartesian(
    xlim = c(.2, .8),
    ylim = c(-.4, .5),
    expand = FALSE
  ) +
  scale_alpha_manual(
    name = NULL,
    breaks = c("all", "notest", "70", "90"),
    values = c("all"=.1, "90"=.25, "70"=.6, "notest" = 1)
  ) +
  scale_y_continuous(NULL, breaks = 0.1+seq(-.4,.4,.1), minor_breaks = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank()
  )
