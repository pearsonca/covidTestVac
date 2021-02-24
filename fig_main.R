
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
})

.args <- if (interactive()) c(
  "main.png"
) else commandArgs(trailingOnly = TRUE)

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

lblsize <- 3
lnsize <- 1

p.coverage <- ggplot(ref.mlt) +
  aes(seropos, value, group = ord, alpha = serostatus, fill = scenario) +
  geom_area() +
  coord_cartesian(xlim = c(0.2, 0.8), ylim = c(0, 0.4), expand = FALSE) +
  scale_alpha_manual(
    values = c("sn"=.5, "sp"=.75),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c(notest="grey5", test70=blues9[5], test90=blues9[9]),
    guide = "none"
  ) +
  geom_line(aes(seropos, y=covaxcoverage, linetype = "notest"), ref, inherit.aes = FALSE, show.legend = F, size = lnsize) +
  geom_line(aes(seropos, y=(costlim70+1)*covaxcoverage, linetype = "test70"), ref, inherit.aes = FALSE, show.legend = F, size = lnsize) +
  geom_line(aes(seropos, y=(costlim90+1)*covaxcoverage, linetype = "test90"), ref, inherit.aes = FALSE, show.legend = F, size = lnsize) +
  annotate(
    "label",
    x=0.5, y= covaxcoverage * .95,
    label = expression(symbol('\257')*"  "*symbol('\257')*'  Coverage without Testing  '*symbol('\257')*"  "*symbol('\257')),
    vjust = "top",
    label.size = 0, fill = alpha(c("white"), 0.5),
    size = lblsize
  ) +
  annotate(
    "label",
    x = (0.8-0.2)*.95 + 0.2, y = (PPDmul(1, .7, 0.5)*covaxcoverage + covaxcoverage)/2,
    label = 'Additional Coverage with 70% Sensitive Test',
    hjust = "right",
    label.size = 0, fill = alpha(c("white"), 0.5),
    size = lblsize
  ) +
  annotate(
    "label",
    x = (0.8-0.2)*.95 + 0.2,
    y = (PPDmul(1, .9, 0.8)+PPDmul(1, .7, 0.8))*covaxcoverage/2,
    label = '...with 90% Sensitive Test',
    hjust = "right",
    label.size = 0, fill = alpha(c("white"), 0.5),
    size = lblsize
  ) +
  annotate(
    "text",
    x = (0.8-0.2)*.05 + .2,
    y = covaxcoverage/2,
    label = 'Seronegative\nIndividuals',
    hjust = "left",
    color = "black",
    size = lblsize
  ) +
  annotate(
    "text",
    x = (0.8-0.2)*.95 + .2,
    y = covaxcoverage/2,
    label = 'Seropositive\nIndividuals',
    hjust = "right",
    color = "white",
    size = lblsize
  ) +
  scale_y_continuous("Total Immunized %", breaks = seq(0,0.5,0.1), minor_breaks = NULL, labels = function(b) sprintf("%i%%", as.integer(b*100))) +
  scale_x_continuous("% Seropositive", minor_breaks = NULL, labels = function(b) sprintf("%i%%", as.integer(b*100))) +
  scale_linetype_manual(
    values = c(notest="solid", test70="longdash", test90="dotted")
  ) +
  theme_minimal() +
  theme(
    plot.margin = margin(r = 1, unit = "line")
  )

p.cost <- ggplot(ref) +
  aes(seropos) +
  geom_blank(aes(linetype = "notest")) +
  geom_ribbon(aes(ymin=0, ymax=costlim70), fill = alpha(blues9[5], .8), show.legend = F) +
  geom_ribbon(aes(ymin=costlim70, ymax=costlim90), fill = alpha(blues9[9], .8), show.legend = F) +
  geom_line(aes(y=costlim70, linetype = "test70"), size = lnsize) +
  geom_line(aes(y=costlim90, linetype = "test90"), size = lnsize) +
  annotate(
    "label",
    x = (0.8-0.2)*.975 + .2,
    y = (PPDmul(1, .7, .6) - 1)/3*2,
    label = '...if sensitivity >70%',
    hjust = "right",
    label.size = 0, fill = alpha(c("white"), 0.5),
    vjust = "top",
    size = lblsize
  ) +
  annotate(
    "label",
    x = (0.8-0.2)*.975 + .2,
    y = (PPDmul(1, .7, .8) + PPDmul(1, .9, .8) - 2)/2,
    label = 'buy tests\nif sensitivity > 90%',
    hjust = "right",
    label.size = 0, fill = alpha(c("white"), 0.5),
    size = lblsize
  ) +
  annotate(
    "text",
    x = .3,
    y = .75+.125,
    label = 'buy vaccine doses',
    hjust = "left",
    size = lblsize
  ) +
  coord_cartesian(xlim = c(0.2, 0.8), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous("% Seropositive", minor_breaks = NULL, labels = function(b) sprintf("%i%%", as.integer(b*100))) +
  scale_y_continuous("Test Cost % of Vaccine Dose Cost", labels = function(b) sprintf("%i%%", as.integer(b*100))) +
  scale_linetype_manual(
    "Test Performance",
    breaks = c("notest", "test70", "test90"),
    labels = c(
      "notest"="(No Test)",
      "test70"="70% Sensitivity",
      "test90"="90% Sensitivity"
    ),
    values = c("solid", "longdash", "dotted"),
    guide = "none"
  ) +
  theme_minimal() +
  theme(
    plot.margin = margin(r = 1, unit = "line")
  )

res.p <- p.coverage / p.cost +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave(tail(.args, 1), res.p, height = 7.5, width = 4, dpi = 600)
