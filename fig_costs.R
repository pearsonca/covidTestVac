
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(emojifont)
  require(showtext)
})

grid.dt <- as.data.table(expand.grid(
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

grid.dt[, protect_mul := expanded_protection_factor(sens, spec, seropos) ]

grid.dt[, cost_mul := cost_per_dose(sens, spec, seropos, test_cost_frac)/protect_mul - 1 ]

shrd <- list(
  scale_x_continuous("Test Sensitivity (True Positive Rate)"),
  scale_y_continuous("Test Specificity (True Negative Rate)"),
  theme_minimal(base_size = 18),
  theme(
    legend.position = "top",
    legend.box.margin = margin(b = -20),
    legend.margin = margin(),
    legend.key.height = unit(10, "pt"),
    legend.text = element_text(margin = margin(t=-15.5), face = "bold")
  ),
  coord_equal(expand = F)
)

refdim <- 2.5

p.people <- ggplot(grid.dt[test_cost_frac == 1.0]) + aes(sens, spec, fill = (protect_mul-1)*100) +
  facet_grid(
    . ~ seropos,
    labeller = labeller(seropos = function(s) sprintf("Sero+ = %s%%", as.numeric(s)*100))
  ) +
  geom_raster() +
  scale_fill_gradient2(
    "% change in people protected",
    midpoint = 0, limits = c(-100, 100), breaks = seq(-75, 75, by=25)
  ) +
  shrd

ggsave(
  "people.png",
   p.people,
   device = "png",
   width = grid.dt[, length(unique(seropos))]*refdim*2/3,
   height = refdim,
   units = "in"
)

pcost <- ggplot(grid.dt) + aes(sens, spec, fill = (cost_mul-1)*100) +
  facet_grid(
    test_cost_frac ~ seropos,
    labeller = labeller(
      seropos = function(s) sprintf("Sero+ = %s%%", as.numeric(s)*100),
      test_cost_frac = function(s) sprintf("T/V = %s", as.numeric(s))
    )
  ) +
  geom_raster() +
  scale_fill_gradient2(
    "% change in cost per person protected",
    breaks = c(-100, 0, 100),
    midpoint = 0,
    low = scales::muted("blue"), high = scales::muted("red")
  ) +
  shrd

ggsave("costs.png",
       pcost,
       width = grid.dt[, length(unique(seropos))]*refdim,
       height = grid.dt[, length(unique(test_cost_frac))],
       dpi = 900, units = "in"
)

ref.coverage <- .2
people.covered <- round(100*ref.coverage)
ref.sero.pos <- .4

waffle.dt <- as.data.table(expand.grid(
  x = 1:10, y = 1:10
))[, p := runif(.N) ][,
  col := c("sero-","sero+")[(p <= ref.sero.pos)+1]
]

waffle.dt[, vaccinated := NA_character_ ]
waffle.dt[1:people.covered, vaccinated := "no-test" ]

add70 <- round(people.covered*expanded_protection_factor(.7, 1, ref.sero.pos))
add90 <- round(people.covered*expanded_protection_factor(.9, 1, ref.sero.pos))

waffle.dt[1:add70, vaccinated := fifelse(is.na(vaccinated),"test-70",vaccinated)]
waffle.dt[1:add90, vaccinated := fifelse(is.na(vaccinated),"test-90",vaccinated)]

serocols <- c("sero+"="dodgerblue", "sero-"="firebrick")

p.prot <- ggplot(waffle.dt) + 
  aes(
    color = col, x = y, y = x, label = fontawesome("fa-user")
  ) +
  geom_raster(
    aes(alpha = vaccinated, fill = !is.na(vaccinated))
  ) +
  geom_text(
    family = "fontawesome-webfont",
    size = 9, key_glyph = "rect"
  ) +
  scale_color_manual(
    name = "Antibody Status",
    labels = c("sero+"="Seropositive", "sero-"="Seronegative"),
    values = serocols
  ) +
  scale_alpha_manual(
    name = "Individuals Protected",
    breaks = c("no-test", "test-70", "test-90"),
    labels = c(
      "no-test"="Coverage w/o Testing",
      "test-70"="w/ 70% Sensitivity Test",
      "test-90"="w/ 90% Sensitivity Test"
    ),
    values = rev(c(1, .7, .5)),
    guide = guide_legend(override.aes = list(fill = "grey18"))
  ) +
  scale_fill_manual(guide = "none", values = c(NA, "grey18")) +
  coord_equal() +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )

ggsave(
  "protection.png",
  p.prot, width = 7, height = 5,
  dpi = 900,
  units = "in"
)
  