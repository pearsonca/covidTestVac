
suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(patchwork)
  require(magrittr)
  require(dplyr)
  require(scales)
  require(forcats)
  require(tibble)
})

.args <- if (interactive()) c(
  "main.png"
) else commandArgs(trailingOnly = TRUE)

PPDmul <- function(TNR, TPR, seropos) TNR + seropos*(1+TPR-TNR)

ref <- data.table(seropos = seq(0,1,by=0.01))

covaxcoverage <- .2

refspec <- .99

ref[, notest_sp := covaxcoverage*seropos ]
ref[, notest_sn := covaxcoverage*(1-seropos) ]
ref[, test70_sp := notest_sp*(PPDmul(refspec, .7, seropos))]
ref[, test70_sn := notest_sn*(PPDmul(refspec, .7, seropos))]
ref[, test90_sp := notest_sp*(PPDmul(refspec, .9, seropos))]
ref[, test90_sn := notest_sn*(PPDmul(refspec, .9, seropos))]

ref.mlt <- melt.data.table(
  ref,
  id.vars = "seropos"
)[, c("scenario", "serostatus") := tstrsplit(variable, "_")]

ref.mlt[,
        ord := 2*c("notest"=2,"test70"=1,"test90"=0)[scenario]+c(sp=0, sn=1)[serostatus]
        ]


lblsize <- 2
lnsize <- 1

cols <- c('seropos', 'value')
ref.mlt2 <- data.table::copy(ref.mlt)[ ,
                                       c("value", "seropos") := list(floor((100)*value),
                                                                     round(100*seropos))
                                       ][(seropos %% 10) == 0]


scale_x <-     scale_x_continuous(labels = scales::percent,
                                  breaks = seq(0,1,by=0.20),
                                  limits = c(0.2, 0.8), expand = expansion(mult = 0.1)) 

scale_palette <- data.frame(color = c(notest = "black", 
                                      test70 = "lightskyblue", 
                                      test90 = "orange"),
                            labels = c("No test",
                                       "70% ≤ Sensitivity < 90%", 
                                       "Sensitivity ≥ 90%")) %>% 
  tibble::rownames_to_column("breaks")


scale_fill <- scale_fill_manual(
  values = setNames(scale_palette$color, scale_palette$breaks),
  labels = scale_palette$labels,
  name   = "Scenario") 

scale_color <- scale_color_manual(
  values = setNames(scale_palette$color, scale_palette$breaks),
  labels = scale_palette$labels,
  name   = "Scenario")

fig1 <- ref.mlt %>%
  dplyr::mutate(Serostatus = forcats::fct_recode(serostatus,
                                                 `Sero-negative` = 'sn',
                                                 `Sero-positive` = 'sp'),
                scenario = factor(scenario)) %>%
  dplyr::bind_rows(.,
                   dplyr::group_by_at(.,
                                      .vars = vars(-c(Serostatus, 
                                                      serostatus, ord, value,
                                                      variable))
                                      ) %>%
                     summarise(value = sum(value)) %>%
                     mutate(Serostatus = "All")) %>%
  ggplot(data = ., aes(y = value, x = seropos)) +
  geom_line(aes(color = scenario)) +
  theme_minimal() +
  ylab("Total immunized, %") + 
  xlab("Percent seropositive") +
  scale_y_continuous(labels = scales::percent, limits = c(0,NA), expand = c(0,0)) +
  scale_x +
  theme(legend.position  = 'bottom',
        legend.direction = 'horizontal',
        legend.spacing.y = unit(0, 'mm')) +
  scale_color +
  facet_grid(cols = vars(Serostatus))#, margins = 'Serostatus')

annotations <- list(data.frame(x = 0.6,
                               y = (PPDmul(refspec, .7, .6) - 1)/2,
                               label = 'prefer tests\nif sensitivity ≥ 70%',
                               scenario = 'test70',
                               hjust = 1),
                    data.frame(x = 0.6,
                               y = (PPDmul(refspec, .7, .8) + PPDmul(refspec, .9, .8) - 2)/3,
                               label = 'prefer tests\nif sensitivity ≥ 90%',
                               scenario = 'test90',
                               hjust = 0.5),
                    data.frame(x = .6,
                               y = .75,
                               label = 'prefer additional\nvaccine doses',
                               scenario = 'notest',
                               hjust = 0)) %>%
  dplyr::bind_rows(.)

ref[, costlim70 := PPDmul(refspec, .7, seropos) - 1 ]
ref[, costlim90 := PPDmul(refspec, .9, seropos) - 1 ]


p.cost <- 
  ggplot(ref) +
  geom_vline(lty = 2, xintercept = 0.4) +
  aes(seropos) +
  #geom_blank(aes(linetype = "notest")) +
  geom_ribbon(aes(ymin=0, ymax=costlim70), 
              fill = 'lightskyblue',
              alpha = 0.5,
              show.legend = FALSE) +
  geom_ribbon(aes(ymin=costlim70, ymax=costlim90), 
              fill = "orange",
              alpha = 0.5,
              show.legend = FALSE) +
  scale_color + 
  scale_fill +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0),
                     "Testing Cost as % of\nVaccination Dose Cost", 
                     labels = function(b) sprintf("%i%%", as.integer(b*100))) +
  scale_x +
  xlab("Percent seropositive") +
  theme_minimal() +
  theme(
    plot.margin = margin(r = 1, unit = "line"),
    legend.position = 'none'
  )

for (i in 1:nrow(annotations)){
  p.cost <- p.cost + 
    geom_label(
      x = annotations[i,"x"],
      y = annotations[i,"y"],
      label = annotations[i,"label"],
      fill = dplyr::filter(scale_palette, breaks == annotations[i,"scenario"]) %>%
        dplyr::pull(color) %>%
        sub(pattern = 'black', replacement = 'white', x = .) ,
      size = lblsize)
}

res.p <- (fig1 + p.cost +
            plot_annotation(tag_levels = "A") +
            plot_layout(guides = "collect", widths = c(3,1))) &
  theme(legend.position = 'bottom')

ggsave(tail(.args, 1), res.p, height = 3, width = 9, dpi = 600)

