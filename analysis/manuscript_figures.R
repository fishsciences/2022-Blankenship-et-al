#-------------------------------------------------------#
# Manuscript plots
# M. Johnston
# Thu Sep 30 15:25:49 2021 ------------------------------
#-------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(fishpals)

#-------------------------------------------------------#
# Covariate effects  plot (Figure 2)

allfx = read.csv("output/allfx.csv")
str(allfx)
colnames(allfx) = c("parm",
  "Mean",
                    "lwr",
                    "Median",
                    "upr",
                    "Model")

alphas = grep(pattern = "Intercept", allfx$parm)
sigmas = grep(pattern = "sigma", allfx$parm)

i = allfx[-c(alphas,sigmas), ]
i = i[-2, ]

i$parm2 = c("Distance (m)",
            # "Eddy diffusivity (K=0.5)",
            # "Distance (m)",
            # "Target spp. (Delta Smelt)",
            # "Side (upstream)",
            #"Tracer concentration (K=1.0)",
            "Distance (m)",
            "Time elapsed since cage removal (hours)",
            "Target spp. (Delta Smelt)",
            "Side (north)",
            # "Eddy diffusivity (K=0.5)",
            # "Distance (m)",
            # "Time deployed*elapsed (hrs)",
            # "Target spp. (Delta Smelt)",
            # "Side (upstream)",
            "Distance (m)",
            "Target spp. (Delta Smelt)",
            "Target spp. (steelhead trout)",
            "Tidal velocity (ft/s)",
            "River stage (ft)",
            "Side (north)")

i = arrange(i, Model, parm2)

pos = position_dodge(width = 0.75)

i %>%
ggplot(aes(y = Median, x = parm2)) +
  geom_point(aes(color= Model,
                 shape = Model,
                 group = Model),
             position = pos,
             size = 3,
             alpha = 0.75) +
  geom_linerange(aes(ymin = lwr,
                   ymax = upr,
                   # xmin = reorder(parm2, Model),
                   # xmax = reorder(parm2, Model),
                   group = Model,
                   color = Model),
               size = 1,
               alpha = 0.75,
               position = pos) +
    coord_flip() +
  geom_hline(aes(yintercept = 0),
             lty = 2,
             size = 0.5) +
  theme_report() +
  fishpals::scale_color_fishpals() +
  labs(y = "Median and 95% Credible Interval \nof estimated effect on ln[eDNA]",
       x = NULL) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal")

ggsave("figs/fx_plot.jpg", height = 6.5, width = 7.5)
