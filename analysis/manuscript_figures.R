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

i |>
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


#------------------------------------------------------------------------------#
# Tidal figures
# M. Johnston edited by M. Espe
# Mon Apr 19 15:20:21 2021 ------------------------------

library(ggplot2)
library(fishpals)
library(dplyr)
pull_CDEC_data = TRUE

#-------------------------------------------------------#
# 5.1
#-------------------------------------------------------#
## data
tidal = readRDS("data_clean/tidal5.1_w_hydro.rds") # created in 5.1_data_prep.R

if(pull_CDEC_data | !file.exists("data_clean/t5.1_cdec_data.rds")){
t5.1 = CDECRetrieve::cdec_query("SGG", 1, dur_code = "E", # river stage in ft
                                start_date = "2019-10-22",
                                end_date = "2019-10-24")

t5.1 = t5.1[t5.1$datetime >= lubridate::ymd_hms("2019-10-22 18:00:00", tz = "America/Los_Angeles") &
              t5.1$datetime <= lubridate::ymd_hms("2019-10-23 20:00:00", tz = "America/Los_Angeles"), ]

saveRDS(t5.1, "data_clean/t5.1_cdec_data.rds")

}else{
  t5.1 = readRDS("data_clean/t5.1_cdec_data.rds")
}

# munge data
#-------------------------------------------------------#
tidal$abs_dist_m = abs(tidal$Distance_m)
mu = mean(tidal$abs_dist_m)
sd = sd(tidal$abs_dist_m)
tidal$stand_dist = scale(tidal$abs_dist_m, mu, sd) # scale distances

# stdcurvs - from elaphos::StdCrvKey.  
stdcrvs = structure(list(
  StdCrvID = c("an-2019-10-29", "ds-2019-10-29"), 
StdCrvAlpha_lnForm = c(18.441, 20.583), 
StdCrvBeta_lnForm = c(-1.545, -1.486)), 
row.names = c(3L, 4L), class = "data.frame")

stdcrvs$Sp_abbv = c("Anchovy", "Delta Smelt")
i = match(tidal$Species, stdcrvs$Sp_abbv)

tidal$std_alpha = stdcrvs$StdCrvAlpha_lnForm[i]
tidal$std_beta = stdcrvs$StdCrvBeta_lnForm[i]

tidal$Distance_m[tidal$Distance_m == -1] <- 1

ts_ix = tidal |> 
  group_by(Distance_m, SampleTime) |> 
  summarize(n_reps = n()) |> 
  ungroup()

ds_plot = tidal |> 
 # filter(Species == "Delta Smelt") |> 
  select(SampleTime, Cq, `K0.5`, Distance_m, Species) |> 
  group_by(Distance_m, SampleTime, Species) |> 
  summarize(npos = sum(Cq < 40),
            k = unique(`K0.5`)) |> 
  ungroup() |> 
  left_join(ts_ix) |> 
  mutate(prop_pos = npos/n_reps) |> 
  mutate(shape = ifelse(npos == 0, "neg", "pos"))

ds_plot = ds_plot |> 
group_by(SampleTime) |> 
  mutate(k= max(k)) |> 
  ungroup() 

ds_plot$k_scaled = (sd(ds_plot$Distance_m)*(scale(ds_plot$k)) + mean(ds_plot$Distance_m)) 
t5.1$stg_scaled = -1*(sd(ds_plot$Distance_m)*(scale(t5.1$parameter_value)) + mean(ds_plot$Distance_m)) 

## Stage
#-------------------------------------------------------#
ggplot(ds_plot, aes(y = Distance_m, x = SampleTime)) +
  geom_point(
    ## geom_raster(
    aes(
      shape = prop_pos > 0,
      ## fill = prop_pos,
      size = prop_pos
    )
  ) +
  geom_line(data = t5.1, 
            size = 0.65,
            lty = 3,
             aes(
                 x = datetime, y = stg_scaled),
             alpha = 0.85) +
  ##   scale_fill_gradient(low = "white", high = "gray20",
  ##                       aesthetics = "fill",
  ##                       ## limits = c(0.0, 0.61),
  ##                       breaks = c(0.0, 0.25, 0.50, 0.6),
                        
  ##                       guide = guide_legend(
  ##                           title = "Proportion Positive",
  ##                           ## override.aes = list(shape = 21, 
  ##                                               ## size = 4),
  ##                           ## breaks = c(0.01, 0.25, 0.50, 0.6)
  ##                         )
  ##                       ) + 
  ## guides(shape = "none", size = "none", color = "none") +
  scale_shape_manual(values = c(1, 20),
                     guide = guide_legend(
                     title = "Positive sample")) + 
  scale_size(breaks = c(0.01, 0.25, 0.5),
             guide = guide_legend(title = "Proportion positive",
                                  breaks = c(0.01, 0.25, 0.5))) + 
  scale_y_reverse(breaks = c(600, 400, 200, 100, 50, 1, -50, -100, -200, -400, -600)) +
  scale_x_datetime(date_breaks = "3 hours",
                   date_labels = "%H:%M") +
  geom_vline(aes(xintercept = lubridate::ymd_hms("2019-10-23 08:05:00", tz = "America/Los_Angeles")),
             size = 0.4,
             lty = 2,
             color = "gray35") +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  labs(x = "Sampling time",
       y = "Distance from livecar (m)") +
  facet_wrap(~Species, nrow = 2)

ggsave("figs/5.1_ds_stg_pts.png", height = 11, width = 8.5, bg = "white")

range(t5.1$parameter_value)

## K (eddy diffusivity)
#-------------------------------------------------------#
tidal$conc = artemis::cq_to_lnconc(tidal$Cq, tidal$std_alpha, tidal$std_beta)

tidal |> 
  #filter(Cq < 40) |> 
  mutate(livecar_status = ifelse(SampleTime >= lubridate::ymd_hms("2019-10-23 08:05:00", tz = "America/Los_Angeles"), "out", "in")) |> 

ggplot(aes(x = livecar_status, y = conc)) +
  geom_jitter(aes(color = Species),
              position = position_jitterdodge(dodge.width = 0.3,
                                              jitter.width = 0.1),
              alpha = 0.2) +
  geom_boxplot(aes(color = Species), 
               outlier.shape = NA,
               alpha = 0.5,
               width = 0.3) +
  theme_report()


#-------------------------------------------------------#
# 5.3
#-------------------------------------------------------#

d = readRDS("data_clean/tidal5.3_cleaned_with_flow.rds")

#-------------------------------------------------------#
# what we want: proportion positive tech reps for each sampling event by target and distance.
# group by: Target, distance, interval; total the n tech reps:

d$Species = factor(d$Target, levels = c("STH", "DS", "Anch"),
                   labels = c("Steelhead", "Delta Smelt", "Anchovy"))
tr_ix = d |> 
  group_by(dist, interval_datetime) |> 
  summarize(n_reps = n()) |> 
  ungroup()

sth_plot = d |> 
   # filter(Target == "STH" ) |> #, Cq < 40) |> 
  select(Target, interval_datetime, Cq, river_stage_ft, velocity_f.s, dist, Species) |> 
  group_by(dist, interval_datetime, Species) |> 
  summarize(npos = sum(Cq < 40),
            river_stage_ft = unique(river_stage_ft),
            velocity_f.s = unique(velocity_f.s)) |> 
  ungroup() |> 
  left_join(tr_ix) |> 
  mutate(prop_pos = npos/n_reps) |> 
  mutate(shape = ifelse(npos == 0, "neg", "pos"))

# standardize vel & river stage across distances:

sth_plot = sth_plot |> 
  group_by(interval_datetime) |> 
  mutate(river_stage_ft = max(river_stage_ft),
         velocity_f.s = max(velocity_f.s)) |> 
  ungroup() 

# if you want to scale one to the other, its easiest done using the max values of either
# e.g. max(y) = ?? * max(x)

flow = CDECRetrieve::cdec_query("SGG", sensor_num = 1, dur_code = "E",
                                start_date = "2020-12-11",
                                end_date = "2020-12-12")

flow = flow[flow$datetime > "2020-12-11 05:15" & flow$datetime < "2020-12-11 21:45", ]

sth_plot$stage_max_scaled = (sd(sth_plot$dist)*(scale(sth_plot$river_stage_ft)) + mean(sth_plot$dist)) 
sth_plot$vel_max_scaled = (sd(sth_plot$dist)*(scale(sth_plot$velocity_f.s)) + mean(sth_plot$dist)) 

#-------------------------------------------------------#
sth_plot$dist[sth_plot$dist == -600] <- -300
sth_plot$dist[sth_plot$dist == 600] <- 300

scaled_dists = c(-250, -100, -50, -25, 25, 50, 100, 250)

flow$stg_scaled = (sd(scaled_dists)*(scale(flow$parameter_value)) + mean(scaled_dists)) 

# scaled river stage
ggplot(sth_plot, 
       aes(y = dist, x = interval_datetime)) +
  geom_point(
    
  ## geom_raster(
    aes(
            ## fill = prop_pos ,
            size = prop_pos,
            shape = prop_pos > 0,
        )
    ) +
   geom_line(data = flow, 
            size = 0.65,
            lty = 3,
             aes(
                 x = datetime, y = stg_scaled),
            alpha = 0.85) +
  ##   scale_fill_gradient(low = "white", high = "gray20",
  ##                       aesthetics = "fill",
  ##                       ## limits = c(0.01, 0.61),
  ##                       ## breaks = c(0.01, 0.25, 0.50, 0.6),
                        
  ##                       guide = guide_legend(
  ##                           title = "Proportion Positive",
  ##                           override.aes = list(shape = 21, 
  ##                                               size = 4),
  ##                           breaks = c(0.01, 0.25, 0.50, 0.6))
  ##                       ) +
  scale_shape_manual(values = c(1, 20),
                     guide = guide_legend(
                       title = "Positive sample")) + 
  scale_size(breaks = c(0.01, 0.1, 0.2),
             guide = guide_legend(title = "Proportion positive",
                                  breaks = c(0.01, 0.25, 0.5, 0.75))) + 
  ## guides(shape = "none", size = "none", color = "none") +
  scale_y_continuous(breaks = c(-300, -100, -50, -25, 0, 25, 50, 100, 300),
                     labels = c("-600",
                                "-100",
                                "-50",
                                "-25",
                                "0",
                                "25",
                                "50",
                                "100",
                                "600")) +
  scale_x_datetime(date_breaks = "3 hours",
                   date_labels = "%H:%M") +
  theme_pub() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13)
  ) +
  labs(x = "Sampling time",
       y = "Distance from livecar (m)") +
  geom_hline(aes(yintercept = 0),
             size = 0.4) +
  facet_wrap(~Species, nrow = 3)

ggsave("figs/5.3_proppos_stage_pts.png", height = 11, width = 8.5, bg = "white")


#-------------------------------------------------------#
# scaled vel

ggplot(sth_plot[sth_plot$prop_pos > 0, ], aes(y = dist, x = interval_datetime)) +
  geom_point(
    shape = 21,
    aes(
        fill = prop_pos ,
        size = 4
    )
  ) +
    scale_fill_gradient(low = "white", high = "gray20",
                        aesthetics = "fill",
                        limits = c(0.01, 0.61),
                        breaks = c(0.01, 0.25, 0.50, 0.6),
                        
                        guide = guide_legend(
                            title = "Proportion Positive",
                            override.aes = list(shape = 21, 
                                                size = 4),
                            breaks = c(0.01, 0.25, 0.50, 0.6))
                        ) + 
  guides(shape = "none", size = "none", color = "none") +
  scale_y_continuous(breaks = unique(d$dist)) +
  theme_report(inner_border = FALSE) +
  scale_x_datetime(date_breaks = "1.5 hours",
                   date_labels = "%H:%M") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text = element_text(size = 8)
  ) +
  labs(x = "Sampling time",
       y = "Distance from livecar (m)") +
  # facet_wrap(~interval_datetime ,
  #            ncol = 1
  #            ) +
  geom_smooth(data = sth_plot[sth_plot$prop_pos > 0, ],
    aes( x = interval_datetime, y = vel_max_scaled),
    size = 0.35,
    color = "black",
    se = TRUE
  ) +
  geom_point(data = sth_plot[sth_plot$prop_pos > 0, ], 
             aes(
                 x = interval_datetime, y = vel_max_scaled),
             shape = 17,
             alpha = 0.5)

ggsave("figs/5.3_proppos_vel.png", height = 10, width = 8.5, bg = "white")


sth_plot[sth_plot$dist == -600 & sth_plot$interval_datetime == max(sth_plot$interval_datetime), ]


tide = d |> 
  select(interval_datetime, dist, river_stage_ft, velocity_f.s, sampling_datetime)

ggplot(tide, aes(x = interval_datetime, y = river_stage_ft)) +
  geom_line() +
  geom_point() +
  scale_x_datetime(date_breaks = "1.5 hours",
                   date_labels = "%H:%M")

#-------------------------------------------------------#
# number of hours past which 0 detections:
str(tidal)
table(tidal$SampleTime)
