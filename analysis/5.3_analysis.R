#-------------------------------------------------------#
# M. Johnston
# Expt 5.3 analysis
# Goal: Model the effects of Distance, target, time since deployment,  and tidal alignment (velocity),  and river stage on ln[eDNA]
#-------------------------------------------------------#
# last edited # Thu Sep 30 09:38:49 2021 ------------------------------
#-------------------------------------------------------#

# script inputs:
#-------------------------------------------------------#
##   - data_clean/tidal5.3_cleaned_with_flow.rds
##   - data_clean/tidal5.1_w_hydro.rds   ## to scale the abs dist column of 5.3

# script outputs:
#-------------------------------------------------------#
##   - results/5.3_model_fit2.rds
##   - results/5.3_model_sidefit.rds

library(artemis)

d = readRDS("data_clean/tidal5.3_cleaned_with_flow.rds")
res_dir = "results" #"~/DropboxCFS/Tidal_results" # where the model fits are saved

# stdcurvs - from elaphos::StdCrvKey.  Apply to 5.3 only
stdcrvs = structure(list(StdCrvID = c("an-2019-10-29", "rbt-2019-07-16",
                                      "ds-2020-04-06"),
                         Species = c("Anchovy",
                                     "Rainbow Trout",
                                     "Delta Smelt"),
                         StdCrvAlpha_lnForm =
                             c(18.441, 21.523, 18.625),
                         StdCrvBeta_lnForm = c(-1.545,
                                               -1.4985, -1.542)),
                    row.names = 32:34, class = "data.frame")

stdcrvs$Sp_abbv = c("Anch", "STH", "DS")
i = match(d$Target, stdcrvs$Sp_abbv)

d$std_alpha = stdcrvs$StdCrvAlpha_lnForm[i]
d$std_beta = stdcrvs$StdCrvBeta_lnForm[i]


# other covariates
d$Side = ifelse(d$dist > 0, "downstream", "upstream")
colSums(is.na(d)) # should be no NAs now


# add abs(distance) columns to both 5.1 and 5.3
# tidal01 is not modeled here, but the abs dist  column of 5.3 is scaled to the distances in 5.1 for better comparison of the distance effects across models
#-------------------------------------------------------#
tidal5.1 = readRDS("data_clean/tidal5.1_w_hydro.rds") # created in analysis/5.1_data_prep.R
tidal5.1$abs_dist_m = abs(tidal5.1$Distance_m)
mu = mean(tidal5.1$abs_dist_m)
sd = sd(tidal5.1$abs_dist_m)

d$abs_dist_m = abs(d$dist)

#-------------------------------------------------------#
# Modeling 5.3
#-------------------------------------------------------#
Ncores = getOption("mc.cores", parallel::detectCores())
refit_models = FALSE

if(refit_models){
m1_old = eDNA_lmer(Cq ~ scale(abs_dist_m, mu, sd) + # mu and sd are scaled to 5.1
                           Target +
                           scale(velocity_f.s) +
                           scale(river_stage_ft) +
                       #    Side +
                           (1|FilterID),
                       data = d,
                       seed = 20,
                       std_curve_alpha = d$std_alpha,
                       std_curve_beta = d$std_beta,
                       cores = Ncores)

summary(m1_old)
saveRDS(m1_old, paste0(res_dir, "/5.3_model_fit2.rds")) }


if(refit_models) {
m1_formula = eDNA_lmer(Cq ~ scale(abs_dist_m, mu, sd) +
                           Target +
                           scale(velocity_f.s) +
                           scale(river_stage_ft) +
                           Side +
                           (1|FilterID),
                       data = d,
                       seed = 20,
                       std_curve_alpha = d$std_alpha,
                       std_curve_beta = d$std_beta,
                       cores = Ncores)

saveRDS(m1_formula, paste0(res_dir, "/5.3_model_sidefit.rds"))}



#-------------------------------------------------------#
# archived models as of
# Wed Jun 16 14:10:16 2021 ------------------------------
#-------------------------------------------------------#
if(FALSE){

  # biomass (not modeled for manuscript b/c confounded by Target):
d$weight_g = NA
d$weight_g[d$Target == "Anch"] <- 125.8
d$weight_g[d$Target == "DS"] <- 132.2
d$weight_g[d$Target == "STH"] <- 901.7


m2 = eDNA_lmer(Cq ~ scale(dist) +
                           factor(weight_g) +
                           scale(velocity_f.s) +
                           scale(river_stage_ft) +
                           (1|FilterID),
                       data = d,
                       std_curve_alpha = d$std_alpha,
                       std_curve_beta = d$std_beta,
                       cores = Ncore)

# Draft model - incomplete; not sure how to do multiple standard curves?
d2 = d[d$Cq < 40, ]

m1_pos = eDNA_lmer(Cq ~ scale(dist) +
                           Target +
                           scale(velocity_f.s) +
                           scale(river_stage_ft) +
                           (1|FilterID),
                       data = d2,
                       std_curve_alpha = d2$std_alpha,
                       std_curve_beta = d2$std_beta,
                       cores = Ncores)
}
