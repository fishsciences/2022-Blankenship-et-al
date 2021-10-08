#-------------------------------------------------------#
# Modeling Tidal 5.1 (Oct 2019 experiment)
# This script fits the models compared and presented in Blankenship et al 2022.
###   Lines 48 - 79 fit covariates models
###   Lines 82 - 141 fit models for question of effect of time after dispersion

# M. Johnston, M. Espe
# Tue Sep 28 15:00:21 2021 ------------------------------

# script inputs
#-------------------------------------------------------#
##   - data_clean/tidal5.1_w_hydro.rds (created in analysis/5.1_data_prep.R)
##   - analysis/functions (for convenience fxns)


# script outputs:
#-------------------------------------------------------#

##   - "modK1.0_tidal01.rds" # covariates model with K = 1
##   - "m0.5_tidal01.rds" # covariates model with K = 0.5
##   - "m0.1_tidal01.rds" # covariates model with K = 0.1
##   - "modK1.0_tidal01_disp.rds" # dispersion model with K = 1
##   - "m0.5_tidal01_disp.rds" # dispersion model with K - 0.5
##   - "m0.1_tidal01_disp.rds" # dispersion model with K = 0.1
##   - "uni_mod.rds" # unidirectional model


library(artemis)
source("analysis/functions.R")

## dirs
data_dir = "data_clean"
res_dir = "output" #"~/DropboxCFS/Tidal_results" - old results for comparison

if(!dir.exists(res_dir)) dir.create(res_dir)

## data
tidal = readRDS(file.path(data_dir, "tidal5.1_w_hydro.rds")) # created in analysis/5.1_data_prep.R
plot(Cq ~ SampleTime, tidal)

# Model
#-------------------------------------------------------#
refit_models = FALSE
Ncores = getOption("mc.cores", parallel::detectCores())

#-------------------------------------------------------#
# 5.1 models
if(refit_models){

modK1.0 = eDNA_lmer(Cq ~ scale(K1.0) + stand_dist + Species + Side + (1|FilterID),
                data = tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
               cores = Ncores)

summary(modK1.0)
saveRDS(modK1.0, file.path(res_dir, "modK1.0_tidal01.rds"))


mod0.5 = eDNA_lmer(Cq ~ scale(K0.5) + stand_dist + Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.5)
saveRDS(mod0.5, file.path(res_dir, "m0.5_tidal01.rds"))

mod0.1 = eDNA_lmer(Cq ~ scale(K0.1) + stand_dist + Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.1)
saveRDS(mod0.1, file.path(res_dir, "m0.1_tidal01.rds")) }

#-------------------------------------------------------#
# Dispersion

# Two ideas:
# 1. increasing time since removal
# 2. increasing time since deploy, then decreasing over time
# deployed 2019-10-22 19:23 PDT
# removed  2019-10-23 08:05 PDT

tidal$since_deploy = difftime(tidal$SampleTime, as.POSIXct("2019-10-22 19:23", tz = "Etc/GMT+7"), units = "hours") #PDT
tidal$since_removal = difftime(tidal$SampleTime, as.POSIXct("2019-10-23 08:05", tz = "Etc/GMT+7"), units = "hours") #PDT

# cringe
tidal$deploy_combo = ifelse(tidal$since_removal > 0, 12 - tidal$since_removal, tidal$since_deploy)
plot(tidal$SampleTime, tidal$deploy_combo)

# dispersion predictor
tidal$disp = tidal$since_removal
tidal$disp[tidal$disp < 0] = 0
plot(tidal$SampleTime, tidal$disp)

tidal$detection = ifelse(tidal$Cq < 40, 1, 0)
table(tidal$detection, tidal$since_removal)

ggplot(tidal, aes(x = factor(as.numeric(since_removal)), y = factor(detection))) +
  geom_jitter(aes(color = Species), width = 0.1, size = 1)

# Models

if(refit_models) {
modK1.0_disp = eDNA_lmer(Cq ~ scale(K1.0) +
                              stand_dist +
                              disp +
                              Species +
                              Side +
                              (1|FilterID),
                         data = tidal,
                         std_curve_alpha = tidal$std_alpha,
                         std_curve_beta = tidal$std_beta,
                         seed = 1234,
                         cores = Ncores)

summary(modK1.0_disp)
saveRDS(modK1.0_disp, file.path(res_dir, "modK1.0_tidal01_disp.rds"))

mod0.5_disp = eDNA_lmer(Cq ~ scale(K0.5) + stand_dist + disp + Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.5_disp)
saveRDS(mod0.5_disp, file.path(res_dir, "m0.5_tidal01_disp.rds"))

mod0.1_disp = eDNA_lmer(Cq ~ scale(K0.1) + stand_dist + disp +  Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.1_disp)
saveRDS(mod0.1_disp, file.path(res_dir, "m0.1_tidal01_disp.rds")) }


#-------------------------------------------------------#
# Compare with CVP and delta distance effects
# Unidirectional data: keep the volume columns for all, because volume is different from uni to tidal

library(elaphos)
cvps = rbind(elaphos::cvp01, elaphos::cvp02)

#delta04 = "ds-2019-03-08"
#otherwise = "ds-2018-09-27"

crvs = elaphos::StdCrvKey[elaphos::StdCrvKey$StdCrvID %in% cvps$StdCrvID,
                          c("StdCrvID", "StdCrvAlpha_lnForm", "StdCrvBeta_lnForm")]

cvps$std_alpha = crvs$StdCrvAlpha_lnForm
cvps$std_beta = crvs$StdCrvBeta_lnForm
len(cvps$FilterID)

if(refit_models){

uni_mod = eDNA_lmer(Cq ~ scale(Distance_m) + scale(Volume_mL) + (1|FilterID),
                    data = cvps,
                std_curve_alpha = cvps$std_alpha,
                std_curve_beta = cvps$std_beta,
                seed = 1234,
                cores = Ncores)

summary(uni_mod)
saveRDS(uni_mod, file.path(res_dir, "uni_mod.rds"))}

if(FALSE){ # old models and qaqc code

# pulsey points: singled out because we were trying to isolate potential contamination; email thread confirmed no contamination, but that the livecar drifted close to the +1 point; points were grouped from there forward
off = dplyr::filter(tidal, DateTimePT == "2019-10-23 04:05:00" & Distance_m == -1)
stopifnot(unique(off$Distance_m) == -1)
stopifnot(unique(off$SampleTime == "2019-10-23 04:05:00"))
plot(off$Cq)
tidal = dplyr::anti_join(tidal, off)

tidal$FilterUnq = paste(tidal$FilterID, tidal$Species)

mod = eDNA_lmer(Cq ~ K1.0 + abs(Distance_m) + Side + factor(Biomass_N) + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_crv_alpha,
                std_curve_beta = tidal$std_crv_beta)

mod2 = eDNA_lmer(Cq ~ K0.1 + factor(Biomass_N) + (1|FilterUnq),
                tidal,
                std_curve_alpha = tidal$std_crv_alpha,
                std_curve_beta = tidal$std_crv_beta, verbose = TRUE)

mod3 = eDNA_lm(Cq ~ K0.1 + Species,
                tidal,
                std_curve_alpha = tidal$std_crv_alpha,
                std_curve_beta = tidal$std_crv_beta, verbose = TRUE)

## Combo
if(refit_models){
# modK1.0_deploy_combo = eDNA_lmer(Cq ~ scale(K1.0) + stand_dist +
#                              deploy_combo +
#                              Species + Side + (1|FilterID),
#                          data = tidal,
#                          std_curve_alpha = tidal$std_alpha,
#                          std_curve_beta = tidal$std_beta,
#                          seed = 1234,
#                          cores = Ncores)
#
# summary(modK1.0_deploy_combo)
# saveRDS(modK1.0_deploy_combo, file.path(res_dir, "modK1.0_tidal01_deploy_combo.rds"))

mod0.5_deploy_combo = eDNA_lmer(Cq ~ scale(K0.5) + stand_dist + deploy_combo + Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.5_deploy_combo)
saveRDS(mod0.5_deploy_combo, file.path(res_dir, "m0.5_tidal01_deploy_combo.rds"))

mod0.1_deploy_combo = eDNA_lmer(Cq ~ scale(K0.1) + stand_dist + deploy_combo +  Species + Side + (1|FilterID),
                tidal,
                std_curve_alpha = tidal$std_alpha,
                std_curve_beta = tidal$std_beta,
                seed = 1234,
                cores = Ncores)

summary(mod0.1_deploy_combo)
saveRDS(mod0.1_deploy_combo, file.path(res_dir, "m0.1_tidal01_deploy_combo.rds"))}

}
