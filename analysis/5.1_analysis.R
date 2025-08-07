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
refit_models = TRUE
Ncores = getOption("mc.cores", parallel::detectCores())

#-------------------------------------------------------#
# 5.1 models
# All K params were dropped
if(refit_models){

  mod = eDNA_lmer(Cq ~ #scale(K1.0) +
                    stand_dist + Species + Side + (1|FilterID),
                  data = tidal,
                  std_curve_alpha = tidal$std_alpha,
                  std_curve_beta = tidal$std_beta,
                  seed = 1234,
                  parallel_chains = Ncores)

summary(mod)
saveRDS(mod, file.path(res_dir, "mod_tidal01.rds"))
}

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
mod_disp = eDNA_lmer(Cq ~ #scale(K1.0) +
                              stand_dist +
                              disp +
                              Species +
                              Side +
                              (1|FilterID),
                         data = tidal,
                         std_curve_alpha = tidal$std_alpha,
                         std_curve_beta = tidal$std_beta,
                         seed = 1234,
                         parallel_chains = Ncores)

summary(mod_disp)
saveRDS(mod_disp, file.path(res_dir, "mod_tidal01_disp.rds"))
}

#-------------------------------------------------------#
# Compare with CVP and delta distance effects
# Unidirectional data: keep the volume columns for all, because volume is different from uni to tidal

cvps = readRDS("data/cvp_unidirectional.rds")

crvs = readRDS("data/unidirectional_std_crvs.rds")

cvps$std_alpha = crvs$StdCrvAlpha_lnForm
cvps$std_beta = as.numeric(crvs$StdCrvBeta_lnForm)
len(cvps$FilterID)

if(refit_models){

uni_mod = eDNA_lmer(Cq ~ scale(Distance_m) + scale(Volume_mL) + (1|FilterID),
                    data = cvps,
                std_curve_alpha = cvps$std_alpha,
                std_curve_beta = cvps$std_beta,
                seed = 1234,
                parallel_chains = Ncores)

summary(uni_mod)
  saveRDS(uni_mod, file.path(res_dir, "uni_mod.rds"))
}

