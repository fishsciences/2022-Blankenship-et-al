# Clean/prep data from tidal experiment 5.1 (Oct 2019)
# Joins 5.1 experimental data with hydro modeling data
# M. Johnston, M. Espe
# Tue Sep 28 14:47:51 2021 ------------------------------

# script inputs:
#-------------------------------------------------------#
##   - hydro data (data/hydro)

# script outputs:
#-------------------------------------------------------#
##   - data_clean/tidal5.1_w_hydro.rds; this is the input to analysis/5.1_analysis.R

source("analysis/functions.R") # data cleaning and hydro joining functions
res_dir = "data_clean" # for storing cleaned data (output of this script)

tidal01 = readRDS("data/tidal01.rds") # formally contained in elaphos

## Combine with hydro model output
fs = list.files("data/hydro", full.names = TRUE)

d = lapply(fs, read.csv)

d = lapply(d, function(x){
    x$Time = as.POSIXct(x$Time, format = "%d-%b-%Y %H:%M:%S", tz = "America/Los_Angeles")
    x$Car = rowMeans(x[,c("NEG1", "LC", "POS1")])
    x
})


hydros = lapply(d, function(x) hydro_vals(tidal01, x))
names(hydros) = gsub("livecar_simulationv3run[0-9](.*)\\.csv$", "\\1", basename(fs))
hydros = as.data.frame(hydros)
tidal5.1 = cbind(tidal01, hydros)

# munge data
#-------------------------------------------------------#
tidal5.1$abs_dist_m = abs(tidal5.1$Distance_m)
mu = mean(tidal5.1$abs_dist_m)
sd = sd(tidal5.1$abs_dist_m)
tidal5.1$stand_dist = scale(tidal5.1$abs_dist_m, mu, sd) # scale distances

# stdcurvs - from elaphos::StdCrvKey.
stdcrvs = structure(list(
  StdCrvID = c("an-2019-10-29", "ds-2019-10-29"),
StdCrvAlpha_lnForm = c(18.441, 20.583),
StdCrvBeta_lnForm = c(-1.545, -1.486)),
row.names = c(3L, 4L), class = "data.frame")

stdcrvs$Sp_abbv = c("Anchovy", "Delta Smelt")
i = match(tidal5.1$Species, stdcrvs$Sp_abbv)

tidal5.1$std_alpha = stdcrvs$StdCrvAlpha_lnForm[i]
tidal5.1$std_beta = stdcrvs$StdCrvBeta_lnForm[i]

if(!dir.exists(res_dir)) dir.create(res_dir)

saveRDS(tidal5.1, file.path(res_dir, "tidal5.1_w_hydro.rds"))
