# Clean/prep data from tidal experiment 5.1 (Oct 2019)
# Joins 5.1 experimental data with hydro modeling data
# M. Johnston, M. Espe
# Tue Sep 28 14:47:51 2021 ------------------------------

# script inputs:
#-------------------------------------------------------#
##   - tidal01 (elaphos)
##   - hydro data (data/hydro)

# script outputs:
#-------------------------------------------------------#
##   - data_clean/tidal5.1_w_hydro.rds; this is the input to analysis/5.1_analysis.R

library(elaphos)
source("analysis/functions.R") # data cleaning and hydro joining functions
res_dir = "data_clean" # for storing cleaned data (output of this script)

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

if(!dir.exists(res_dir)) dir.create(res_dir)

saveRDS(tidal5.1, file.path(res_dir, "tidal5.1_w_hydro.rds"))
