# Clean/prep data from tidal experiment 5.1 (Oct 2019)
# Joins 5.1 experimental data with hydro modeling data
# M. Johnston, M. Espe
# Wed Sep 29 15:30:43 2021 ------------------------------

# script inputs:
#-------------------------------------------------------#
##   - 5.3 metadata (matches 5.3_filecompile_with_Pre from S&T repo)
##   - 5.3 qPCR data

# script outputs:
#-------------------------------------------------------#
##   - data_clean/tidal5.3_cleaned.rds
##   - data_clean/tidal5.3_cleaned_with_flow.rds

library(data.table)
#-------------------------------------------------------#

# Sample metadata
#-------------------------------------------------------#
meta = readxl::read_excel("data/sampling_metadata/Expt_5.3_metadata_plate_layouts_original.xlsx",
                          sheet = "metadata_MEJ")

#meta = meta[meta$include == "y", ] # Include is a column added by Von to indicate control/not control

meta = meta[!(meta$sampler %in% c("Pre", "Post")), ] # exclude all pre/post samples for now

keep = c("FilterID", "distance_m", "actual_time", "sampling_interval",
         "date_collected", "volume_ml", "n_ds", "n_anchovy", "n_rbt", "sampler")

meta = data.frame(meta[ , keep])

colnames(meta) <- c("FilterID", "dist", "actual_time",
                    "sampling_interval", "date", "vol",
                    "n_ds", "n_anchovy", "n_rbt", "sampler")

table(meta$n_anchovy)
table(meta$n_rbt)
table(meta$n_ds)

cntrls = grep("fdct_|exct_", meta$FilterID, value = TRUE)

meta = meta[!(meta$FilterID %in% cntrls), ]

meta$date = as.Date(meta$date)

# Raw PCR data
#-------------------------------------------------------#
if(FALSE){ # depends on elaphos package, which is not generally available
cq_files = list.files("data/Expt5.3_xls/",
                      pattern = ".xls", full.names = TRUE)


cqs = lapply(cq_files, FUN = elaphos::read_qpcr, Experiment = "Tidal_5.3")
cqs = do.call(rbind, cqs)
saveRDS(cqs, file = "data/exp5.3_cq.rds")
}
cqs = readRDS("data/exp5.3_cq.rds")
#-------------------------------------------------------#


df1 = cqs[ , c("FilterID", "TechRep", "Cq", "Target")]


sum(meta$FilterID %in% df1$FilterID)
meta$FilterID[!(meta$FilterID %in% df1$FilterID)] # should be 0 now

# filterIDs that do not have metadata/were "post" samples:
rawids = unique(df1$FilterID)
nn = rawids[!(rawids %in% meta$FilterID)] # exct/pre/post

# remove them
df1 = df1[!(df1$FilterID %in% nn), ]

# check
identical(sort(unique(df1$FilterID)), sort(meta$FilterID))

# merged dataframe
df2 = merge(df1, meta, by = "FilterID", all.x = TRUE)

# Make datetime column in PT
df2$time = hms::as_hms(df2$actual_time) # actual sampling time

df2$sampling_datetime = lubridate::force_tz(as.POSIXct(paste(df2$date, df2$time, " ")), "Etc/GMT+8")

# target sampling time:
df2$interval_datetime = lubridate::force_tz(lubridate::ymd_hms(paste(as.character(df2$date),
                                                             df2$sampling_interval,
                                                              " ")), "Etc/GMT+8")

saveRDS(df2, "data_clean/tidal5.3_cleaned.rds")

#-------------------------------------------------------#
# Merge with flow data (from CDEC)

flow = readxl::read_excel("data/sampling_metadata/Expt_5.3_metadata_plate_layouts_original.xlsx",
                          sheet = "flow")

# target sampling time:
flow$sampling_datetime = lubridate::force_tz(as.POSIXct(paste(flow$date, hms::as_hms(flow$time), " ")),
                                           "Etc/GMT+8")
flow_data = as.data.table(data.frame(flow))
dets = as.data.table(df2)

setkeyv(dets, "sampling_datetime")
setkeyv(flow_data, "sampling_datetime")

dets_w_flows = flow_data[ , c("sampling_datetime", "velocity_f.s", "river_stage_ft")][dets, roll = "nearest"]
dets_w_flows= as.data.frame(dets_w_flows)

livecar_start = lubridate::force_tz(as.POSIXct("2020-12-11 09:30:00"), "Etc/GMT+8")
# livecar_end = lubridate::force_tz(as.POSIXct("2020-12-11 17:00:00"), "Etc/GMT+8") # uncomment for post-samples
dets_w_flows$time_elapsed = as.numeric(dets_w_flows$sampling_datetime - livecar_start) # units: minutes

dets_w_flows = dets_w_flows[ , c("sampling_datetime", "velocity_f.s", "river_stage_ft",
                                 "FilterID", "TechRep", "Cq",
                                 "Target", "dist", "date",
                                 "vol", "interval_datetime", "time_elapsed")]

saveRDS(dets_w_flows, "data_clean/tidal5.3_cleaned_with_flow.rds")






