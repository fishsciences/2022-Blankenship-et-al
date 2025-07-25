# Prep data for analysis of the effect of distance in tidal vs. unidirectional systems ("Q2" in manuscript)
# M. Espe
# 28 Apr 2021
#-------------------------------------------------------#

library(artemis)
library(elaphos)

# Data
all_uni = readRDS("data/all_unidirectional_expts.rds")
cvps = rbind(cvp01, cvp02)

d5.3 = readRDS("data_clean/tidal5.3_cleaned_with_flow.rds")
d5.1 = readRDS("data_clean/tidal5.1_w_hydro.rds")


#----------------------------------------#
# Munge tidal

# Add biomass to 5.3
d5.3$Biomass_N = as.numeric(as.character(factor(d5.3$Target, labels = c(25,25,15))))

# Add StdCrv to 5.3 (taken from 5.3_analysis.R
stdcrvs_tmp = data.frame(StdCrvID = c("an-2019-10-29", "rbt-2019-07-16",
                                      "ds-2020-04-06"),
                         Species = c("Anchovy",
                                     "Rainbow Trout",
                                     "Delta Smelt"),
                         Sp_abbv = c("Anch", "STH", "DS"))
i = match(d5.3$Target, stdcrvs_tmp$Sp_abbv)
d5.3$StdCrvID = stdcrvs_tmp$StdCrvID[i]

# need to rename some columns
cnms5.3 = colnames(d5.3)

cnms5.3[cnms5.3 == "dist"] = "Distance_m"
cnms5.3[cnms5.3 == "vol"] = "Volume_mL"
cnms5.3[cnms5.3 == "date"] = "Date"
colnames(d5.3) = cnms5.3

# make a "Target" column in 5.1
d5.1$Target = as.character(factor(d5.1$Species, labels = c("Anch", "DS")))
# and Species in 5.3
d5.3$Species = as.character(factor(d5.3$Target, labels = c("Anchovy", "Delta Smelt", "Rainbow Trout")))

#Combine
tar = intersect(colnames(d5.1), colnames(d5.3))
tar # Just checking this includes everything we want

all_tidal = rbind(d5.3[,tar], d5.1[,tar])
all_tidal$Experiment = rep(c("5.3", "5.1"), c(nrow(d5.3), nrow(d5.1)))

#----------------------------------------#
# munge uni
tar = intersect(colnames(cvps), colnames(all_tidal))

all_data = rbind(cvps[,tar], all_tidal[,tar])

#----------------------------------------#
# get std curves

i = match(all_data$StdCrvID, StdCrvKey$StdCrvID)
all_data$std_crv_a = StdCrvKey$StdCrvAlpha_lnForm[i]
all_data$std_crv_b = StdCrvKey$StdCrvBeta_lnForm[i]

stopifnot(all(complete.cases(all_data))) # check for NA

saveRDS(all_data, "data_clean/all_data_combined.rds")

#-------------------------------------------------------#
