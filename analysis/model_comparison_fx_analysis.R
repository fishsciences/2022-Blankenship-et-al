# M. Johnston
# model comparison for unidirectional and tidal 5.1 and 5.3
# and getting effects estimates for use in tidal draft manuscript
# Wed Jun 16 14:14:10 2021 ------------------------------

library(loo)
res_dir = "output" #"~/DropboxCFS/Tidal_results" - old results for comparison
#res_dir = "~/DropboxCFS/Tidal_results/"

# Load fit model objects
#-------------------------------------------------------#
# unidirectional
um = readRDS(file.path(res_dir, "uni_mod.rds"))

# Tidal Expt 5.1 (October 2019)
m5.1 = readRDS(file.path(res_dir, "mod_tidal01.rds"))

m5.1_disp = readRDS(file.path(res_dir, "mod_tidal01_disp.rds"))

# Tidal Expt 5.3 (Dec 2020)
m5.3_old = readRDS(file.path(res_dir, "5.3_model_fit2.rds"))

m5.3 = readRDS(file.path(res_dir, "5.3_model_sidefit.rds")) # this was the better fit according to loo_compare on # Wed Jun 30 2021 ------------------------------

if(FALSE){ # K parameters dropped
# Model comparison, Tidal Expt 5.1 (October 2019)
# rough guide is the elpd_diff has to be more than 2x the se_diff for one model to be preferred
# first compare 2 sets separately (to see if one K value is better than another), then compare all six together to see if there is a best overall model fit
#-------------------------------------------------------#
lhi = loo::loo(m5.1_hi)
lmed = loo::loo(m5.1_med)
llo = loo::loo(m5.1_lo)

m5.1_comp = loo::loo_compare(lhi, lmed, llo)
m5.1_comp = data.frame(
                      cbind(Model = c("K = 1.0 (modK1.0_tidal01.rds)",
                                       "K = 0.5 (m0.5_tidal01.rds)",
                                       "K = 0.1 (m0.1_tidal01.rds)"),
                             m5.1_comp)
                       )

m5.1_comp = apply(m5.1_comp[ , c(2:9)], 2, function(x) round(as.numeric(x), 2))

write.csv("output/m5.1_comp.csv", row.names = TRUE)
#-------------------------------------------------------#
# dispersal, expt 5.1

K1_disp = readRDS(file.path(res_dir, "modK1.0_tidal01_disp.rds"))
lhi_disp = loo::loo(K1_disp)

K05_disp = readRDS(file.path(res_dir, "m0.5_tidal01_disp.rds"))
lmed_disp = loo::loo(K05_disp)

Kpt1_disp = readRDS(file.path(res_dir, "m0.1_tidal01_disp.rds"))
llo_disp = loo::loo(Kpt1_disp)

disp = loo::loo_compare(lhi_disp, lmed_disp, llo_disp) #

m5.1_disp_comp = data.frame(
                      cbind(Model = c("K = 1.0 (modK1.0_tidal01_disp.rds)",
                                       "K = 0.5 (m0.5_tidal01_disp.rds)",
                                       "K = 0.1 (m0.1_tidal01_disp.rds)"),
                             disp)
                       )

  write.csv(m5.1_disp_comp, "output/m5.1_disp_comp.csv", row.names = TRUE)
  }
#-------------------------------------------------------#
# full comparison

m5.1_allmodels = loo::loo_compare(loo(m5.1), loo(m5.1_disp))
write.csv(m5.1_allmodels, "output/m5.1_allmodels.csv", row.names = TRUE)


# Model comparison, Expt 5.3 (Dec 2020)
#-------------------------------------------------------#
lc = loo::loo_compare(artemis:::loo.eDNA_model(m5.3),
                      artemis:::loo.eDNA_model(m5.3_old))

write.csv(print(data.frame(lc), digits = 2),
          "output/m5.3_comp.csv", row.names = TRUE)

lc # m5.3 (5.3_model_sidefit.rds) is the better fit as of
# Wed Jun 30 13:45:12 2021 ------------------------------
# estimates in the paper should use that one



# Write effects to table
#-------------------------------------------------------#

allfx = rbind(cbind(summary(um), Model = "Unidirectional"),
             # cbind(summary(m5.1_med), Model = "October 2019"),
              cbind(summary(m5.1_disp), Model = "October 2019"),
              # cbind(summary(k0.5_combo), Model = "Dispersal*Time"),
              cbind(summary(m5.3), Model = "December 2020")
              )


allfx[ , 1:4] = round(allfx[, 1:4], digits =3)

write.csv(allfx, "output/allfx.csv", row.names = TRUE)


# Q2 analysis effects sizes for manuscript
#-------------------------------------------------------#
q2 = readRDS(file.path(res_dir, "Q2_final_model_fit.rds"))
summary(q2)

write.csv(summary(q2), "output/q2_fx.csv")
