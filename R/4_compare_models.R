# this script compares catalytic model fits

# it depends on:
#  1. the DIC function in ~utils.R
#  2. the stan models fits in ~fits/real_data/ (NOTE - these fits are too large to upload to 
#     github so please regenerate using scripts ~R/2_fit_real_data.R followed by ~R/3_process_fits.R)

# it generates:
#  Table S1 for the manuscript, made up of:
#   a. ~tables/DIC1.csv using core assumptions about sensitivity and specificity of Ab tests
#   b. ~tables/DIC2.csv using alternative assumptions about sensitivity and specificity of Ab tests

source("dependencies.R")

# read in fits

fit1bb_exp_1 <- readRDS("fits/real_data/sens_spec_1/exp/fit1bb.rds")
fit2bb_exp_1 <- readRDS("fits/real_data/sens_spec_1/exp/fit2bb.rds")
fit3bb_exp_1 <- readRDS("fits/real_data/sens_spec_1/exp/fit3bb.rds")
fit4bb_exp_1 <- readRDS("fits/real_data/sens_spec_1/exp/fit4bb.rds")


# DIC across chains

# using core assumptions about sensitivity

DIC_fit1bb_exp_1 <- DIC_across_chains(fit = fit1bb_exp_1, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit2bb_exp_1 <- DIC_across_chains(fit = fit2bb_exp_1, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit3bb_exp_1 <- DIC_across_chains(fit = fit3bb_exp_1, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit4bb_exp_1 <- DIC_across_chains(fit = fit4bb_exp_1, no_chains = 4, total_length = 100000, burnin = 5000)


# generate Table X in manuscript

tab1 <- data.frame(model = c("model 1: seroconversion",
                             "model 2: seroconversion + seroreversion",
                             "model 3: seroconversion + mAbs",
                             "model 4: seroconversion + seroreversion + mAbs"),
                   DIC = round(c(DIC_fit1bb_exp_1$av_DIC$`mean(DIC)`,
                           DIC_fit2bb_exp_1$av_DIC$`mean(DIC)`,
                           DIC_fit3bb_exp_1$av_DIC$`mean(DIC)`,
                           DIC_fit4bb_exp_1$av_DIC$`mean(DIC)`), 1),
                   Dbar = round(c(DIC_fit1bb_exp_1$av_DIC$`mean(D_bar)`,
                            DIC_fit2bb_exp_1$av_DIC$`mean(D_bar)`,
                            DIC_fit3bb_exp_1$av_DIC$`mean(D_bar)`,
                            DIC_fit4bb_exp_1$av_DIC$`mean(D_bar)`), 1),
                   pD = round(c(DIC_fit1bb_exp_1$av_DIC$`mean(pD)`,
                          DIC_fit2bb_exp_1$av_DIC$`mean(pD)`,
                          DIC_fit3bb_exp_1$av_DIC$`mean(pD)`,
                          DIC_fit4bb_exp_1$av_DIC$`mean(pD)`), 1))
write.csv(file = here("tables/DIC1.csv"), tab1)


# assuming lower sensitivity of NTs

# read in fits

fit1bb_exp_2 <- readRDS("fits/real_data/sens_spec_2/exp/fit1bb.rds")
fit2bb_exp_2 <- readRDS("fits/real_data/sens_spec_2/exp/fit2bb.rds")
fit3bb_exp_2 <- readRDS("fits/real_data/sens_spec_2/exp/fit3bb.rds")
fit4bb_exp_2 <- readRDS("fits/real_data/sens_spec_2/exp/fit4bb.rds")

DIC_fit1bb_exp_2 <- DIC_across_chains(fit = fit1bb_exp_2, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit2bb_exp_2 <- DIC_across_chains(fit = fit2bb_exp_2, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit3bb_exp_2 <- DIC_across_chains(fit = fit3bb_exp_2, no_chains = 4, total_length = 100000, burnin = 5000)
DIC_fit4bb_exp_2 <- DIC_across_chains(fit = fit4bb_exp_2, no_chains = 4, total_length = 100000, burnin = 5000)


# generate Table SX in Supplementary material

tab2 <- data.frame(model = c("model 1: seroconversion",
                             "model 2: seroconversion + seroreversion",
                             "model 3: seroconversion + mAbs",
                             "model 4: seroconversion + seroreversion + mAbs"),
                   DIC = round(c(DIC_fit1bb_exp_2$av_DIC$`mean(DIC)`,
                           DIC_fit2bb_exp_2$av_DIC$`mean(DIC)`,
                           DIC_fit3bb_exp_2$av_DIC$`mean(DIC)`,
                           DIC_fit4bb_exp_2$av_DIC$`mean(DIC)`), 1),
                   Dbar = round(c(DIC_fit1bb_exp_2$av_DIC$`mean(D_bar)`,
                            DIC_fit2bb_exp_2$av_DIC$`mean(D_bar)`,
                            DIC_fit3bb_exp_2$av_DIC$`mean(D_bar)`,
                            DIC_fit4bb_exp_2$av_DIC$`mean(D_bar)`), 1),
                   pD = round(c(DIC_fit1bb_exp_2$av_DIC$`mean(pD)`,
                          DIC_fit2bb_exp_2$av_DIC$`mean(pD)`,
                          DIC_fit3bb_exp_2$av_DIC$`mean(pD)`,
                          DIC_fit4bb_exp_2$av_DIC$`mean(pD)`), 1))
write.csv(file = here("tables/DIC2.csv"), tab2)

