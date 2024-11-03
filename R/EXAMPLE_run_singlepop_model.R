
# this script is an example script that can be adapted by the user to run their
# own simulations using the stochastic single-population model of MERS-CoV transmission 
# and vaccination in dromedary camels.

# Run this if you want to assume a well-mixed homogeneous camel population. If you want to model population
# structures i.e. where camels have more contact with those in the same herd/area and less contact with those 
# in neighbouring herds/areas you can use "EXAMPLE_run_metapop_model.R" (this takes much longer to run)

# This example depends on the model "
#     - dynamic-odin-models/single_patch_model_vax.R"
#     - the estimated rate of waning of maternal antibodies from the catalytic modelling
#       "fits/processed_real/exp/sens_spec_1/fit4bb.rds" but you can replace this with any value
#       of your choice

# This script generates a data frame "out" with 1 row per day and one column per model output per stochastic run
# so for the default example for 30 years with all 814 outputs and 100 stochastic model runs the dimensions of "out" are 10801x81400

#######################
## load dependencies ##
#######################

# first lets load the necessary dependencies
library(here) # this helps make file paths specific to your directory so the same code can be used by different people
source(here("dependencies.R"))

# load the single population transmission model with vaccination
sp_model <- odin::odin("dynamic-odin-models/single_patch_model_vax.R", 
                            verbose = FALSE, skip_cache = TRUE)

#########################
## specify time period ##
#########################

# input the time period that you wish to simulate transmission in (in days)
# NOTE - I tend to use 360 as an approximation for a year as this makes it easier to align ageing of camels and monthly interventions
time_period <- 10*360 # e.g ten years
t <- seq(0:(time_period + 20*360)) # add this time_period to the set up period of 20 years (10 years with background FoI on, 10 years to reach equilibrium)

########################
## specify parameters ##
########################

##############
# demography #
##############

N_age <- 49 # number of age classes - DO NOT CHANGE

N_0 <- 1000000 # set the total size of initial population

# mean birth rate (per camel per day) 
alpha <- 0.000565 # default 0.000565: 90% female * 50% reproductive age * 45.2% fecundity 
## based on doi: 10.9755/ejfa.v25i4.15491 and Abbas et al. 2000 Revue Élev. Méd. vét. Pays trop., 2000, 53 (3) : 293-298

# seasonality of births
seasonality <- 1 # set the relative strength of seasonality of births 
# this informs the daily birthrate: N_0 * alpha * (1 + (seasonality * (cos(2 * pi * tt / 360)))) in the model
# 1 = full strength reflecting that observed in KSA, 0 = not seasonal)

# the age dependent removal rate - chosen to balance birthrate for demographic equilibrium
mu_1st_yr <- 0.0011 # default 0.0011 death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # default 0.0011 death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # default 0.0003603 death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # default 0.0003603 death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # default 0.0003603 death rate in adulthood (>4 years)

mu <- vector(length = N_age)
mu[1:12] <- mu_1st_yr
mu[13:24] <- mu_2nd_yr
mu[25:36] <- mu_3rd_yr
mu[37:48] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

#############################
# transmission and recovery #
#############################

# the transmission rate
beta <- 7/14 # default 7/14 - corresponds to moderate transmission setting, R0 = 7, when gamma = 14

# the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 14
gamma <- 1/duration_infection

###########################
# relative infectiousness #
###########################

# default values are based on https://www.nature.com/articles/s41598-019-52730-4 
# and the assumption that viral load is proportional to infectiousness

# relative infectiousness of previously infected animals
reduced_shed <- 0.01 
# relative infectiousness of vaccinated animals
v_shed <- reduced_shed
# relative infectiousness of previously infected vaccinated animals
v_reduced_shed <- 0.0015

############
# immunity #
############

mAb_susc <- 0 # default = 0, proportion of susceptibility experienced by calves with mAbs
waning <- 1/30 # average rate at which complete immunity wanes following infection (/day)
Ab_susc <- 0.75 # relative susceptibility following complete immunity waning
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
mAb_waning <- par_esti["sigma_m", "mode_ms"] / 360 # read in mAbs waning from catalytic modelling


###############
# vaccination #
###############

# note if you want vaccination off, set coverage = 0
# otherwise vaccination starts immediately at the begining of your simulation period, after the 20 year set-up period

# vaccine impact
v_gamma <- gamma # default same as unvaccinated -  rate of recovery from infection in vaccinated animals
v_sigma <- waning # default same as unvaccinated - rate of waning of complete infection induced immunity in vaccinated animals
v_sigma_m <- mAb_waning # default same as unvaccinated rate of waning of mAbs in vaccinated animals
v_mAb_susc <- mAb_susc # default same as unvaccinated - relative susceptibility of vaccinated animals with mAbs
v_susc <- Ab_susc # default same as natural infection
v_Ab_susc <- 0.5 # default 0.5
rho <- c(1/(10*360)) # default 1/ 10 years - the rate at which vaccine induced affects wane

# vaccine implementation
age_targ_idx <- 6 # target the 6 month olds (NOTE - indexes 1:48 denote 1:48 months old respectively, and 49 denotes all animals >4 years)
coverage <- 0.8 # default 0.8 - vaccination coverage of the target age group

vaxp <- rep(0, 49)
vaxp[age_targ_idx] <- coverage # setting up your specified coverage and target age - do not edit


########################
# initialise the model #
########################

x <- sp_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = waning, sigma_m = mAb_waning, Ab_susc = Ab_susc, 
                       mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu = mu, N_0 = N_0,
                       importation_rate = importation_rate, imp_t = imp_t, delta = seasonality, ind1 = ind1, ind2 = ind2,
                       v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_susc  = v_susc, v_mAb_susc = v_mAb_susc, 
                       v_Ab_susc = v_Ab_susc, v_shed = v_shed, v_reduced_shed = v_reduced_shed, 
                       vaxp = vaxp, rho = rho,
                       foi_bg_usr = foi_bg_usr)


##############################################
# set the number of stochastic runs you want #
##############################################

n_runs <- 100 ## (e.g. this will run the model 100 times in parallel)

###############################
# select the outputs you want #
###############################

outputs <- 1:814 # default selects all model outputs stratified by age, otherwise specify the indexes of the outputs you want 
# e.g. if you want just the total number of infections across all ages, Itot, that is 606, if you want total incidence, that is 613.
# if you run once with all variables and then look at names(out) you will be able to use this to see when indexes you want.

#################
# run the model #
#################

# example takes 3.5 minutes on my local machine
out <- as.data.frame(replicate(n_runs, x$run(t)[,outputs]))

#################################################
# example plot - number of infections over time #
#################################################

# NOTE - discarding model outputs from 0-3600 days as this is the set up time with background FoI

outPlot <- out[3601:10801, c(grepl(names(out), pattern = "Itot"))]
outPlot$time <- as.numeric(rownames(outPlot))
outPlot <- outPlot%>%
  pivot_longer(names_to = "run", values_to = "infections", cols = 1:100)%>%
  mutate(run = str_extract(run, pattern = "[0-9]+"))


ggplot(outPlot)+
  geom_line(aes(x = (time - 3600)/360, y = infections, group = as.factor(run)), col = "forestgreen", alpha = 0.2)+
  xlab("time (years)")+
  geom_vline(xintercept = 10, col = "black", lwd = 1, lty = 2)+ #time vaccination introduced (at t = 7201)
  theme_bw()


