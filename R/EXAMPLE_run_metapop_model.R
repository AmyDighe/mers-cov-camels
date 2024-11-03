
# this script is an example script that can be adapted by the user to run their
# own simulations using the stochastic metapopulation model of MERS-CoV transmission 
# and vaccination in dromedary camels.

# In order to explore the effect population structure has on dynamics, we developed
# a rudimentary structured population model where sub-populations or patches are arranged 
# over a grid (Figure S5). Individuals are most likely to be in contact with other 
# individuals in the same patch, less likely to meet individuals in neighbouring patches,
# and do not meet individuals in distant patches. Patches experience different degrees
# of connection depending on whether they are in a corner, on edge or within the middle
# of the grid. In reality, herds will also have varying levels of connection to other
# herds or larger sub-populations. 

# This example depends on the model "
#     - dynamic-odin-models/metapop_model_vax.R"
#     - the estimated rate of waning of maternal antibodies from the catalytic modelling
#       "fits/processed_real/exp/sens_spec_1/fit4bb.rds" but you can replace this with any value
#       of your choice

# Note that this example runs the model for a single set of parameters, but odin.dust has fantastic
# ways of easily and efficiently updating parameters for comparing different parameter sets 
# (e.g. conducting sensitivity analyses, running under different transmission settings or comparing 
# vaccine coverage). An examples can be seen in 11_estimate_persistence_with_vax.R, and for more detailed
# instructions see the documentation here: https://mrc-ide.github.io/odin.dust/articles/sir_models.html

# This script generates a data frame "out" with 1 row per day and one column per model output per stochastic run
# so for the default example for 30 years with all 814 outputs and 100 stochastic model runs the dimensions of "out" are 10801x81400

#######################
## load dependencies ##
#######################

library(here) # this helps make file paths specific to your directory so the same code can be used by different people
source(here("dependencies.R"))
source(here("R/stochastic_init.R")) # function to initiate model at demographic equilibrium

###################
## compile model ##
###################
mp_model <- odin.dust::odin_dust("dynamic-odin-models/metapopvax.R", verbose = FALSE, skip_cache = TRUE)

#########################
## specify time period ##
#########################

# input the time period that you wish to simulate transmission in (in days)
# NOTE - I tend to use 360 as an approximation for a year as this makes it easier to align ageing of camels and monthly interventions
time_period <- 10*360 # e.g ten years
t <- seq(0:(time_period + 20*360)) # add this time_period to the set up period of 20 years (10 years with background FoI on, 10 years to reach equilibrium)

#######################
## specify structure ##
#######################

n_r <- 5 # number of rows in grid of sub-populations
n_c <- 5 # number of cols in grid of sub-populations
N_patch <- n_r*n_c # total number of sub-populations

# correction for balancing external outside of sub-population vs internal within
# sub-population force of infection
correction_ex <- matrix(c(2, rep(3, n_c - 2), 2,
                          rep(c(3, rep(4, n_c - 2), 3), n_r - 2),
                          2, rep(3, n_c - 2), 2), nrow = n_r, ncol = n_c)


########################
## specify parameters ##
########################

##############
# demography #
##############
N_0 <- 1000000 # set the total size of initial population
patch_pop <- round(N_0 / N_patch, 0)

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


######################
# initiate the model #
######################

# number of stochastic runs (run in parallel)
n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha = alpha, delta = seasonality, 
                      N_0 = patch_pop,
                      mu = mu, 
                      N_age = N_age, 
                      n_r = n_r, n_c = n_c)

storage.mode(S_ini_p) <- "double"

x <- mp_model$new(pars = list(N_age = N_age, nr = n_r, nc = n_c, 
                             alpha = alpha, 
                             beta = beta,
                             gamma = gamma, 
                             sigma = waning,
                             sigma_m = mAb_waning, 
                             Ab_susc = Ab_susc,
                             mAb_susc = mAb_susc, 
                             reduced_shed = reduced_shed, 
                             mu = mu, 
                             N_0 = patch_pop, 
                             importation_rate = importation_rate, 
                             imp_t = imp_t, 
                             delta =  seasonality, 
                             connectivity = 0.01,
                             correction_ex = correction_ex,
                             S_ini_p = S_ini_p,
                             foi_bg_usr = foi_bg_usr,
                             v_gamma = v_gamma, 
                             v_sigma = v_sigma, 
                             v_sigma_m = v_sigma_m, 
                             v_susc  = v_susc, 
                             v_mAb_susc = v_mAb_susc, 
                             v_Ab_susc = v_Ab_susc, 
                             v_shed = v_shed, 
                             v_reduced_shed = v_reduced_shed, 
                             vaxp = rep(0, 49), # no vaccination at initiation
                             rho = rho
), 
time = 1, 
n_particles = n_particles, 
n_threads = 2L, 
seed = 1L)

x$set_index(c(x$info()$index$Itot)) # outputting total number of infections
# add any other outputs of interest to you here. e.g. if you were also keen to
# see the total number of infections per sub-population, you could add 
# x$info()$index$Itot_patch into c()

# time steps to run the model (in days)
# I'm splitting this into two as I want to be able to start vaccination during steps_2
steps <- seq(1, 20*360, by = 10)# up to the point of vaccination
steps_2 <- seq(20*360, 30*360, by = 10) # after vaccination

# run the model for the first steps, with no vaccination
tictoc::tic()
out <- x$simulate(steps) # this takes ~6 minutes on my local machine
tictoc::toc()

saved_state <- x$state()

Itot_mat_initial <- (t(out[1,,]))

# update parameters and run the model forwards from the saved_state for steps 2, 
# with vaccination

x$update_state(pars = list(N_age = N_age, nr = n_r, nc = n_c, 
                           alpha = alpha, 
                           beta = beta,
                           gamma = gamma, 
                           sigma = waning,
                           sigma_m = mAb_waning, 
                           Ab_susc = Ab_susc,
                           mAb_susc = mAb_susc, 
                           reduced_shed = reduced_shed, 
                           mu = mu, 
                           N_0 = patch_pop, 
                           importation_rate = importation_rate, 
                           imp_t = imp_t, 
                           delta =  seasonality, 
                           connectivity = 0.01,
                           correction_ex = correction_ex,
                           S_ini_p = S_ini_p,
                           foi_bg_usr = foi_bg_usr,
                           v_gamma = v_gamma, 
                           v_sigma = v_sigma, 
                           v_sigma_m = v_sigma_m, 
                           v_susc  = v_susc, 
                           v_mAb_susc = v_mAb_susc, 
                           v_Ab_susc = v_Ab_susc, 
                           v_shed = v_shed, 
                           v_reduced_shed = v_reduced_shed, 
                           vaxp = vaxp, # vaccination of your target age group
                           rho = rho), 
               time = 20*360,
               state = saved_state,
               set_initial_state = FALSE)

state <- x$simulate(steps_2) # takes ~3 minutes on my local machine

Itot_mat_postVax <- (t(state[1,,]))

#################################################
# example plot - number of infections over time #
#################################################

# reshape data for plotting
ItotPlot <- as.data.frame(rbind(Itot_mat_initial, Itot_mat_postVax[2:361,])) # get rid of the starting state which is already in Itot_mat_initial
ItotPlot$time <- seq(1, 30*360, by = 10)

ItotPlot <- ItotPlot %>%
  pivot_longer(cols = 1:100, names_to = "run", values_to = "infections")%>%
  mutate(run = str_extract(run, pattern = "[0-9]+"))

ggplot(ItotPlot)+
  geom_line(aes(x = (time - 3600)/360, y = infections, group = as.factor(run)), col = "forestgreen", alpha = 0.2)+
  xlab("time (years)")+
  geom_vline(xintercept = 10, col = "black", lwd = 1, lty = 2)+ #time vaccination introduced (at t = 7201)
  theme_bw()
