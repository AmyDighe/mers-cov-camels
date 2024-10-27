# This script estimates the critical community size - the population size above 
# which the extinction of MERS-CoV transmission by chance becomes unlikely, which
# we define as the population size needed for transmission to persist in <50% of
# stochastic model runs. The drop off from persistence to extinction across all model runs
# happens very quickly over a relatively small change in population size.

# it depends on:
#  1. the meta population model "dynamic-odin-models/metapop.R"
#  2. the estimated rate of waning of mAbs from the best fitting catalytic model
#     "fits/processed_real/exp/sens_spec_1/fit4bb.rds

# it generates:
#  1. "generated_data/persistence.rds" - the proportion of stochastic model runs 
#      in which transmission persists used to estimate the CCS using interpolation
#      in `generate_Fig2_Table2.R`

source("dependencies.R")

# compile metapopulation model using dust
metadust <- odin.dust::odin_dust("dynamic-odin-models/metapop.R")

# stochastic initialisation near to demographic equilibrium
source("R/stochastic_init.R") # function for starting model at zoographic equilibrium

#############################
## user defined variables: ##
#############################

n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of subpops

# correction for balancing external vs internal foi
correction_ex <- matrix(c(2, 3, 3, 3, 2,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          3, 4, 4, 4, 3,
                          2, 3, 3, 3, 2), nrow = 5, ncol = 5)

# connectivity between patches
connectivity <- 0.01

# input the time period that you wish to run the model for (in days)
time_period <- 35*360 
t <- seq(0:time_period)

# read in mAbs waning from catalytic modelling
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
omega <- par_esti["sigma_m", "mode_ms"] / 360

# set up parameter grid
shedding <- c("shedding_0.01" = 0.01,
              "shedding_0.25" = 0.25, 
              "shedding_0.50" = 0.5)
beta <- list("shedding_0.01" = c(3.5, 7, 14)/14,
             "shedding_0.25" = c(2.0, 3.5, 5.0)/14,
             "shedding_0.50" = c(1.75, 2.3, 3.0)/14)
waning <- c(1/30, 1)
susc <- c(0.75, 1)
seasonality <- c(0.5, 1)

# pop same as single patch CCS script divided by 25 (it is fed into the model as the per-patch population)
pop <- list("shedding_0.01" = c(c(5, 7.5) %o% 10^3, seq(1,7.5,0.5) %o% 10^4, 1 %o% 10^(5:6))/25,
            "shedding_0.25" = c(5 %o% 10^2, seq(1,9, 1) %o% 10^3, seq(1,5, 0.5) %o% 10^4, 1 %o% 10^5)/25,
            "shedding_0.50" = c(5 %o% 10^2, seq(1,9, 1) %o% 10^3, c(seq(1,5, 0.5), 7.5) %o% 10^4, 1 %o% 10^5)/25)

# pars for each shedding assumption
par_grid_core_0.01mp <- expand.grid(pop = pop[["shedding_0.01"]],
                                  shedding = shedding["shedding_0.01"],
                                  beta = beta[["shedding_0.01"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_core_0.25mp <- expand.grid(pop = pop[["shedding_0.25"]],
                                  shedding = shedding["shedding_0.25"],
                                  beta = beta[["shedding_0.25"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_core_0.50mp <- expand.grid(pop = pop[["shedding_0.50"]],
                                  shedding = shedding["shedding_0.50"],
                                  beta = beta[["shedding_0.50"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_coremp <- rbind(par_grid_core_0.01mp,
                       par_grid_core_0.25mp,
                       par_grid_core_0.50mp)

###################
## run the model ##
###################

# benchmarked at about 10 minutes per run --> long run time overall ~24 hours per shedding
# could maybe just sample every 100 steps for persistence

# initiate model
n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha, delta = par_grid_core_mp$seasonality[1], 
                      N_0 = par_grid_core_mp$pop[1],
                      mu, N_age, n_r = n_r, n_c = n_c)

storage.mode(S_ini_p) <- "double"


msirs_model <- metadust$new(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                               alpha = alpha, 
                                               beta = par_grid_core_mp$beta[1],
                                               gamma = gamma, 
                                               sigma = par_grid_core_mp$waning[1],
                                               sigma_m = omega, 
                                               Ab_susc = par_grid_core_mp$susc[1],
                                               mAb_susc = mAb_susc, 
                                               reduced_shed = par_grid_core_mp$shedding[1], mu = mu, 
                                               N_0 = par_grid_core_mp$pop[1], 
                                               importation_rate = importation_rate, 
                                               imp_t = imp_t, 
                                               delta =  par_grid_core_mp$seasonality[1], 
                                               connectivity = connectivity,
                                               correction_ex = correction_ex,
                                               S_ini_p = S_ini_p,
                                               foi_bg_usr = foi_bg_usr), 
                                   step = 0, 
                                   n_particles = n_particles, 
                                   n_threads = 2L, 
                                   seed = 1L)

msirs_model$set_index(msirs_model$info()$index$Itot) # just extract Itot


# PERSISTENCE ANALYSIS

## number of model runs
n_particles <- 100L

## stochastic initialisation
S_ini_list <- list(length = dim(par_grid_coremp)[1]) ## stochastic initialization
for(i in 1:(dim(par_grid_coremp)[1])){
  S_ini_p <- stoch_init(alpha, delta = par_grid_coremp$seasonality[i],
                        N_0 = par_grid_coremp$pop[i],
                        mu, N_age, n_r = n_r, n_c = n_c)
  storage.mode(S_ini_p) <- "double"
  S_ini_list[[i]] <- S_ini_p
}

  persist <- vector(length = (dim(par_grid_coremp)[1]))
  steps <- seq(0, 12700, by = 100)
  
## set up the model
  for(i in 1:(dim(par_grid_coremp)[1])){
    msirs_model$update_state(pars = list(N_age = N_age, n_r = n_r, n_c = n_c, 
                                         alpha = alpha, 
                                         beta = par_grid_coremp$beta[i],
                                         gamma = gamma, 
                                         sigma = par_grid_coremp$waning[i],
                                         sigma_m = omega, 
                                         Ab_susc = par_grid_coremp$susc[i],
                                         mAb_susc = mAb_susc, 
                                         reduced_shed = par_grid_coremp$shedding[i], mu = mu, 
                                         N_0 = par_grid_coremp$pop[i], 
                                         importation_rate = importation_rate, 
                                         imp_t = imp_t, 
                                         delta =  par_grid_coremp$seasonality[i], 
                                         connectivity = connectivity,
                                         correction_ex = correction_ex,
                                         S_ini_p = S_ini_list[[i]],
                                         foi_bg_usr = foi_bg_usr), step = 0)
    
    msirs_model$set_index(msirs_model$info()$index$Itot) # just extract Itot
    state <- msirs_model$simulate(steps)
    out <- as.data.frame(t(drop(state[1,,])))
    
    persist[i] <- sum(out[(360*35)/100, ] > 0)
    print(paste("i = ", i, "persist = ",persist[i], sep = " "))
    
  }
  
  saveRDS(file = "generated_data/persistence.rds", 
          object = persist)
  
# save outputs
par_grid_core_mp$persist <- persist
saveRDS(par_grid_core_mp, file = here("generated_data", "persistence_mp.rds"))
