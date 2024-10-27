# This script estimates the critical community size - the population size above 
# which the extinction of MERS-CoV transmission by chance becomes unlikely, which
# we define as the population size needed for transmission to persist in <50% of
# stochastic model runs. The drop off from persistence to extinction across all model runs
# happens very quickly over a relatively small change in population size.

# it depends on:
#  1. the single patch transmission model "dynamic-odin-models/single_patch_model.R"
#  2. the estimated rate of waning of mAbs from the best fitting catalytic model
#     "fits/processed_real/exp/sens_spec_1/fit4bb.rds

# it generates:
#  1. "generated_data/persistence.rds" - the proportion of stochastic model runs 
#      in which transmission persists used to estimate the CCS using interpolation
#      in `generate_Fig2.R`

source("dependencies.R")


# load single patch dynamic transmission model
sir_model <- odin::odin("dynamic-odin-models/single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input the time period that you wish to run the model for (in days)
time_period <- 35*360 
t <- seq(0:time_period)

# read in mAbs waning from catalytic modelling
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
omega <- par_esti["sigma_m", "mode_ms"] / 360

# set up parameter grid to run over mutliple scenarios and population sizes
shedding <- c("shedding_0.01" = 0.01,
              "shedding_0.25" = 0.25, 
              "shedding_0.50" = 0.5)
beta <- list("shedding_0.01" = c(3.5, 7, 14)/14, # R0 values selected based on Fig S4
             "shedding_0.25" = c(2.0, 3.5, 5.0)/14, # R0 values selected based on Fig S4
             "shedding_0.50" = c(1.75, 2.3, 3.0)/14) # R0 values selected based on Fig S4
waning <- c(1/30, 1) # two alternative rates of waning of antibodies following infection
susc <- c(0.75, 1) # two alternative rates of susceptibility of recovered animals
seasonality <- c(0.5, 1) # two alternative strengths of seasonality

# select populations sizes - depends on beta to give maximum resolution around the CCS
pop <- list("shedding_0.01" = c(c(5, 7.5) %o% 10^3, seq(1,7.5,0.5) %o% 10^4, 1 %o% 10^(5:6)),
             "shedding_0.25" = c(5 %o% 10^2, seq(1,9, 1) %o% 10^3, seq(1,5, 0.5) %o% 10^4, 1 %o% 10^5),
             "shedding_0.50" = c(5 %o% 10^2, seq(1,9, 1) %o% 10^3, c(seq(1,5, 0.5), 7.5) %o% 10^4, 1 %o% 10^5))

par_grid_core_0.01 <- expand.grid(pop = pop[["shedding_0.01"]],
                                  shedding = shedding["shedding_0.01"],
                                  beta = beta[["shedding_0.01"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_core_0.25 <- expand.grid(pop = pop[["shedding_0.25"]],
                                  shedding = shedding["shedding_0.25"],
                                  beta = beta[["shedding_0.25"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_core_0.50 <- expand.grid(pop = pop[["shedding_0.50"]],
                                  shedding = shedding["shedding_0.50"],
                                  beta = beta[["shedding_0.50"]],
                                  waning = waning[1], 
                                  susc = susc[1],
                                  seasonality = seasonality)
par_grid_core <- rbind(par_grid_core_0.01,
                       par_grid_core_0.25,
                       par_grid_core_0.50)


persist <- vector(length = (dim(par_grid_core)[1]))

# run model for each parameter combination
for(i in 1:(dim(par_grid_core)[1])){

x <- sir_model$new(alpha = alpha, beta = par_grid_core$beta[i], gamma = gamma, 
                   sigma = par_grid_core$waning[i], sigma_m = omega, 
                   Ab_susc = par_grid_core$susc[i], mAb_susc = mAb_susc, 
                   reduced_shed = par_grid_core$shedding[i], mu = mu,
                   N_0 = par_grid_core$pop[i], importation_rate = importation_rate, 
                   imp_t = 1, delta = par_grid_core$seasonality[i], 
                   ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)


out <- as.data.frame(replicate(100, x$run(t)[, 354])) # run 100x and isolate total number infected (col 354)

persist[i] <- sum(out[360*35, ] > 0) # 35yrs after model initiation (25yrs after background FOI switched off)

print(paste0(i, "/", (dim(par_grid_core)[1])))

}

par_grid_core$persist <- persist

# save output
saveRDS(file = "generated_data/persistence.rds", object = par_grid_core)
