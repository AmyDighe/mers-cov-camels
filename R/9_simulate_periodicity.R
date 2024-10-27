# This script estimates the periodicity of MERS-CoV infections in camels based
# on our assumption that births follow the same seasonality as observed in KSA
# for a range of R0 values reflecting our range of FoI estimates across the different
# camel populations studied. 

# The periodicity is determined by assessing which lag (>100 days) maximises 
# the autocorrelation coefficient for the simulated infection time series 
# and confirmed by looking at the time series by eye.

# this script depends on:
#  1. the dynamic transmission model "dynamic-odin-models/single_patch_model.R"
#  2. the estimated rate of waning of mAbs from the best fitting catalytic model
#     "fits/processed_real/exp/sens_spec_1/fit4bb.rds
#  3. The acf function in the R stats package

# it generates:
#  1. "generated_data/period.rds" used in `generate_Fig3.R`


source("dependencies.R")

## Extract periodicity of infections for persistent runs
sir_model <- odin::odin("dynamic-odin-models/single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

# input the time period that you wish to run the model for (in days)
time_period <- 35*360 
t <- seq(0:time_period)

# read in mAbs waning from catalytic modelling
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
omega <- par_esti["sigma_m", "mode_ms"] / 360

# read in parameter grid with persistence results
results_sp <- readRDS(here("generated_data", "persistence.rds"))
results_sp <- results_sp %>% filter(persist>50) # filter out those which only persisted <50% of the time

# include new combination assuming births are NOT seasonal for comparison
# only include versions that match with combinations that persist when seasonality = 1
par_grid_period <- rbind(results_sp %>% filter(seasonality == 1), 
                         results_sp %>% filter(seasonality == 1) %>% 
                           mutate(seasonality = 0, persist = NA), 
                         results_sp %>% filter(seasonality == 1) %>% 
                           mutate(seasonality = 0.5, persist = NA))

persist <- vector(length = (dim(par_grid_period)[1]))
period <- matrix(ncol = (dim(par_grid_period)[1]), nrow = 100)

# run the model for the differnt scenarios reflected in par_grid
for(i in 1:(dim(par_grid_period)[1])){
  
  x <- sir_model$new(alpha = alpha, beta = par_grid_period$beta[i], gamma = gamma, 
                     sigma = par_grid_period$waning[i], sigma_m = omega, 
                     Ab_susc = par_grid_period$susc[i], mAb_susc = mAb_susc, 
                     reduced_shed = par_grid_period$shedding[i], mu = mu,
                     N_0 = par_grid_period$pop[i], importation_rate = importation_rate, 
                     imp_t = 1, delta = par_grid_period$seasonality[i], 
                     ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)
  
  
  out <- as.data.frame(replicate(100, x$run(t)[, 354])) # run 100x and isolate total number infected (col 354)
  
  persist[i] <- sum(out[360*35, ] > 0) # 35yrs after model initiation (25yrs after background FOI switched off)
  
  if(persist[i] > 50){
    
    idx_persist <- which(out[360*35, ] > 0)
    out_persist <- out[,idx_persist]
    
    for(j in 1:(persist[i])){
      ACF <- acf(out_persist[(360*25):(360*35), j], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[180:(360*5)], lag = ACF$lag[180:(360*5)])
      # write something in here about making sure acf is above the significance level
      period[j,i] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
  } else {
    period[,i] <- NA
  }
  
  print(paste0(i, "/", (dim(par_grid_period)[1])))
  
}

par_grid_period$persist2 <- persist

# save model output
saveRDS(file = "generated_data/period.rds", object = list(par_grid = par_grid_period,
                                                                           period = period))
