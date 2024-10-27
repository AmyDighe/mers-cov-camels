
################################################
# fit to simulated data to assess par recovery #
################################################
simdat_1bb <- readRDS(file = "data/simulated/1bb_exp.rds")
simdat_2bb <- readRDS(file = "data/simulated/2bb_exp.rds")
simdat_3bb <- readRDS(file = "data/simulated/3bb_exp.rds")
simdat_4bb <- readRDS(file = "data/simulated/4bb_exp.rds")

AGE_L <- readRDS("data/real/AGE_L.rds") # STUDY X AGE class matrix of lower age bound
AGE_U <- readRDS("data/real/AGE_U.rds") # STUDY X AGE class matrix of upper age bound
N_CAMELS <- readRDS("data/simulated/N_CAMELS.rds")

## MODEL 1 beta-binomial 
fit1bb_smry <- vector(mode = "list", length = length(simdat_1bb))
fit1bb_sim <- vector(mode = "list", length = length(simdat_1bb))

for(i in 1:length(simdat_1bb)){
  fit1bb_sim[[i]] <- stan(
    file = here::here("catalytic-stan-models/model1bb_exp.stan"),
    data = list(
      S = nrow(simdat_1bb[[i]]),
      A =  ncol(simdat_1bb[[i]]),
      N = N_CAMELS,
      pos = simdat_1bb[[i]],
      age1 = AGE_L,
      age2 = AGE_U,
      sigma_m = 12,
      sigma_r = 0,
      sens = rep(0.99, 23),
      spec = rep(1, 23),
      mabs = -1,
      mu_0 = mu_0,
      mu = mu
    ),
    chains = 4,
    iter = 10000,
    warmup = 2000,
    verbose = TRUE#,
    ##init = init_fun
    ##control = list(adapt_delta = 0.99)
  )
  fit1bb_smry[[i]] <- as.data.frame(rstan::summary(fit1bb_sim[[i]])$summary)
}
saveRDS(fit1bb_smry, file = "fits/sim_data/1bbsim.rds")
saveRDS(fit1bb_sim, file = "fits/sim_data/1bbsim_fit.rds")

## MODEL 2 beta-binomial
fit2bb_smry <- vector(mode = "list", length = length(simdat_2bb))
fit2bb_sim <- vector(mode = "list", length = length(simdat_2bb))

for(i in 1:length(simdat_2bb)){
  fit2bb_sim[[i]] <- stan(
    file = here::here("catalytic-stan-models/model2bb_exp.stan"),
    data = list(
      S = nrow(simdat_2bb[[i]]),
      A =  ncol(simdat_2bb[[i]]),
      N = N_CAMELS,
      pos = simdat_2bb[[i]],
      age1 = AGE_L,
      age2 = AGE_U,
      sigma_m = 12,
      sens = rep(0.99, 23),
      spec = rep(1, 23),
      mabs = -1,
      mu_0 = mu_0,
      mu = mu
    ),
    chains = 4,
    iter = 10000,
    warmup = 2000,
    verbose = TRUE
    ##control = list(adapt_delta = 0.99) 
  )
  fit2bb_smry[[i]] <- as.data.frame(rstan::summary(fit2bb_sim[[i]])$summary)
}
saveRDS(fit2bb_smry, file = "fits/sim_data/2bbsim.rds")
saveRDS(fit2bb_sim, file = "fits/sim_data/2bbsim_fit.rds")

## MODEL 3 beta-binomial
fit3bb_smry <- vector(mode = "list", length = length(simdat_3bb))
fit3bb_sim <- vector(mode = "list", length = length(simdat_3bb))

for(i in 1:length(simdat_3bb)){
  fit3bb_sim[[i]] <- stan(
    file = here::here("catalytic-stan-models/model3bb_exp.stan"),
    data = list(
      S = nrow(simdat_3bb[[i]]),
      A =  ncol(simdat_3bb[[i]]),
      N = N_CAMELS,
      pos = simdat_3bb[[i]],
      age1 = AGE_L,
      age2 = AGE_U,
      sigma_r = 0,
      sens = rep(0.99, 23),
      spec = rep(1, 23),
      mabs = 1,
      mu_0 = mu_0,
      mu = mu
    ),
    chains = 4,
    iter = 10000,
    warmup = 2000,
    verbose = TRUE,
    ##control = list(max_treedepth = 15)
    ##control = list(adapt_delta = 0.99)
  ) 
  fit3bb_smry[[i]] <- as.data.frame(rstan::summary(fit3bb_sim[[i]])$summary)
}
saveRDS(fit3bb_smry, file = "fits/sim_data/3bbsim.rds")
saveRDS(fit3bb_sim, file = "fits/sim_data/3bbsim_fit.rds")

## MODEL 4 beta-binomial

fit4bb_smry <- vector(mode = "list", length = length(simdat_4bb))
fit4bb_sim <- vector(mode = "list", length = length(simdat_4bb))

for(i in 1:length(simdat_4bb)){
  fit4bb_sim[[i]] <- stan(
    file = here::here("catalytic-stan-models/model4bb_exp.stan"),
    data = list(
      S = nrow(simdat_4bb[[i]]),
      A =  ncol(simdat_4bb[[i]]),
      N = N_CAMELS,
      pos = simdat_4bb[[i]],
      age1 = AGE_L,
      age2 = AGE_U,
      sens = rep(0.99, 23),
      spec = rep(1, 23),
      mabs = 1,
      mu_0 = mu_0,
      mu = mu
      # B = length(age1_fine),
      # age1_fine = age1_fine,
      # age2_fine = age2_fine
    ),
    chains = 4,
    iter = 10000,
    warmup = 2000,
    verbose = TRUE
  )
  
  fit4bb_smry[[i]] <- as.data.frame(rstan::summary(fit4bb_sim[[i]])$summary)
}

saveRDS(fit4bb_smry, file = "fits/sim_data/4bbsim.rds")
saveRDS(fit4bb_sim, file = "fits/sim_data/4bbsim_fit.rds")
