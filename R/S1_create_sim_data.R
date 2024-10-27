
#####################
#####################
## SIMULATING DATA ##
#####################
#####################

# the number and size of data-sets and the age classes are chosen to
# be the same as the real data

###############
###############
## par grids ##
###############
###############

k <- c(0.01, 0.1, 1.0)
sigma_m <- c(1/0.25, 1/0.5, 1/0.75)
sigma_r <- c(1/5, 1/10, 1/15)

# pars_global <- data.frame(k = k, omega = sigma_m, sigma_r = sigma_r)
# pars_global_long <- pivot_longer(pars_global, 1:3, names_to = "par",
#                                  values_to = "value")
# saveRDS(file = "data/simulated/pars_global.rds", pars_global_long)
# 
 foi_idx <- c("foi1", "foi2", "foi3")
# foi <- vector(mode = "list", length = 3)
# for(i in 1:3){
# foi[[i]] <- sample(c(runif(n = 15, min = 0.5, max = 1.0),
#          runif(n = 8, min = 1.1, max = 5.0)))}
# names(foi) <- foi_idx
# saveRDS(file = "data/simulated/foi.rds", foi)

# number of camels - use 1000 and same number of empty age classes as the real data
# N_CAMELS <- readRDS("data/real/N_CAMELS.rds") # STUDY X AGE class matrix of N
# N_CAMELS_SIM <-apply(N_CAMELS,2,function(x) ifelse(x > 0, 1000, x))
# saveRDS(file = here("data/simulated/N_CAMELS.rds"), N_CAMELS_SIM)

## read in saved parameters
pars_global_long <- readRDS(file = "data/simulated/pars_global.rds")
foi <- readRDS("data/simulated/foi.rds")
N_CAMELS_SIM <- readRDS("data/simulated/N_CAMELS.rds")
AGE_L <- readRDS("data/real/AGE_L.rds") # STUDY X AGE class matrix of lower age bound
AGE_U <- readRDS("data/real/AGE_U.rds") # STUDY X AGE class matrix of upper age bound

########################################################################################
## beta-binomial EXPONENTIAL AGE DISTRIBUTION ##
########################################################################################

AGE_U[which(AGE_U == 2.0)] <- 1.99999 # important as mortality boundary set to exactly 2yrs

###########
# model 1 #
###########
pars_1bb_exp <- expand.grid(foi = foi_idx,
                            k = k)
pars_1bb_exp$id <- 1:(dim(pars_1bb_exp)[1])
pars_1bb_exp$foi <- as.character(pars_1bb_exp$foi)
pars_1bb_exp <- purrr::map_dfr(seq_len(3), ~pars_1bb_exp)


simdat_1bb_exp <- vector(mode = "list", length = dim(pars_1bb_exp)[1])

for(i in 1:(dim(pars_1bb_exp)[1])){
  
  simdat_1bb_exp[[i]] <- sim_data_betabinom_exp(n_datasets = dim(AGE_U)[1], n_ages = dim(AGE_U)[2], 
                                                gamma = foi[[pars_1bb_exp$foi[i]]], sigma = 0, 
                                                omega = 12, mabs = -1, #doesn't matter what omega is
                                                N_camels = N_CAMELS_SIM, 
                                                age_upper = AGE_U, age_lower = AGE_L, od = pars_1bb_exp$k[i],
                                                sens = 0.99, spec = 1, mu_0 = mu_0, mu = mu)$simulated
  
}

saveRDS(file = "data/simulated/1bb_exp.rds", simdat_1bb_exp)
saveRDS(file = "data/simulated/pars_1bb_exp.rds", pars_1bb_exp)

###########
# model 2 #
###########
pars_2bb_exp <- rbind(expand.grid(sigma_r = sigma_r,
                                  foi = foi_idx[2]),
                      expand.grid(sigma_r = sigma_r[2],
                                  foi = foi_idx[-2]))
pars_2bb_exp <- purrr::map_dfr(seq_len(3), ~pars_2bb_exp) %>%
  mutate(k= rep(k, each = 5))
pars_2bb_exp$id <- 1:(dim(pars_2bb_exp)[1])
pars_2bb_exp$foi <- as.character(pars_2bb_exp$foi)
pars_2bb_exp <- purrr::map_dfr(seq_len(3), ~pars_2bb_exp)

simdat_2bb_exp <- vector(mode = "list", length = dim(pars_2bb_exp)[1])

for(i in 1:(dim(pars_2bb_exp)[1])){
  
  simdat_2bb_exp[[i]] <- sim_data_betabinom_exp(n_datasets = dim(AGE_U)[1], n_ages = dim(AGE_U)[2], 
                                                gamma = foi[[pars_2bb_exp$foi[i]]], sigma = pars_2bb_exp$sigma_r[i], 
                                                omega = 12, mabs = -1, #doesn't matter what omega is
                                                N_camels = N_CAMELS_SIM, 
                                                age_upper = AGE_U, age_lower = AGE_L, od = pars_2bb_exp$k[i],
                                                sens = 0.99, spec = 1, mu_0 = mu_0, mu = mu)$simulated
  
}

saveRDS(file = "data/simulated/2bb_exp.rds", simdat_2bb_exp)
saveRDS(file = "data/simulated/pars_2bb_exp.rds", pars_2bb_exp)

###########
# model 3 #
###########
pars_3bb_exp <- rbind(expand.grid(sigma_m = sigma_m,
                                  foi = foi_idx[2]),
                      expand.grid(sigma_m = sigma_m[2],
                                  foi = foi_idx[-2]))
pars_3bb_exp <- purrr::map_dfr(seq_len(3), ~pars_3bb_exp) %>%
  mutate(k= rep(k, each = 5))
pars_3bb_exp$id <- 1:(dim(pars_3bb_exp)[1])
pars_3bb_exp$foi <- as.character(pars_3bb_exp$foi)
pars_3bb_exp <- purrr::map_dfr(seq_len(3), ~pars_3bb_exp)

simdat_3bb_exp <- vector(mode = "list", length = dim(pars_3bb_exp)[1])

for(i in 1:(dim(pars_3bb_exp)[1])){
  
  simdat_3bb_exp[[i]] <- sim_data_betabinom_exp(n_datasets = dim(AGE_U)[1], n_ages = dim(AGE_U)[2], 
                                                gamma = foi[[pars_3bb_exp$foi[i]]], sigma = 0, 
                                                omega = pars_3bb_exp$sigma_m[i], mabs = 1,
                                                N_camels = N_CAMELS_SIM, 
                                                age_upper = AGE_U, age_lower = AGE_L, od = pars_3bb_exp$k[i],
                                                sens = 0.99, spec = 1, mu_0 = mu_0, mu = mu)$simulated
  
}

saveRDS(file = "data/simulated/3bb_exp.rds", simdat_3bb_exp)
saveRDS(file = "data/simulated/pars_3bb_exp.rds", pars_3bb_exp)

###########
# model 4 #
###########
pars_4bb_exp <- rbind(expand.grid(sigma_m = sigma_m,
                                  sigma_r = sigma_r[2],
                                  foi = foi_idx[2]),
                      expand.grid(sigma_m = sigma_m[2],
                                  sigma_r = sigma_r[-2],
                                  foi = foi_idx[2]), 
                      expand.grid(sigma_m = sigma_m[2],
                                  sigma_r = sigma_r[2],
                                  foi = foi_idx[-2]))
pars_4bb_exp <- purrr::map_dfr(seq_len(3), ~pars_4bb_exp) %>%
  mutate(k= rep(k, each = 7))
pars_4bb_exp$id <- 1:(dim(pars_4bb_exp)[1])
pars_4bb_exp$foi <- as.character(pars_4bb_exp$foi)
pars_4bb_exp <- purrr::map_dfr(seq_len(3), ~pars_4bb_exp)

simdat_4bb_exp <- vector(mode = "list", length = dim(pars_4bb_exp)[1])

for(i in 1:(dim(pars_4bb_exp)[1])){
  
  simdat_4bb_exp[[i]] <- sim_data_betabinom_exp(n_datasets = dim(AGE_U)[1], n_ages = dim(AGE_U)[2], 
                                                gamma = foi[[pars_4bb_exp$foi[i]]], sigma = pars_4bb_exp$sigma_r[i], 
                                                omega = pars_4bb_exp$sigma_m[i], mabs = 1,
                                                N_camels = N_CAMELS_SIM, 
                                                age_upper = AGE_U, age_lower = AGE_L, od = pars_4bb_exp$k[i], 
                                                sens = 0.99, spec = 1, mu_0 = mu_0, mu = mu)$simulated
  
}

saveRDS(file = "data/simulated/pars_4bb_exp.rds", pars_4bb_exp)
saveRDS(file = "data/simulated/4bb_exp.rds", simdat_4bb_exp)

