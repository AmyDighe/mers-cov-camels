######################################################################
# leave each study out in turn to assess impact on overall estimates #
######################################################################

# this script looks at the impact each seroprevalence dataset has on the overall 
# estimates of the force of infection and antibody waning, by leaving each out in turn. 
# The output is used to produce FigS2 and S3 in `generate_figS2S3.R`.

# read in reshaped data (output of 1_process_raw_data.R)
data_sero <- readRDS("data/real/data_sero.rds") 
n_datasets <- length(unique(data_sero$STUDY_COUNTRY))
data <- list()
data_reshaped <- list()
STUDY_TEST_TYPE <- list()
STUDY_TEST_TYPE_1 <- list()

for(i in 1:n_datasets){
  
  data[[i]] <- data_sero %>% filter(STUDY_COUNTRY!= unique(STUDY_COUNTRY)[i])
  
  data_reshaped[[i]] <- reshape_data(data[[i]])
  
  STUDY_TEST_TYPE[[i]] <- unique(data_reshaped[[i]]$data_sero %>% dplyr::select(STUDY_COUNTRY, TEST_TYPE))
  STUDY_TEST_TYPE[[i]]$id <- 1:nrow(STUDY_TEST_TYPE[[i]])
  
  STUDY_TEST_TYPE_1[[i]] <- merge(STUDY_TEST_TYPE[[i]], TEST_SPEC_SENS_1)
  STUDY_TEST_TYPE_1[[i]] <- STUDY_TEST_TYPE_1[[i]][order(STUDY_TEST_TYPE_1[[i]]$id),]
}

## MODEL 1
fit1bb_smry <- vector(mode = "list", length = 23)
fit1bb_exp <- vector(mode = "list", length = 23)

for(i in 1:n_datasets){
  fit1bb_exp[[i]] <- stan(
    file = here::here("stan-models/model1bb_exp.stan"),
    data = list(
      S = nrow(data_reshaped[[i]]$SEROPOS),
      A =  ncol(data_reshaped[[i]]$SEROPOS),
      N = data_reshaped[[i]]$N_CAMELS,
      pos = data_reshaped[[i]]$SEROPOS,
      age1 = data_reshaped[[i]]$AGE_L,
      age2 = data_reshaped[[i]]$AGE_U,
      sigma_m = 12,
      sigma_r = 0,
      sens = STUDY_TEST_TYPE_1[[i]]$SENS,
      spec = STUDY_TEST_TYPE_1[[i]]$SPEC,
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
  fit1bb_smry[[i]] <- as.data.frame(rstan::summary(fit1bb_exp[[i]])$summary)
}
saveRDS(fit1bb_smry, file = "fits/real_data/sens_spec_1/exp/LOO/1bb.rds")
saveRDS(fit1bb_exp, file = "fits/real_data/sens_spec_1/exp/LOO/1bb_fit.rds")


## MODEL 2
fit2bb_smry <- vector(mode = "list", length = 23)
fit2bb_exp <- vector(mode = "list", length = 23)

for(i in 1:n_datasets){
  fit2bb_exp[[i]] <- stan(
    file = here::here("stan-models/model2bb_exp.stan"),
    data = list(
      S = nrow(data_reshaped[[i]]$SEROPOS),
      A =  ncol(data_reshaped[[i]]$SEROPOS),
      N = data_reshaped[[i]]$N_CAMELS,
      pos = data_reshaped[[i]]$SEROPOS,
      age1 = data_reshaped[[i]]$AGE_L,
      age2 = data_reshaped[[i]]$AGE_U,
      sigma_m = 12,
      sens = STUDY_TEST_TYPE_1[[i]]$SENS,
      spec = STUDY_TEST_TYPE_1[[i]]$SPEC,
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
  fit2bb_smry[[i]] <- as.data.frame(rstan::summary(fit2bb_exp[[i]])$summary)
}
saveRDS(fit2bb_smry, file = "fits/real_data/sens_spec_1/exp/LOO/2bb.rds")
saveRDS(fit2bb_exp, file = "fits/real_data/sens_spec_1/exp/LOO/2bb_fit.rds")


## MODEL 3 beta-binomial
fit3bb_smry <- vector(mode = "list", length = 23)
fit3bb_exp <- vector(mode = "list", length = 23)

for(i in 1:n_datasets){
  fit3bb_exp[[i]] <- stan(
    file = here::here("stan-models/model3bb_exp.stan"),
    data = list(
      S = nrow(data_reshaped[[i]]$SEROPOS),
      A =  ncol(data_reshaped[[i]]$SEROPOS),
      N = data_reshaped[[i]]$N_CAMELS,
      pos = data_reshaped[[i]]$SEROPOS,
      age1 = data_reshaped[[i]]$AGE_L,
      age2 = data_reshaped[[i]]$AGE_U,
      sigma_r = 0,
      sens = STUDY_TEST_TYPE_1[[i]]$SENS,
      spec = STUDY_TEST_TYPE_1[[i]]$SPEC,
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
  fit3bb_smry[[i]] <- as.data.frame(rstan::summary(fit3bb_exp[[i]])$summary)
}
saveRDS(fit3bb_smry, file = "fits/real_data/sens_spec_1/exp/LOO/3bb.rds")
saveRDS(fit3bb_exp, file = "fits/real_data/sens_spec_1/exp/LOO/3bb_fit.rds")


## MODEL 4 beta-binomial
fit4bb_smry <- list()
fit4bb_exp <- vector(mode = "list", length = 23)

for(i in 1:n_datasets){
  fit4bb_exp[[i]] <- stan(
    file = here::here("stan-models/model4bb_exp.stan"),
    data = list(
      S = nrow(data_reshaped[[i]]$SEROPOS),
      A =  ncol(data_reshaped[[i]]$SEROPOS),
      N = data_reshaped[[i]]$N_CAMELS,
      pos = data_reshaped[[i]]$SEROPOS,
      age1 = data_reshaped[[i]]$AGE_L,
      age2 = data_reshaped[[i]]$AGE_U,
      sens = STUDY_TEST_TYPE_1[[i]]$SENS,
      spec = STUDY_TEST_TYPE_1[[i]]$SPEC,
      mabs = 1,
      mu_0 = mu_0,
      mu = mu
      # B = length(age1_fine),
      # age1_fine = age1_fine,
      # age2_fine = age2_fine
    ),
    chains = 4,
    iter = 5000,
    warmup = 2000,
    verbose = TRUE
  )
  
  fit4bb_smry[[i]] <- as.data.frame(rstan::summary(fit4bb_exp[[i]])$summary)
}

saveRDS(fit4bb_smry, file = "fits/real_data/sens_spec_1/exp/LOO/4bb.rds")
saveRDS(fit4bb_exp, file = "fits/real_data/sens_spec_1/exp/LOO/4bb_fit.rds")
