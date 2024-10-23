##############################
# 1. fixed global parameters #
##############################

######################################
## test sensitivity and specificity ##
######################################
data_sero <- readRDS("data/real/data_sero.rds") # dataframe STUDY * VARIABLES
TEST_TYPES <- unique(data_sero %>% dplyr::select(TEST_TYPE))

# assuming high sensitivity and specificity
TEST_SPEC_SENS_1 <- data.frame(TEST_TYPES, 
                               SENS = c(rep(0.980, 4), 0.980*0.980),
                               SPEC = c(rep(0.995, 2), rep(0.985, 2), 
                                        0.995 + 0.985 - (0.995 * 0.985)))
# assuming lower sensitivity of neutralisation based tests
TEST_SPEC_SENS_2 <- data.frame(TEST_TYPES, 
                               SENS = c(rep(0.850, 2), rep(0.980, 2),
                                        0.850*0.980),
                               SPEC = c(rep(0.995, 2), rep(0.985,2), 
                                        0.995 + 0.985 - (0.995 * 0.985)))
#############################################
## mortality rates for catalytic modelling ##
#############################################
mu_0 <- 0.0011 * 360 # daily mortality rate for individuals <2 yrs
mu <- 0.0003603 * 360 # daily mortality rate for individuals >2 yrs

####################
## age resolution ##
####################
# used in generated quantities - to produce predicted seroprev. for fine ages
# for analysing visual cataytic model fits (Figure X)
age1_fine <- seq(0, 20, by = 0.1)
age2_fine <- age1_fine + 0.099 

#####################################
## dynamic model  fixed parameters ##
#####################################

# number of age classes
N_age <- 49

# mean birth rate (per camel per day) 
alpha <- 0.000565 # 90% female * 50% reproductive age * 45.2% fecundity

# the average duration of the infectious period (in days) 
duration_infection <- 14 # default = 14
gamma <- 1/duration_infection

# proportion of susceptibility experienced by calves with mAbs
mAb_susc <- 0 # default = 0

# the age dependent removal rate - balance birthrate
mu_1st_yr <- 0.0011 # death rate for 1st year of life = 40% removal
mu_2nd_yr <- 0.0011 # death rate for 2nd year of life
mu_3rd_yr <- 0.0003603 # death rate for 3rd year of life = 14% removal
mu_4th_yr <- 0.0003603 # death rate for 4th year of life
mu_adult_over_4 <- 0.0003603 # death rate in adulthood (>4 years)

mu <- vector(length = N_age)
mu[1:12] <- mu_1st_yr
mu[13:24] <- mu_2nd_yr
mu[25:36] <- mu_3rd_yr
mu[37:48] <- mu_4th_yr
mu[N_age] <- mu_adult_over_4

# importation rate for introducing infectious individuals
importation_rate <- 0

# background force of infection for first ten years
foi_bg_usr <- 0.000015
imp_t <- 1

##############################################
# stochastic initialisation of dynamic model #
##############################################

# index for summing births for each age class
ind1 <- rep(0,12)
ind2 <- rep(0,12)

for(y in 2:13){
  ind1[y-1] <- 360 - ((y - 1) * 30) + 1 
  ind2[y-1] <- 360 - ((y - 2) * 30)
}

# repeating for 4 years to cover all age classes
ind1 <- rep(ind1, 4)
ind2 <- rep(ind2, 4)


#############################
## study names for figures ##
#############################

pretty_names <- c(  "Bangladesh (Islam et al 2015)",
                    "Egypt (Ali et al 2017b)", 
                    "Egypt (Kandeil et al 2019)", 
                    "Ethiopia (Reusken et al 2014)",
                    "Iraq (Kandeil et al 2019)",
                    "Iraq (Thwiny et al 2018)",
                    "Jordan (Kandeil et al 2019)",
                    "Jordan (van Doremalen et al 2017)",
                    "Kenya (Deem et al 2015)",
                    "Kenya (Munyua et al 2017)",
                    "Kenya (Ommeh et al 2018)",
                    "KSA (Alagaili et al 2014 - 1992)",
                    "KSA (Alagaili et al 2014 - 2013)",
                    "KSA (Harrath et al 2018)",
                    "KSA (Hemida et al 2013)",
                    "KSA (Kandeil et al 2019)",
                    "Pakistan (Saqib et al 2016)",
                    "Pakistan (Zohaib et al 2018)",
                    "Senegal (Kandeil et al 2019)",
                    "Tunisia (Kandeil et al 2019)",
                    "Tunisia (Reusken et al 2014)",
                    "UAE (Wernery et al 2015)",
                    "Uganda (Kandeil et al 2019)") 


################ 
# 2. functions #
################

# get CIs of data points

ci_lower <- function(x, n){
  
  lower <- binom::binom.confint(x = x, n = n, method = "exact")$lower
  return(lower)
}

ci_upper <- function(x, n){
  
  upper <- binom::binom.confint(x = x, n = n, method = "exact")$upper
  return(upper)
}

# reshape data for stan

reshape_data <- function(data_sero){
  
  # restructure the data into matrices for stan
  
  n_datasets <- length(unique(data_sero$STUDY_COUNTRY))
  n_ages <-  max(dplyr::count(data_sero, STUDY, COUNTRY)$n)
  
  AGE_L <- matrix(data = 21, nrow = n_datasets, ncol = n_ages)
  AGE_U <- matrix(data = 22, nrow = n_datasets, ncol = n_ages)
  N_CAMELS <- matrix(data = 0, nrow = n_datasets, ncol = n_ages)
  SEROPOS <- matrix(data = 0, nrow = n_datasets, ncol = n_ages)
  row_names <- unique(data_sero$STUDY_COUNTRY)
  rownames(AGE_L) <- row_names
  rownames(AGE_U) <- row_names
  rownames(N_CAMELS) <- row_names
  rownames(SEROPOS) <- row_names
  col_names <- 1:n_ages
  colnames(AGE_L) <- col_names
  colnames(AGE_U) <- col_names
  colnames(N_CAMELS) <- col_names
  colnames(SEROPOS) <- col_names
  
  # fill out these matrices with the ragged arrays 
  study_country <- unique(data_sero$STUDY_COUNTRY)
  
  for(s in 1:n_datasets){
    # lower age
    c_s_agel <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$LOW_AGE
    names(c_s_agel) <- seq(1: length(c_s_agel))
    
    AGE_L[study_country[s], names(c_s_agel)] <- c_s_agel
    
    # upper age
    c_s_ageu <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$UPP_AGE
    names(c_s_ageu) <- seq(1: length(c_s_ageu))
    
    AGE_U[study_country[s], names(c_s_ageu)] <- c_s_ageu
    
    # N_camels
    c_s_ncam <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$SERO_N
    names(c_s_ncam) <- seq(1: length(c_s_ncam))
    
    N_CAMELS[study_country[s], names(c_s_ncam)] <- c_s_ncam
    
    # sero_pos
    c_s_spos <- (data_sero %>% filter(STUDY_COUNTRY == study_country[s]))$SERO_POS
    names(c_s_spos) <- seq(1: length(c_s_spos))
    
    SEROPOS[study_country[s], names(c_s_spos)] <- c_s_spos
    
  }
  
  # add seroprevalnce and ci
  
  data_sero$AGE_MID <- data_sero$LOW_AGE + (data_sero$UPP_AGE - data_sero$LOW_AGE)/2
  
  data_sero$seroprevalence <- data_sero$SERO_POS / data_sero$SERO_N
  
  data_sero$ci_low <- apply(data_sero[,c("SERO_POS", "SERO_N")], 1, 
                            function(x) ci_lower(x[1], x[2]))
  
  data_sero$ci_upp <- apply(data_sero[,c("SERO_POS", "SERO_N")], 1, 
                            function(x) ci_upper(x[1], x[2]))
  
  return(list(data_sero = data_sero,
              SEROPOS = SEROPOS,
              AGE_L = AGE_L,
              AGE_U = AGE_U,
              N_CAMELS = N_CAMELS))
  
}

# model comparison using DIC
DIC <- function(fit){ # based on Spiegelhalter 2002
  
  fits <- rstan::extract(fit)
  ll <- fits$log_lik # get array of ll per iteration per study per age class (total iterations include all chains)
  ll <- rowSums(ll) #sum across studies & ages (zeros from ragged arrays shouldn't matter)
  
  D_bar <- -2*mean(ll) # mean deviance - by this point there are no zeros etc. as they are absorbed
  idx_mode <- which.max(fits$lp__) # index at which log posterior is maximal
  pD_mode <- D_bar - (-2*ll[idx_mode])
  
  DIC_mode <- D_bar + pD_mode
  
  return(list(D_bar = D_bar,
              pD_mode = pD_mode,
              DIC_mode = DIC_mode))
  
}

DIC_across_chains <- function(fit, no_chains, total_length, burnin){
  # extract fit and isolate log likelihood
  fit <- rstan::extract(fit)
  ll <- fit$log_lik
  ll_per_iteration <- rowSums(ll) # sum across all data sets
  
  chain_length <-  total_length - burnin
  
  ll_mat <- matrix(ll_per_iteration, nrow = chain_length, ncol = no_chains, byrow = F)
  lp_mat <- matrix(fit$lp__, nrow = chain_length, ncol = no_chains, byrow = F)
  
  D_bar <- vector(length = no_chains)
  pD <- vector(length = no_chains)
  DIC <- vector(length = no_chains)
  idx_mode <- vector(length = no_chains)
  
  for(i in 1:no_chains){
    D_bar[i] <- -2*mean(ll_mat[, i])
    idx_mode[i] <- which.max(lp_mat[, i])
    pD[i] <- D_bar[i] - (-2*ll_mat[idx_mode[i], i])
    DIC[i] <- D_bar[i] + pD[i]
  }
  
  return(list(chain_DIC = tibble(DIC, D_bar, pD),
              av_DIC = tibble(mean(DIC), mean(D_bar), mean(pD))))
  
}

# generate summary of parameter estimates including mode

gen_summary <- function(fit){
  
  fit_sum <- as.data.frame((rstan::summary(fit))$summary)
  fit <- as.data.frame(rstan::extract(fit)) %>%
    rename(od = k)
  
  mode_est_default <- vector(length = dim(fit_sum)[1])
  mode_est_ms <- vector(length = dim(fit_sum)[1])
  
  name_select <- colnames(fit)[grep(colnames(fit), pattern = "foi|sigma|od")]
  
  mode_est_default[1:(length(name_select))] <- apply(X = as.matrix(fit[,name_select]),
                                                     MARGIN = 2, FUN = mlv, simplify = TRUE)
  
  mode_est_ms[1:(length(name_select))] <- apply(X = as.matrix(fit[, name_select]),
                                                MARGIN = 2, FUN = mlv, simplify = TRUE, method = "meanshift")
  
  fit_sum$mode <- mode_est_default
  fit_sum$mode_ms <- mode_est_ms
  
  return(fit_sum)
}

extract_foi <- function(fit){
  foifit <- fit[grep("foi", rownames(fit)), c("mode_ms", "2.5%", "97.5%")]%>%
    rename(low = '2.5%',
           high = '97.5%')
  foifit <- foifit %>%
    mutate(foi = paste0(round(mode_ms, 2), " (", round(low, 2), 
                        ", ", round(high, 2), ")"))
  
  if("sigma_r" %in% rownames(fit)){
    sigma_r <- paste0(round(fit["sigma_r", "mode_ms"], 2), " (", round(fit["sigma_r", "2.5%"], 2), 
                ", ", round(fit["sigma_r", "97.5%"], 2), ")")
  } else {
    sigma_r <- NA
  }
  
  if("sigma_m" %in% rownames(fit)){
    sigma_m <- paste0(round(fit["sigma_m", "mode_ms"], 2), " (", round(fit["sigma_m", "2.5%"], 2), 
                ", ", round(fit["sigma_m", "97.5%"], 2), ")")
  } else {
    sigma_m <- NA
  }
  
  k <- paste0(round(fit["k", "mode_ms"], 2), " (", round(fit["k", "2.5%"], 2), 
              ", ", round(fit["k", "97.5%"], 2), ")")
  
  fit_ci <- c(foifit$foi, sigma_r, sigma_m, k)
  
  return(fit_ci)
  
}


# expected proportion seropositive due to past infection
pprev4_exp <- function(foi, sigma_r, sigma_m, M, age1, age2, mu_0, mu){
  
  if(age1<2 && age2<2){
    
    pp = mu_0 / (exp(-mu_0 * age1) - exp(-mu_0 * age2)) * (
      foi/(foi + sigma_r) * (
        (exp(-mu_0 * age1) - exp(-mu_0 * age2)) / mu_0 - 
          (exp(-(foi + sigma_r + mu_0)*age1) - exp(-(foi + sigma_r + mu_0)*age2))/(foi + sigma_r + mu_0)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu_0) * age1) - exp(-(sigma_m + mu_0) * age2)) / (sigma_m + mu_0) - 
                (exp(-(foi + sigma_r + mu_0) * age1) - exp(-(foi + sigma_r + mu_0) * age2)) / (foi + sigma_r + mu_0)
            )
          ) 
    )
    
  } else if(age1>=2 && age2>=2){
    
    pp = mu / (exp(-mu * age1) - exp(-mu * age2)) * (
      foi/(foi + sigma_r) * (
        (exp(-mu * age1) - exp(-mu * age2)) / mu - 
          (exp(-(foi + sigma_r + mu)*age1) - exp(-(foi + sigma_r + mu)*age2))/(foi + sigma_r + mu)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu) * age1) - exp(-(sigma_m + mu) * age2)) / (sigma_m + mu) - 
                (exp(-(foi + sigma_r + mu) * age1) - exp(-(foi + sigma_r + mu) * age2)) / (foi + sigma_r + mu)
            )
          ) 
    )
    
  } else if(age1<2 && age2>=2){
    
    F1 = mu * (exp(-mu_0 * age1) - exp(-2 * mu_0));
    
    F2 = mu_0 * exp(-2 * mu_0) * (1 - exp(-mu_0 * (age2 - 2)))
    
    pp1 = mu / (exp(-mu * age1) - exp(-mu * 2)) * (
      foi/(foi + sigma_r) * (
        (exp(-mu * age1) - exp(-mu * 2)) / mu - 
          (exp(-(foi + sigma_r + mu)*age1) - exp(-(foi + sigma_r + mu)*2))/(foi + sigma_r + mu)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu) * age1) - exp(-(sigma_m + mu) * 2)) / (sigma_m + mu) - 
                (exp(-(foi + sigma_r + mu) * age1) - exp(-(foi + sigma_r + mu) * 2)) / (foi + sigma_r + mu)
            )
          ) 
    )
    
    pp2 = mu / (exp(-mu * 2) - exp(-mu * age2)) * (
      foi/(foi + sigma_r) * (
        (exp(-mu * 2) - exp(-mu * age2)) / mu - 
          (exp(-(foi + sigma_r + mu)*2) - exp(-(foi + sigma_r + mu)*age2))/(foi + sigma_r + mu)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu) * 2) - exp(-(sigma_m + mu) * age2)) / (sigma_m + mu) - 
                (exp(-(foi + sigma_r + mu) * 2) - exp(-(foi + sigma_r + mu) * age2)) / (foi + sigma_r + mu)
            )
          ) 
    )
    
    pp = (F1 * pp1 + F2 * pp2) / (F1 + F2);
    
  }
  
  return(pp)
}

# expected proportion with maternal antibodies 
pmAbs_exp <- function(M, sigma_m, age2, age1, mu_0, mu){
  
  if(age1<2 && age2<2){
    
    mp = (((mu_0 * M) / (sigma_m + mu_0)) * (
      exp(-(mu_0 + sigma_m)* age1) - exp(-(mu_0 + sigma_m) * age2)))/ 
      (exp(-mu_0 * age1) - exp(-mu_0 * age2))
    
  } else if(age1>=2 && age2>=2){
    
    mp = (((mu * M) / (sigma_m + mu)) * (
      exp(-(mu + sigma_m)* age1) - exp(-(mu + sigma_m) * age2)))/ 
      (exp(-mu * age1) - exp(-mu * age2))
    
  } else if(age1<2 && age2>=2){
    
    F1 = mu * (exp(-mu_0 * age1) - exp(-2 * mu_0))
    
    F2 = mu_0 * exp(-2 * mu_0) * (1 - exp(-mu_0 * (age2 - 2)))
    
    mp1 = (((mu_0 * M) / (sigma_m + mu_0)) * (
      exp(-(mu_0 + sigma_m)* age1) - exp(-(mu_0 + sigma_m) * 2)))/ 
      (exp(-mu_0 * age1) - exp(-mu_0 * 2))
    
    mp2 = (((mu * M) / (sigma_m + mu)) * (
      exp(-(mu + sigma_m)* 2) - exp(-(mu + sigma_m) * age2)))/ 
      (exp(-mu * 2) - exp(-mu * age2))
    
    mp = (F1 * mp1 + F2 * mp2) / (F1 + F2)
    
  }
  
  return(mp)
}

# total expected observed seroprevalence
total_pprev4_exp <- function(foi, sigma_r, sigma_m, M, age1, age2, mu_0, mu, sens, spec){
  
  total_true <- pprev4_exp(foi, sigma_r, sigma_m, M, age1, age2, mu_0, mu) + 
    pmAbs_exp(M, sigma_m, age2, age1, mu_0, mu)
  total_obs <- sens * total_true + (1 - spec) * (1 - total_true)
}


# assuming overdispersion and exponential age distribution
# simulate data-sets
# assign prevalence of Abs and mAbs

sim_data_betabinom_exp <- function (n_datasets, n_ages, gamma, sigma, omega, mabs,
                                    N_camels, age_upper, age_lower, od, sens, spec, 
                                    mu_0, mu){
  
  M_initial <- vector(length = n_datasets)
  pred_prev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pred_mAb <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  pos_data <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  obs_pprev <- matrix(data = NA, nrow = n_datasets, ncol = n_ages, byrow = T)
  
  for(s in 1:n_datasets){
    
    if(mabs == 1){
      
      M_initial[s] <-total_pprev4_exp(foi = gamma[s],
                                      age2 = 10,
                                      age1 = 4,
                                      sigma_r = sigma,
                                      sigma_m = omega,
                                      M = 0,
                                      mu_0 = mu_0,
                                      mu = mu,
                                      sens = sens,
                                      spec = spec) 
      
    } else {
      
      M_initial[s] <- 0
    }
    
    for(a in 1:n_ages){
      pred_prev[s,a] <- pprev4_exp(foi = gamma[s],
                                   age2 = age_upper[s,a], 
                                   age1 = age_lower[s,a],
                                   sigma_r = sigma, 
                                   sigma_m = omega, 
                                   M = M_initial[s],
                                   mu_0 = mu_0,
                                   mu = mu
      )
      pred_mAb[s,a]<- pmAbs_exp(M = M_initial[s], 
                                sigma_m = omega, 
                                age2 = age_upper[s,a], 
                                age1 = age_lower[s,a],
                                mu_0 = mu_0,
                                mu = mu
      )
      obs_pprev[s,a] <- sens * (pred_prev[s,a] + pred_mAb[s,a]) + 
        (1 - spec) * (1 - (pred_prev[s,a] + pred_mAb[s,a]))
      if(N_camels[s,a] > 0){
        pos_data[s,a] <- rbetabinom(n = 1, 
                                    size = N_camels[s,a], 
                                    prob = obs_pprev[s,a],
                                    theta = (((N_camels[s,a]) - 1)/od) - 1)
      } else {
        pos_data[s,a] <- 0
      }
    }
  }
  
  return(list(M_initial = M_initial,
              pmAbs = pred_mAb,
              pAbs = pred_prev,
              obs_seroprev = obs_pprev,
              simulated = pos_data))
}

# plot fig 1A for manuscript
plot_fig1a <- function(fit, data){
  
  # data
  data <- data%>%
    dplyr::select(STUDY_COUNTRY, LOW_AGE, UPP_AGE, seroprevalence, ci_low, ci_upp)%>%
    rename(seroprev.c_data = seroprevalence,
           seroprev.u_data = ci_upp,
           seroprev.l_data = ci_low)
  
  REGION <- unique(data_sero$REGION_COUNTRY_STUDY)
  
  # model prediction
  predic <- fit[grep("seroprevalence", rownames(fit)),]
  predic$REGION <- rep(REGION, each = ncol(SEROPOS))
  
  # got to filter out the cushioning (stan requires matrix size consistency)
  predic$N <- as.vector(t(N_CAMELS))
  predic <- predic %>%
    dplyr::filter(N>0) %>%
    rename(seroprev.c_predic = mean,
           seroprev.u_predic = '97.5%',
           seroprev.l_predic = '2.5%')%>%
    dplyr::select(seroprev.c_predic,
                  seroprev.u_predic,
                  seroprev.l_predic,
                  REGION)
  rownames(predic) <- NULL
  
  # merge data and model
  data_predic <- cbind(data, predic)
  
  #convert to long form for plotting
  
  # step one seperate both bound and model vs data
  d_p_long <- pivot_longer(data_predic, 4:9, names_to = c("bound", "key"),
                           names_sep = "_",
                           values_to = "seroprevalence")
  
  # step two: need bound seperated back out for geom_rect/geom_seg
  d_p <- pivot_wider(d_p_long, names_from = bound,
                     values_from = seroprevalence)
  
  ## adding on fine scale age class predictions
  predic_fine <- fit[grep("fine_seroprev", rownames(fit)),]
  predic_fine$REGION <- rep(REGION, each = length(age1_fine))
  predic_fine <- predic_fine %>%
    rename(seroprev.c = mean,
           seroprev.u = '97.5%',
           seroprev.l = '2.5%')%>%
    mutate(LOW_AGE = rep(age1_fine, 23),
           UPP_AGE = rep(age2_fine, 23))%>%
    dplyr::select(seroprev.c,
                  seroprev.u,
                  seroprev.l,
                  LOW_AGE,
                  UPP_AGE,
                  REGION)
  rownames(predic) <- NULL
  
  # region labels for plot
  region.labs <- c("Egypt (18)", 
                   "Egypt (16)", 
                   "Ethiopia (19)",
                   "Kenya (20)",
                   "Kenya (21)",
                   "Kenya (22)",
                   "Senegal (16)",
                   "Tunisia (16)",
                   "Tunisia (19)",
                   "Uganda (16)",
                   "Iraq (16)",
                   "Iraq (23)",
                   "Jordan (16)",
                   "Jordan (24)",
                   "KSA, 1992 (25)",
                   "KSA, 2013 (25)",
                   "KSA (17)",
                   "KSA (26)",
                   "KSA (16)",
                   "UAE (11)",
                   "Bangladesh (27)",
                   "Pakistan (28)",
                   "Pakistan (29)")
  names(region.labs) <- sort(REGION)
  
  # plot 1
  p1 <- ggplot(d_p)+
    geom_segment(aes(x = LOW_AGE, xend = UPP_AGE,
                     y = seroprev.c, yend = seroprev.c,
                     col = key), lwd = 1)+
    geom_rect(aes(xmin  = LOW_AGE, xmax = UPP_AGE,
                  ymin = seroprev.l, ymax = seroprev.u,
                  fill = key), alpha = 0.2)+
    scale_colour_manual(values = c("red", "blue"), name = NULL,
                        labels = c("data", "predicted"))+
    scale_fill_manual(values = c("red", "blue"), name = NULL,
                      labels = c("data", "predicted"))+
    facet_wrap(~REGION, labeller = labeller(REGION = region.labs))+
    theme_minimal()+
    geom_line(data = predic_fine,
              aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                  y = seroprev.c), col = "blue")+
    geom_ribbon(data = predic_fine,
                aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                    ymin = seroprev.l, ymax = seroprev.u), fill = "blue",
                alpha = 0.2)+
    xlab("Age (years)")+
    ylab("seroprevalence (proportion)")+
    theme(text = element_text(size = 45))
  
  # plot 2 (alternative)
  p2 <- ggplot(d_p)+
    geom_segment(aes(x = LOW_AGE, xend = UPP_AGE,
                     y = seroprev.c, yend = seroprev.c,
                     col = key), alpha = 0.2, size = 2)+
    geom_errorbar(aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                      ymin = seroprev.l, ymax = seroprev.u,
                      colour = key), alpha = 0.8, lwd = 1)+
    geom_point(aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                   y = seroprev.c, col = key), alpha = 0.8, size = 2.5)+
    scale_colour_manual(values = c("red", "blue"), name = NULL,
                        labels = c("data", "predicted"))+
    scale_fill_manual(values = c("red", "blue"), name = NULL,
                      labels = c("data", "predicted"))+
    facet_wrap(~REGION, labeller = labeller(REGION = region.labs))+
    theme_minimal()+
    geom_line(data = predic_fine,
              aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                  y = seroprev.c), col = "blue")+
    geom_ribbon(data = predic_fine,
                aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                    ymin = seroprev.l, ymax = seroprev.u), fill = "blue",
                alpha = 0.1)+
    xlab(" ")+
    ylab(" ")+
    theme(text = element_text(size = 35),
          legend.position = c(0.8,0.1),
          panel.spacing.y = unit(2, "lines"))
  
  # plot 3 (alternative)
  p3 <- ggplot(d_p)+
    geom_point(aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                   y = seroprev.c,
                   col = key), alpha = 0.8)+
    geom_errorbar(aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                      ymin = seroprev.l, ymax = seroprev.u,
                      colour = key), alpha = 0.8, lwd = 1)+
    scale_colour_manual(values = c("red", "blue"), name = NULL,
                        labels = c("data", "predicted"))+
    scale_fill_manual(values = c("red", "blue"), name = NULL,
                      labels = c("data", "predicted"))+
    facet_wrap(~REGION, labeller = labeller(REGION = region.labs))+
    theme_minimal()+
    geom_line(data = predic_fine,
              aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                  y = seroprev.c), col = "blue")+
    geom_ribbon(data = predic_fine,
                aes(x = LOW_AGE + ((UPP_AGE - LOW_AGE)/2),
                    ymin = seroprev.l, ymax = seroprev.u), fill = "blue",
                alpha = 0.2)+
    xlab("Age (years)")+
    ylab("seroprevalence (proportion)")+
    theme(text = element_text(size = 25)) 
  
  print(p2)}

# plot figure 1B for the manuscript

plot_fig1b <- function(fit, central_stat){
  
  fit_foi <- fit[grep("foi", rownames(fit)),
                 c(central_stat, "2.5%", "25%", "75%", "97.5%")] %>%
    rename(lquart = "25%", uquart = "75%", low = "2.5%", high = "97.5%")
  names(fit_foi)[which(names(fit_foi) == central_stat)] <- "central_stat"
  fit_foi$REGION <- unique(data_sero$REGION_COUNTRY_STUDY)
  fit_foi <- fit_foi %>%
    mutate(group = case_when(grepl(REGION, pattern = "africa") ~ "africa",
                             grepl(REGION, pattern = "south_asia") ~ "south_asia",
                             grepl(REGION, pattern = "middle_east") ~ "middle_east",
                             TRUE ~ as.character(NA)))
  
  p1 <- ggplot(data = fit_foi)+
    geom_segment(aes(x = log(lquart), xend = log(uquart), y = REGION, yend = REGION, col = group), 
                 size = 3, alpha = 0.6)+
    geom_segment(aes(x = log(low), xend = log(high), y = REGION, yend = REGION, col = group), 
                 alpha = 0.6, size = 1)+
    geom_point(aes(x = log(central_stat), y = REGION, fill = group, col = group), 
               shape = 21, size = 3)+
    scale_colour_manual(values = c("springgreen3", "deepskyblue1", "maroon2"))+
    scale_fill_manual(values = c("springgreen3", "deepskyblue1", "maroon2"))+
    scale_y_discrete(limits = rev)+
    theme_classic()+
    ylab(" ")+
    xlab("log(Annual FoI)")+
    geom_vline(xintercept = 0, lty = 2, col = "grey")+
    theme(text = element_text(size = 30),
          axis.text.y = element_blank(),
          plot.margin = margin(1,1,1,1, unit = "cm"),
          legend.position = "none")
  
  # simplified alternative
  p2 <- ggplot(data = fit_foi)+
    geom_segment(aes(x = log(low), xend = log(high), y = REGION, yend = REGION, col = group), 
                 size = 5, alpha = 0.3)+
    geom_point(aes(x = log(central_stat), y = REGION, col = group, fill = group), 
               shape = 21, size = 4)+
    scale_colour_manual(values = c("springgreen3", "deepskyblue1", "maroon2"))+
    scale_fill_manual(values = c("springgreen3", "deepskyblue1", "maroon2"))+
    scale_y_discrete(limits = rev)+
    theme_classic()+
    ylab(" ")+
    xlab("log(Annual FoI)")+
    geom_vline(xintercept = 0, lty = 2, col = "grey")+
    theme(text = element_text(size = 45),
          axis.text.y = element_blank(),
          plot.margin = margin(1,1,1,1, unit = "cm"),
          legend.position = "none")
  
  print(p2)}

###################
### R0 analyses ###
###################

## interpolating the beta for the catalytic foi 
## converting foi to R0 with CIs

foi_to_R0_simple <- function(beta_vector, foi_vector, foi_cat, duration_infection){
  interpol_beta <- approx(x = foi_vector, y = beta_vector,
                          xout = foi_cat, method = "linear",
                          rule = 2)
  R0 <- interpol_beta$y * duration_infection
  return(list(R0 = R0,
              beta = interpol_beta$y))
}

R0_table <- function(beta_vector, foi_vector, foi_df, duration_infection){
  R0_tab <- data.frame(study = foi_df$study,
                       R0 = foi_to_R0_simple(beta_vector,
                                             foi_vector,
                                             foi_df$mode_ms,
                                             duration_infection)$R0,
                       R0_lower = foi_to_R0_simple(beta_vector,
                                                   foi_vector,
                                                   foi_df$`2.5`,
                                                   duration_infection)$R0,
                       R0_upper = foi_to_R0_simple(beta_vector,
                                                   foi_vector,
                                                   foi_df$`97.5`,
                                                   duration_infection)$R0)
  return(R0_tab)
}


# combine CrI and central values into one column for comparison table

R0_col <- function(R0){
  R0[2:4] <- round(R0[2:4], digits = 2)
  output <- paste(R0$R0," (", R0$R0_lower, ", ", R0$R0_upper, ")", sep = "")
  return(output)
}

###############
# periodicity #
###############

# function to format output - annual examples
# n = the number of variables you saved when you ran the model
format_out <- function(filepath, time_period, par_grid, n){
  out <- readRDS(filepath)
  idx <- readr::parse_number(filepath)
  out$t <- 0:time_period
  n_var = n*100
  out <- out %>%
    pivot_longer(cols = c(1:n_var),
                 names_to = c("var", "run_no"), names_sep = "\\.",
                 values_to = "values")
  
  out$seasonality <- par_grid$seasonality[idx]
  out$shedding <- par_grid$shedding[idx]
  out$R0 <- par_grid$R0[idx]
  
  return(out)
}

# function to format output - EXTRA-annual examples
# n = the number of variables you saved when you ran the model
format_out_extra <- function(filepath, time_period, par_grid, n){
  out <- readRDS(filepath)
  idx <- readr::parse_number(filepath)
  out$t <- 0:time_period
  n_var = n*100
  out <- out %>%
    pivot_longer(cols = c(1:n_var),
                 names_to = "run_no",
                 values_to = "values")
  
  out$seasonality <- par_grid$seasonality[idx]
  out$shedding <- par_grid$shedding[idx]
  out$R0 <- par_grid$R0[idx]
  
  return(out)
}
