# This script runs the meta-population transmission model to simulate vaccination
# of camels against MERS-CoV and estimates the coverage needed to interrupt transmission
# for different vaccine coverages and under different efficacy scenario and transmission intensities

# it depends on:
#  1. the transmission model "dynamic-odin-models/metapopvax.R"
#  2. the estimated rate of waning of mAbs from the best fitting catalytic model
#     "fits/processed_real/exp/sens_spec_1/fit4bb.rds
#  3. ~R/pars_vax.rds generated in 10_estimate_optimal_age_FigS6S7.R

# it generates:
#  1. Very large model output files summarising persistence of transmission
#     across multiple stochastic model runs for different vaccine coverages 
#     under different efficacy scenarios and transmission settings and 
#     saves them to ~generated_data/persistence - used to create tables 2, 3, S5 and S6


############################################################
## PERSISTENCE ANALYSIS - META-POP MODEL WITH VACCINATION ##
############################################################

###################
## compile model ##
###################
metavax <- odin.dust::odin_dust("dynamic-odin-models/metapopvax.R", verbose = FALSE, skip_cache = TRUE)

source("R/stochastic_init.R") # function to initiate model at zoographic equilibrium

#############################
## user defined variables: ##
#############################

## POPULATION STRUCTURE #######################################################

n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of sub-pops
N_patch <- n_r*n_c

# correction for balancing ex vs int foi
correction_ex <- matrix(c(2, rep(3, n_c - 2), 2,
                          rep(c(3, rep(4, n_c - 2), 3), n_r - 2),
                          2, rep(3, n_c - 2), 2), nrow = n_r, ncol = n_c)

## EPIDEMIOLOGICAL PARAMETERS ##################################################

# infection derived immunity parameters
waning <- 1/30 # average rate at which complete immunity wanes following infection (/day)
Ab_susc <- 0.75 # relative susceptibility following complete immunity waning
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
sigma_m <- par_esti["sigma_m", "mode_ms"] / 360 # read in mAbs waning from catalytic modelling

# set a level of seasonality for births 0-1
# (1 being strongly seasonal)
seasonality <-  1 

# infectiousness and susceptibility parameters under 3 vaccine efficacy scenarios, 
# 2 differing assumed relationships between RNA shedding and infectiousness 
# and 3 different transmission intensities 
pars_vax <- readRDS("R/pars_vax.rds")

# population size
patch_pop <- c(2000000/N_patch)

# other vaccine parameters
v_gamma <- gamma
v_mAb_susc <- mAb_susc
v_sigma <- waning
v_sigma_m <- sigma_m

# vaccination coverage vector
coverage <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
 

################################################
## initialise the model using core parameters ##
################################################

n_particles <- 100L

## stochastic initialization
S_ini_p <- stoch_init(alpha = alpha, delta = seasonality, 
                      N_0 = patch_pop[1],
                      mu = mu, 
                      N_age = N_age, 
                      n_r = n_r, n_c = n_c)

storage.mode(S_ini_p) <- "double"

x <- metavax$new(pars = list(N_age = N_age, nr = n_r, nc = n_c, 
                             alpha = alpha, 
                             beta = pars_vax$beta[1],
                             gamma = gamma, 
                             sigma = waning,
                             sigma_m = sigma_m, 
                             Ab_susc = Ab_susc,
                             mAb_susc = mAb_susc, 
                             reduced_shed = pars_vax$reduced_shed[1], 
                             mu = mu, 
                             N_0 = patch_pop[1], 
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
                             v_susc  = pars_vax$v_susc[1], 
                             v_mAb_susc = v_mAb_susc, 
                             v_Ab_susc = pars_vax$v_Ab_susc[1], 
                             v_shed = pars_vax$v_shed[1], 
                             v_reduced_shed = pars_vax$v_reduced_shed[1], 
                             vaxp = rep(coverage[1], N_age), 
                             rho = pars_vax$rho[1]
), 
step = 0, 
n_particles = n_particles, 
n_threads = 2L, 
seed = 1L)

x$set_index(c(x$info()$index$Itot, x$info()$index$incidence, x$info()$index$w_incidence, x$info()$index$Itot_patch, x$info()$index$Re))


##########################
## PERSISTENCE ANALYSIS ##
##########################
p_t <- 35*360 # time to assess persistence

# time steps to run the model 
# (in days)
steps <- seq(0, 20*360, by = 10)# up to the point of vaccination
steps_2 <- seq(20*360, p_t, by = 10) #after vaccination (different coverages)

yr_idx_l <- (seq(20*360, 34*360, by = 360) - 7200)/10
yr_idx_u <- yr_idx_l + 35
patch_years <- array(NA, dim = c(N_patch, length(yr_idx_l), length(coverage)))
an_Re <- array(NA, dim = c(dim(pars_vax)[1], length(yr_idx_l), length(coverage), 100))
an_inc <- an_Re
an_winc <- an_Re


persist <- matrix(NA, nrow = dim(pars_vax)[1], ncol = length(coverage))
inc <- persist
winc <- persist
Re <- persist

for(i in 1:length(patch_pop)){
## stochastic initialization
S_ini_p <- stoch_init(alpha = alpha, delta = seasonality, 
                      N_0 = patch_pop[i],
                      mu = mu, 
                      N_age = N_age, 
                      n_r = n_r, n_c = n_c)
storage.mode(S_ini_p) <- "double"

for(j in c(51:54)){
  #tictoc::tic() 
  x$update_state(pars = list(N_age = N_age, nr = n_r, nc = n_c, 
                             alpha = alpha, 
                             beta = pars_vax$beta[j],
                             gamma = gamma, 
                             sigma = waning,
                             sigma_m = sigma_m, 
                             Ab_susc = Ab_susc,
                             mAb_susc = mAb_susc, 
                             reduced_shed = pars_vax$reduced_shed[j], 
                             mu = mu, 
                             N_0 = patch_pop[i], 
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
                             v_susc  = pars_vax$v_susc[j], 
                             v_mAb_susc = v_mAb_susc, 
                             v_Ab_susc = pars_vax$v_Ab_susc[j], 
                             v_shed = pars_vax$v_shed[j], 
                             v_reduced_shed = pars_vax$v_reduced_shed[j], 
                             vaxp = rep(0, 49), 
                             rho = pars_vax$rho[j]), 
                 step = 0)
  out <- x$simulate(steps)
  saved_state <- x$state()
  
  
  for(k in 1:length(coverage)){
    x$update_state(pars = list(N_age = N_age, nr = n_r, nc = n_c, 
                               alpha = alpha, 
                               beta = pars_vax$beta[j],
                               gamma = gamma, 
                               sigma = waning,
                               sigma_m = sigma_m, 
                               Ab_susc = Ab_susc,
                               mAb_susc = mAb_susc, 
                               reduced_shed = pars_vax$reduced_shed[j], 
                               mu = mu, 
                               N_0 = patch_pop[i], 
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
                               v_susc  = pars_vax$v_susc[j], 
                               v_mAb_susc = v_mAb_susc, 
                               v_Ab_susc = pars_vax$v_Ab_susc[j], 
                               v_shed = pars_vax$v_shed[j], 
                               v_reduced_shed = pars_vax$v_reduced_shed[j], 
                               vaxp = c(rep(0, 5), coverage[k], rep(0, 43)), 
                               rho = pars_vax$rho[j]), 
                   step = 20*360,
                   state = saved_state,
                   set_initial_state = FALSE)
    
    state <- x$simulate(steps_2)
    Itot_mat <- (t(state[1,,]))
    persist[j,k] <- sum(Itot_mat[(p_t - (20*360))/10, ] > 0)
    
    Itot_patch <- state[4:28, , ]
    
    # collect annual measures
    for(m in 1:15){
      # how many runs did the patch experienced an infection, for each year after vaccination
      patch_years[ , m, k] <- rowSums(rowSums(Itot_patch[ , , (yr_idx_l[m] : yr_idx_u[m])], dim = 2) >0)
      # annual incidence each year following vaccine introduction
      an_inc[j,m,k,] <- rowSums(state[2,,(yr_idx_l[m] : yr_idx_u[m])])
      # annual weighted incidence each year following vaccine introduction
      # an_winc[j,m,k,] <- rowSums(state[3,,(yr_idx_l[m] : yr_idx_u[m])])
      # annual average Re for each run following vaccine introduction
      # an_Re[j,m,k,] <- rowMeans(state[29,,(yr_idx_l[m] : yr_idx_u[m])])
      }
    
    # mean cumulative incidence in the 10 years post vaccine introduction
    inc[j,k] <- mean(colSums(t(state[2, , 0:(((30*360) - (20*360))/10)])))
    # mean annual weighted incidence in the 10 years post vaccine introduction
    winc[j,k] <- mean(colSums(t(state[3, , 0:(((30*360) - (20*360))/10)])))
    
    print(paste0("i = ", i, ", j = ", j, ", k = ", k))
  }
  saveRDS(persist[j,], file = paste0("generated_data/persistence/av_persist_2mil", j, "_", i, ".rds"))
  saveRDS(inc[j,], file = paste0("generated_data/persistence/av_inc_2mil", j, "_", i, ".rds"))
  saveRDS(winc[j,], file = paste0("generated_data/persistence/av_winc_2mil", j, "_", i, ".rds"))
  saveRDS(patch_years, file = paste0("generated_data/persistence/patch_years_2mil", j, "_", i, ".rds"))
  saveRDS(an_inc[j,,,], file = paste0("generated_data/persistence/an_inc_2mil", j, "_", i, ".rds"))
  #tictoc::toc() 
} 
}