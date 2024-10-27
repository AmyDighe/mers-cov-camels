####################################################
## Meta-population Model 2024 - odin.dust version ##
####################################################

# This stochastic model is an extension of metapop.R that simulates vaccination
# as well as MERS-CoV transmission in a structured population,
# where the total population is split into a grid of 25 equally sized sub populations.
# sub populations can contribute to the foi in their neighbouring patches
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is approximated as 360 days 
# It is written using odin.dust

N_age <- user(49) #number of age classes
nr <- user(min = 3) # number of rows in grid of sub-pops
nc <- user(min = 3) # number of cols in grid of sub-pops

##############################################################################################################
# RATES OF TRANSITION BETWEEN DISEASE STATES
# all rates are daily
# user-defined, defaults shown here in the parentheses
##############################################################################################################

################
## birth rate ##
################
alpha <- user()

#############################
## infection rates S --> I ##
#############################
beta <- user() #base rate

# adjusted beta values
mAb_susc <- user() # proportion of susceptibility experienced if maternal antibodies (mAbs) present
v_mAb_susc <- user() # proportion of susceptibility experienced if vaccinated AND maternal antibodies (mAbs) present
Ab_susc <- user() # proportion of susceptibility experienced if previously infected
v_susc <- user() # proportion of susceptibility experienced if vaccinated (and naive to natural infection)
v_Ab_susc <- user() # proportion of susceptibility experienced if vaccinated AND previously infected
reduced_shed <- user() # proportion of infectiousness seen in reinfections. default = no difference
v_reduced_shed <- user() # proportion of infectiousness seen in vaccinated first infections
v_shed <- user() # proportion of infectiousness seen in vaccinated reinfections

###############################
## meta-population structure ##
###############################

connectivity <- user() ##connection strength between patches (0<x<1)

## with edge effects

## 1 | 5 | 6 | 7 |2
## --+---+---+---+--
## 14| 17|18 |19 |8
## --+---+---+---+--
## 15| 20| 21| 22|9
## --+---+---+---+--
## 16|23 | 24| 25|10
## --+---+---+---+--
## 4 | 11| 12| 13|3

# contribution of infections in neighbouring patches

# patch population sizes (sum over age classes)
N_patch[ , ] <- sum(N[ , i, j])

# proportion infectious naive per patch
IN_patch[ , ] <- sum(I[ , i, j])/N_patch[i, j] # where i = rows and j = cols

# corners
external_I1[1,   1] <- IN_patch[i, j + 1] + IN_patch[i + 1, j]
external_I1[1,  nc] <- IN_patch[i, j - 1] + IN_patch[i + 1, j]
external_I1[nr,  1] <- IN_patch[i, j + 1] + IN_patch[i - 1, j]
external_I1[nr, nc] <- IN_patch[i, j - 1] + IN_patch[i - 1, j]

# edges
external_I1[2:(nr - 1),  1] <- IN_patch[i, j + 1] + IN_patch[i + 1, j] + IN_patch[i - 1, j]
external_I1[2:(nr - 1), nc] <- IN_patch[i, j - 1] + IN_patch[i + 1, j] + IN_patch[i - 1, j]
external_I1[1,  2:(nc - 1)] <- IN_patch[1, j - 1] + IN_patch[1, j + 1] + IN_patch[i + 1, j]
external_I1[nr, 2:(nc - 1)] <- IN_patch[1, j - 1] + IN_patch[1, j + 1] + IN_patch[i - 1, j]

# middles
external_I1[2:(nr - 1), 2:(nc - 1)] <-
  IN_patch[i, j + 1] + IN_patch[i, j - 1] + IN_patch[i + 1, j] + IN_patch[i - 1, j]

# same for reinfections
# proportion infectious reinfections per patch
I2N_patch[ , ] <- sum(I2[ , i, j])/N_patch[i, j] # where i = rows and j = cols

# corners
external_I2[1,   1] <- I2N_patch[i, j + 1] + I2N_patch[i + 1, j]
external_I2[1,  nc] <- I2N_patch[i, j - 1] + I2N_patch[i + 1, j]
external_I2[nr,  1] <- I2N_patch[i, j + 1] + I2N_patch[i - 1, j]
external_I2[nr, nc] <- I2N_patch[i, j - 1] + I2N_patch[i - 1, j]

# edges
external_I2[2:(nr - 1),  1] <- I2N_patch[i, j + 1] + I2N_patch[i + 1, j] + I2N_patch[i - 1, j]
external_I2[2:(nr - 1), nc] <- I2N_patch[i, j - 1] + I2N_patch[i + 1, j] + I2N_patch[i - 1, j]
external_I2[1,  2:(nc - 1)] <- I2N_patch[1, j - 1] + I2N_patch[1, j + 1] + I2N_patch[i + 1, j]
external_I2[nr, 2:(nc - 1)] <- I2N_patch[1, j - 1] + I2N_patch[1, j + 1] + I2N_patch[i - 1, j]

# middles
external_I2[2:(nr - 1), 2:(nc - 1)] <-
  I2N_patch[i, j + 1] + I2N_patch[i, j - 1] + I2N_patch[i + 1, j] + I2N_patch[i - 1, j]

## and now for vaccinated infectious individuals

# proportion infectious vaccinated per patch
vIN_patch[ , ] <- sum(vI[ , i, j])/N_patch[i, j] # where i = rows and j = cols

# corners
external_vI1[1,   1] <- vIN_patch[i, j + 1] + vIN_patch[i + 1, j]
external_vI1[1,  nc] <- vIN_patch[i, j - 1] + vIN_patch[i + 1, j]
external_vI1[nr,  1] <- vIN_patch[i, j + 1] + vIN_patch[i - 1, j]
external_vI1[nr, nc] <- vIN_patch[i, j - 1] + vIN_patch[i - 1, j]

# edges
external_vI1[2:(nr - 1),  1] <- vIN_patch[i, j + 1] + vIN_patch[i + 1, j] + vIN_patch[i - 1, j]
external_vI1[2:(nr - 1), nc] <- vIN_patch[i, j - 1] + vIN_patch[i + 1, j] + vIN_patch[i - 1, j]
external_vI1[1,  2:(nc - 1)] <- vIN_patch[1, j - 1] + vIN_patch[1, j + 1] + vIN_patch[i + 1, j]
external_vI1[nr, 2:(nc - 1)] <- vIN_patch[1, j - 1] + vIN_patch[1, j + 1] + vIN_patch[i - 1, j]

# middles
external_vI1[2:(nr - 1), 2:(nc - 1)] <-
  vIN_patch[i, j + 1] + vIN_patch[i, j - 1] + vIN_patch[i + 1, j] + vIN_patch[i - 1, j]

# same for vaccinated reinfections
# proportion infectious vaccinated reinfections per patch
vI2N_patch[ , ] <- sum(vI2[ , i, j])/N_patch[i, j] # where i = rows and j = cols

# corners
external_vI2[1,   1] <- vI2N_patch[i, j + 1] + vI2N_patch[i + 1, j]
external_vI2[1,  nc] <- vI2N_patch[i, j - 1] + vI2N_patch[i + 1, j]
external_vI2[nr,  1] <- vI2N_patch[i, j + 1] + vI2N_patch[i - 1, j]
external_vI2[nr, nc] <- vI2N_patch[i, j - 1] + vI2N_patch[i - 1, j]

# edges
external_vI2[2:(nr - 1),  1] <- vI2N_patch[i, j + 1] + vI2N_patch[i + 1, j] + vI2N_patch[i - 1, j]
external_vI2[2:(nr - 1), nc] <- vI2N_patch[i, j - 1] + vI2N_patch[i + 1, j] + vI2N_patch[i - 1, j]
external_vI2[1,  2:(nc - 1)] <- vI2N_patch[1, j - 1] + vI2N_patch[1, j + 1] + vI2N_patch[i + 1, j]
external_vI2[nr, 2:(nc - 1)] <- vI2N_patch[1, j - 1] + vI2N_patch[1, j + 1] + vI2N_patch[i - 1, j]

# middles
external_vI2[2:(nr - 1), 2:(nc - 1)] <-
  vI2N_patch[i, j + 1] + vI2N_patch[i, j - 1] + vI2N_patch[i + 1, j] + vI2N_patch[i - 1, j]

# background foi for introduction
foi_bg_usr <- user()
foi_bg <- if(tt < 3600) foi_bg_usr else 0

# balancing external and internal (within patch) foi
correction_ex[,] <- user()

# frequency dependent rate of infection
rate_internal_infection[ , ] <- beta *  (sum(I[ , i, j]) / N_patch[i, j])  + 
  beta * reduced_shed * (sum(I2[ , i, j]) / N_patch[i, j]) + 
  beta * v_shed * (sum(vI[ , i, j]) / N_patch[i, j]) + 
  beta * v_reduced_shed * (sum(vI2[ , i, j]) / N_patch[i, j])

rate_external_infection[ , ] <- beta * external_I1[i, j] + 
  beta * reduced_shed * external_I2[i, j] + 
  beta * v_shed * external_vI1[i, j] +
  beta * v_reduced_shed * external_vI2[i, j]

rate_infection[ , ] <- (1 - correction_ex[i, j] * connectivity) * rate_internal_infection[i, j] + 
  connectivity * rate_external_infection[i, j] + foi_bg

rate_infection_vaccinated[ , ] <- rate_infection[ i, j] * v_susc
rate_infection_mAb[ , ] <- rate_infection[i, j] * mAb_susc
rate_infection_mAb_vaccinated[ , ] <- rate_infection[i, j] * v_mAb_susc
rate_reinfection[ , ] <- rate_infection[i, j] * Ab_susc
rate_reinfection_vaccinated[ , ] <- rate_infection[i, j] * v_Ab_susc

#####################
## mortality rates ##
#####################
# user-defined age-dependent mortality rate
mu[] <- user()

###########################
## recovery rate I --> R ## where R is non-infectious and completely immune to further infection
###########################
gamma <- user(1/14) # (/day) 
v_gamma <- user(1/14)

####################################################
## rate at which complete immunity wanes R --> S2 ## 
####################################################
sigma <- user() # (/day) to be taken from catalytic work eventually
v_sigma <- user()

###############################################################
## rate at which maternally-acquired immunity wanes Sm --> S ##
###############################################################
sigma_m <- user() # (/day)
v_sigma_m <- user()

###########################
## proportion vaccinated ## (by age)
###########################
vaxp[] <- user()

##################################################
## rate at which vaccine induced immunity wanes ##
##################################################
rho <- user() #(/day) to be varied as currently unknown

##############################################################################################################
# CONVERTING THESE RATES --> PROBABILITIES
# the above boil down to 7 rates which are converted to probabilities below
##############################################################################################################
## for unvaccinated individuals
p_infection[ , ] <- 1 - exp(-rate_infection[i, j])
p_infection_mAb[ , ] <- 1 - exp(-rate_infection_mAb[i, j])
p_reinfection[ , ] <- 1 - exp(-rate_reinfection[i, j])
p_mu[] <- 1 - exp(-mu[i]) # prob death
p_gamma <- 1 - exp(-gamma) # prob recovery
p_sigma <- 1 - exp(-sigma) # prob waned
p_sigma_m <- 1 - exp(-sigma_m) #prob mAbs waned
## for vaccinated individuals
p_v_infection[ , ] <- 1 - exp(-rate_infection_vaccinated[i, j])
p_v_infection_mAb[ , ] <- 1 - exp(-rate_infection_mAb_vaccinated[i, j])
p_v_reinfection[ , ] <- 1 - exp(-rate_reinfection_vaccinated[i, j])
p_v_gamma <- 1 - exp(-v_gamma) # prob recovery if vaccinated
p_v_sigma <- 1 - exp(-v_sigma) # prob susceptible again if vaccinated
p_v_sigma_m <- 1 - exp(-v_sigma_m) #prob mAbs waned if vaccinated
p_rho <- 1 - exp(-rho) # prob vaccine induced immunity wanes

##############################################################################################################
# OUTFLOWS
# 1. infection, recovery, waning immunity and death
# 2. ageing
##############################################################################################################

#compartments are:
# Sm (and vSm for vaccinated individuals) = protected by mAbs
# S (vS) = susceptible
# I (vI) = infected and infectious
# R (vR) = recovered and completely immune
# S2 (vS2) = immunity has waned to some degree
# I2 (vI2) = infected and infectious for the 2nd+ time

# probability of leaving each compartment for any reason 
# (other than ageing which is dealt with later)
p_Sm[ , , ] <- 1 - exp(- (sigma_m + rate_infection_mAb[j, k] + mu[i])) 
p_S[ , , ] <- 1 - exp(- (rate_infection[j, k] + mu[i])) 
p_I[] <- 1 - exp(- (gamma + mu[i])) 
p_R[] <- 1 - exp(- (sigma + mu[i])) 
p_S2[ , , ] <- 1 - exp(- (rate_reinfection[j, k] + mu[i])) 
p_I2[] <- 1 - exp(- (gamma + mu[i]))
## for vaccinated individuals
p_vSm[ , , ] <- 1 - exp(- (v_sigma_m + rate_infection_mAb_vaccinated[j, k] + mu[i] + rho)) 
p_vS[ , , ] <- 1 - exp(- (rate_infection_vaccinated[ j, k] + mu[i] + rho)) 
p_vI[] <- 1 - exp(- (v_gamma + mu[i] + rho)) 
p_vR[] <- 1 - exp(- (v_sigma + mu[i] + rho)) 
p_vS2[ , , ] <- 1 - exp(- (rate_reinfection_vaccinated[j, k] + mu[i] + rho)) 
p_vI2[] <- 1 - exp(- (v_gamma + mu[i] + rho))

# outflows due to infection, recovery or death or vaccination
outflow_Sm[ , , ] <- rbinom(Sm[i, j, k], prob = p_Sm[i, j, k])
outflow_S[ , , ] <- rbinom(S[i, j, k], prob = p_S[i, j, k])
outflow_I[ , , ] <- rbinom(I[i, j, k], prob = p_I[i])
outflow_R[ , , ] <- rbinom(R[i, j, k], prob = p_R[i])
outflow_S2[ , , ] <- rbinom(S2[i, j, k], prob = p_S2[i, j, k])
outflow_I2[ , , ] <- rbinom(I2[i, j, k], prob = p_I2[i])
## for vaccinated individuals
outflow_vSm[ , , ] <- rbinom(vSm[i, j, k], prob = p_vSm[i, j, k])
outflow_vS[ , , ] <- rbinom(vS[i, j, k], prob = p_vS[i, j, k])
outflow_vI[ , , ] <- rbinom(vI[i, j, k], prob = p_vI[i])
outflow_vR[ , , ] <- rbinom(vR[i, j, k], prob = p_vR[i])
outflow_vS2[ , , ] <- rbinom(vS2[i, j, k], prob = p_vS2[i, j, k])
outflow_vI2[ , , ] <- rbinom(vI2[i, j, k], prob = p_vI2[i])

###############################################################################################################
# INFLOWS
# 1. where do these outflows go? 
# 2. births
# 3. importations
###############################################################################################################

##################################################
## new infections, recoveries and waning events ##
##################################################

#normalising the probabilities 
norm_p_sigma_m[, , ] <- p_sigma_m/(p_sigma_m + p_infection_mAb[j, k] + p_mu[i])
norm_p_infection_mAb[ , , ] <- p_infection_mAb[j, k]/(p_sigma_m + p_infection_mAb[j, k] + p_mu[i])
norm_p_infection[ , , ] <- p_infection[j, k]/(p_infection[j, k] + p_mu[i])
norm_p_reinfection[ , , ] <- p_reinfection[j, k]/(p_reinfection[j, k] + p_mu[i])
norm_p_gamma[] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[] <- p_sigma/(p_sigma + p_mu[i])

norm_p_v_sigma_m[ , , ] <- p_v_sigma_m/(p_v_sigma_m + p_v_infection_mAb[j, k] + p_mu[i] + p_rho)
norm_p_v_infection_mAb[ , , ] <- p_v_infection_mAb[j, k]/(p_v_sigma_m + p_v_infection_mAb[j, k] + p_mu[i] + p_rho)
norm_p_v_infection[ , , ] <- p_v_infection[j, k]/(p_v_infection[j, k] + p_mu[i] + p_rho)
norm_p_v_reinfection[ , , ] <- p_v_reinfection[j, k]/(p_v_reinfection[j, k] + p_mu[i] + p_rho)
norm_p_v_gamma[] <- p_v_gamma/(p_v_gamma + p_mu[i] + p_rho)
norm_p_v_sigma[] <- p_v_sigma/(p_v_sigma + p_mu[i] + p_rho)

norm_p_rho_vM[ , , ] <- p_rho/(p_v_sigma_m + p_v_infection_mAb[j, k] + p_mu[i] + p_rho)
norm_p_rho_vS[ , , ] <- p_rho/(p_v_infection[j, k] + p_mu[i] + p_rho)
norm_p_rho_vI[] <- p_rho/(p_v_gamma + p_mu[i] + p_rho)
norm_p_rho_vR[] <- p_rho/(p_v_sigma + p_mu[i] + p_rho)
norm_p_rho_vS2[ , , ] <- p_rho/(p_v_reinfection[j, k] + p_mu[i] + p_rho)
norm_p_rho_vI2[] <- p_rho/(p_v_gamma + p_mu[i] + p_rho)

# number of new infections, vaccinations, recoveries and newly susceptible
new_waned_mAb[ , , ] <- rbinom(outflow_Sm[i, j, k], prob = norm_p_sigma_m[i, j, k])
new_infections_mAb[ , , ] <- rbinom(outflow_Sm[i, j, k], prob = norm_p_infection_mAb[i, j, k])
new_infections[ , , ] <- rbinom(outflow_S[i, j, k], prob = norm_p_infection[i, j, k])
new_recoveries[ , , ] <- rbinom(outflow_I[i, j, k], prob = norm_p_gamma[i])
new_waned[ , , ] <- rbinom(outflow_R[i, j, k], prob = norm_p_sigma[i])
new_reinfections[ , , ] <- rbinom(outflow_S2[i, j, k], prob = norm_p_reinfection[i, j, k])
new_recoveries_two[ , , ] <- rbinom(outflow_I2[i, j, k], prob = norm_p_gamma[i])

v_new_waned_mAb[ , , ] <- rbinom(outflow_vSm[i, j, k], prob = norm_p_v_sigma_m[i, j, k])
v_new_infections_mAb[ , , ] <- rbinom(outflow_vSm[i, j, k], prob = norm_p_v_infection_mAb[i, j, k])
v_new_infections[ , , ] <- rbinom(outflow_vS[i, j, k], prob = norm_p_v_infection[i, j, k])
v_new_recoveries[ , , ] <- rbinom(outflow_vI[i, j, k], prob = norm_p_v_gamma[i])
v_new_waned[ , , ] <- rbinom(outflow_vR[i, j, k], prob = norm_p_v_sigma[i])
v_new_reinfections[ , , ] <- rbinom(outflow_vS2[i, j, k], prob = norm_p_v_reinfection[i, j, k])
v_new_recoveries_two[ , , ] <- rbinom(outflow_vI2[i, j, k], prob = norm_p_v_gamma[i])

new_waned_vax_M[ , , ] <- rbinom(outflow_vSm[i, j, k], prob = norm_p_rho_vM[i, j, k])
new_waned_vax_S[ , , ] <- rbinom(outflow_vS[i, j, k], prob = norm_p_rho_vS[i, j, k])
new_waned_vax_I[ , , ] <- rbinom(outflow_vI[i, j, k], prob = norm_p_rho_vI[i])
new_waned_vax_R[ , , ] <- rbinom(outflow_vR[i, j, k], prob = norm_p_rho_vR[i])
new_waned_vax_S2[ , , ] <- rbinom(outflow_vS2[i, j, k], prob = norm_p_rho_vS2[i, j, k])
new_waned_vax_I2[ , , ] <- rbinom(outflow_vI2[i, j, k], prob = norm_p_rho_vI2[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)
N_0 <- user() # user-defined initial population size per patch

birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial patch population size
new_births[ , ] <- rpois(birth_rate) #per day per patch
births_protected[ , ] <- rbinom(new_births[i, j], seropoz_A4[i, j]) # protected by mAbs, patch specific
births_not_protected[ , ] <- new_births[i, j] - births_protected[i, j] #  NOT protected by mAbs, patch specific

#########################
## importation process ##
#########################

importation_rate <- user(0)
imported_cases <- rpois(importation_rate) #per day
imp_t <- user() # a user defined time at which cases are imported

###############################################################################################################
# EQUATIONS for movement of individuals between disease states and age classes
# time-step = 1 day
###############################################################################################################

## STEP 1 - disease state changes
## all importations (whether using rate or pulse) occur into age class 25 (~2 years old)

# no need for tt %% here as no inflows
new_Sm[1, , ] <- Sm[1, j, k] - outflow_Sm[1, j, k] + births_protected[j,k] + new_waned_vax_M[1, j, k]
new_Sm[2:N_age, , ] <- Sm[i, j, k] - outflow_Sm[i, j, k] + new_waned_vax_M[i, j, k]

new_S[1, , ] <- S[1, j, k] - outflow_S[1, j, k] + new_waned_mAb[1, j, k] + births_not_protected[j, k] + new_waned_vax_S[1, j, k] 
new_S[2:N_age, , ] <- S[i, j, k] - outflow_S[i, j, k] + new_waned_mAb[i, j, k] + new_waned_vax_S[i, j, k]

new_I[ , , ] <- I[i, j, k] - outflow_I[i, j, k] + new_infections[i, j, k] + new_infections_mAb[i, j, k] + new_waned_vax_I[i, j, k]
new_I[25, nr -2, nc - 2] <- I[i, j, k] - outflow_I[i, j, k] + new_infections[i, j, k] + new_infections_mAb[i, j, k] + new_waned_vax_I[i, j, k] + imported_cases

new_R[ , , ] <- R[i, j, k] - outflow_R[i, j, k] + new_recoveries[i, j, k] + new_recoveries_two[i, j, k] + new_waned_vax_R[i, j, k]
new_S2[ , , ] <- S2[i, j, k] - outflow_S2[i, j, k] + new_waned[i, j, k] + new_waned_vax_S2[i, j, k]
new_I2[ , , ] <- I2[i, j, k] - outflow_I2[i, j, k] + new_reinfections[i, j, k] + new_waned_vax_I2[i, j, k]

# for vaccinated individuals
new_vSm[ , , ] <- vSm[i, j, k] - outflow_vSm[i, j, k]
new_vS[ , , ] <- vS[i, j, k] - outflow_vS[i, j, k] + v_new_waned_mAb[i, j, k]
new_vI[ , , ] <-  vI[i, j, k] - outflow_vI[i, j, k] + v_new_infections[i, j, k] + v_new_infections_mAb[i, j, k]
new_vR[ , , ] <- vR[i, j, k] - outflow_vR[i, j, k] + v_new_recoveries[i, j, k] + v_new_recoveries_two[i, j, k]
new_vS2[ , , ] <- vS2[i, j, k] - outflow_vS2[i, j, k] + v_new_waned[i, j, k]
new_vI2[ , , ] <- vI2[i, j, k] - outflow_vI2[i, j, k] + v_new_reinfections[i, j, k]


## STEP 2 update with ageing & vaccination

vax[] <- if(tt > (7199)) vaxp[i] else 0

vaxd_Sm[ , , ] <- rbinom(new_Sm[i, j, k], prob = vax[i])
vaxd_S[ , , ] <- rbinom(new_S[i, j, k], prob = vax[i])
vaxd_I[ , , ] <- rbinom(new_I[i, j, k], prob = vax[i])
vaxd_R[ , , ] <- rbinom(new_R[i, j, k], prob = vax[i])
vaxd_S2[ , , ] <- rbinom(new_S2[i, j, k], prob = vax[i])
vaxd_I2[ , , ] <- rbinom(new_I2[i, j, k], prob = vax[i])

update(Sm[1, , ]) <- if(tt %% 30 == 0) 0 else new_Sm[1, j, k]
update(Sm[2:48, , ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j, k] - vaxd_Sm[i - 1, j, k] else new_Sm[i, j, k]
update(Sm[N_age, , ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j, k] - vaxd_Sm[i - 1, j, k] + new_Sm[i, j, k] - vaxd_Sm[i, j, k] else new_Sm[i, j, k]

update(S[1, , ]) <- if(tt %% 30 == 0) 0 else new_S[1, j, k]
update(S[2:48, , ]) <- if(tt %% 30 == 0) new_S[i - 1, j, k] - vaxd_S[i - 1, j, k] else new_S[i, j, k]
update(S[N_age, , ]) <- if(tt %% 30 == 0) new_S[i - 1, j, k] - vaxd_S[i - 1, j, k] + new_S[i, j, k] - vaxd_S[i, j, k] else new_S[i, j, k]

update(I[1, , ]) <-  if(tt %% 30 == 0) 0 else new_I[1, j, k]
update(I[2:48, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] - vaxd_I[i - 1, j, k] else new_I[i, j, k]
update(I[25, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] - vaxd_I[i - 1, j, k] else if(tt == imp_t) 1 + new_I[i, j, k] else new_I[i, j, k]
update(I[N_age, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] - vaxd_I[i - 1, j, k] + new_I[i, j, k] - vaxd_I[i, j, k] else new_I[i, j, k]

update(R[1, , ]) <- if(tt %% 30 == 0) 0 else new_R[1, j, k]
update(R[2:48, , ]) <- if(tt %% 30 == 0) new_R[i - 1, j, k] - vaxd_R[i - 1, j, k] else new_R[i, j, k]
update(R[N_age, , ]) <- if(tt %% 30 == 0) new_R[i - 1, j, k] - vaxd_R[i - 1, j, k] + new_R[i, j, k] - vaxd_R[i, j, k] else new_R[i, j, k]

update(S2[1, , ]) <- if(tt %% 30 == 0) 0 else new_S2[1, j, k]
update(S2[2:48, , ]) <- if(tt %% 30 == 0) new_S2[i - 1, j, k] - vaxd_S2[i - 1, j, k] else new_S2[i, j, k]
update(S2[N_age, , ]) <- if(tt %% 30 == 0) new_S2[i - 1, j, k] - vaxd_S2[i - 1, j, k] + new_S2[i, j, k] - vaxd_S2[i, j, k] else new_S2[i, j, k]

update(I2[1, , ]) <- if(tt %% 30 == 0) 0 else new_I2[1, j, k]
update(I2[2:48, , ]) <- if(tt %% 30 == 0) new_I2[i - 1, j, k] - vaxd_I2[i - 1, j, k] else new_I2[i, j, k]
update(I2[N_age, , ]) <- if(tt %% 30 == 0) new_I2[i - 1, j, k] - vaxd_I2[i - 1, j, k] + new_I2[i, j, k] - vaxd_I2[i, j, k] else new_I2[i, j, k]

## and for vaccinated individuals
update(vSm[1, , ]) <- if(tt %% 30 == 0) 0 else new_vSm[1, j, k]
update(vSm[2:48, , ]) <- if(tt %% 30 == 0) vaxd_Sm[i - 1, j, k] + new_vSm[i - 1, j, k] else new_vSm[i, j, k]
update(vSm[N_age, , ]) <- if(tt %% 30 == 0) vaxd_Sm[i - 1, j, k] + new_vSm[i - 1, j, k] + vaxd_Sm[i, j, k] + new_vSm[i, j, k] else new_vSm[i, j, k]

update(vS[1, , ]) <- if(tt %% 30 == 0) 0 else new_vS[1, j, k]
update(vS[2:48, , ]) <- if(tt %% 30 == 0) vaxd_S[i - 1, j, k] + new_vS[i - 1, j, k] else new_vS[i, j, k]
update(vS[N_age, , ]) <- if(tt %% 30 == 0) vaxd_S[i - 1, j, k] + new_vS[i - 1, j, k] + vaxd_S[i, j, k] + new_vS[i, j, k] else new_vS[i, j, k]

update(vI[1, , ]) <-  if(tt %% 30 == 0) 0 else new_vI[1, j, k]
update(vI[2:48, , ]) <- if(tt %% 30 == 0) vaxd_I[i - 1, j, k] + new_vI[i - 1, j, k] else new_vI[i, j, k]
update(vI[N_age, , ]) <- if(tt %% 30 == 0) vaxd_I[i - 1, j, k] + new_vI[i - 1, j, k] + vaxd_I[i, j, k] + new_vI[i, j, k] else new_vI[i, j, k]

update(vR[1, , ]) <- if(tt %% 30 == 0) 0 else new_vR[1, j, k]
update(vR[2:48, , ]) <- if(tt %% 30 == 0) vaxd_R[i - 1, j, k] + new_vR[i - 1, j, k] else new_vR[i, j, k]
update(vR[N_age, , ]) <- if(tt %% 30 == 0) vaxd_R[i - 1, j, k] + new_vR[i - 1, j, k] + vaxd_R[i, j, k] + new_vR[i, j, k] else new_vR[i, j, k]

update(vS2[1, , ]) <- if(tt %% 30 == 0) 0 else new_vS2[1, j, k]
update(vS2[2:48, , ]) <- if(tt %% 30 == 0) vaxd_S2[i - 1, j, k] + new_vS2[i - 1, j, k] else new_vS2[i, j, k]
update(vS2[N_age, , ]) <- if(tt %% 30 == 0) vaxd_S2[i - 1, j, k] + new_vS2[i - 1, j, k] + vaxd_S2[i, j, k] + new_vS2[i, j, k] else new_vS2[i, j, k]

update(vI2[1, , ]) <- if(tt %% 30 == 0) 0 else new_vI2[1, j, k]
update(vI2[2:48, , ]) <- if(tt %% 30 == 0) vaxd_I2[i - 1, j, k] + new_vI2[i - 1, j, k] else new_vI2[i, j, k]
update(vI2[N_age, , ]) <- if(tt %% 30 == 0) vaxd_I2[i - 1, j, k] + new_vI2[i - 1, j, k] + vaxd_I2[i, j, k] + new_vI2[i, j, k] else new_vI2[i, j, k]

#update(seroprevalence[1:N_age]) <- (I[i] + R[i] + S2[i] + I2[i]) / (Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i])

update(seropoz_A4[ , ]) <- (sum(S2[N_age, i, j]) + sum(I[N_age, i, j]) + sum(I2[N_age, i, j]) + sum(R[N_age, i, j]) + 
                              sum(vS2[N_age, i, j]) + sum(vI[N_age, i, j]) + sum(vI2[N_age, i, j]) + sum(vR[N_age, i, j]))/ 
  (sum(Sm[N_age, i, j]) + sum(S[N_age, i, j]) + sum(S2[N_age, i, j]) + sum(I[N_age, i, j]) + sum(I2[N_age, i, j]) + sum(R[N_age, i, j]) + 
     sum(vSm[N_age, i, j]) + sum(vS[N_age, i, j]) + sum(vS2[N_age, i, j]) + sum(vI[N_age, i, j]) + sum(vI2[N_age, i, j]) + sum(vR[N_age, i, j]))

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

N[ , , ] <- Sm[i, j, k] + S[i, j, k] + I[i, j, k] + R[i, j, k] + S2[i, j, k] + I2[i, j, k] + 
  vSm[i, j, k] + vS[i, j, k] + vI[i, j, k] + vR[i, j, k] + vS2[i, j, k] + vI2[i, j, k]

update(Itot) <- sum(I[ , , ]) + sum(I2[ , , ]) + sum(vI[ , , ]) + sum(vI2[ , , ])
update(Ntot) <- sum(Sm[, ,]) + sum(S[ ,  ,  ]) + sum(I[ ,  ,  ])+
  sum(R[ ,  ,  ]) + sum(S2[ ,  , ]) + sum(I2[ , , ]) + 
  sum(vSm[, ,]) + sum(vS[ ,  ,  ]) + sum(vI[ ,  ,  ])+
  sum(vR[ ,  ,  ]) + sum(vS2[ ,  , ]) + sum(vI2[ , , ])
update(Itot_patch[,]) <- sum(I[, i, j]) + sum(I2[ , i, j]) + sum(vI[ , i, j]) + sum(vI2[ , i, j])

# incidence
update(incidence) <- sum(new_infections[ ,  ,]) + sum(new_reinfections[ , , ]) + sum(v_new_infections[ , , ]) + sum(v_new_reinfections[ , , ])
update(w_incidence) <- sum(new_infections[ , , ]) + reduced_shed * sum(new_reinfections[ , , ]) + 
  v_shed * sum(v_new_infections[ , , ]) + v_reduced_shed * sum(v_new_reinfections[ , , ])

update(new_inf) <- sum(new_infections[,,])
update(new_rinf) <- sum(new_reinfections[,,])
update(new_vinf) <- sum(v_new_infections[,,])
update(new_vrinf) <- sum(v_new_reinfections[,,])

# Re
update(Re) <- (beta / gamma) * ((sum(S[,,]) + reduced_shed * Ab_susc * sum(S2[,,]) + v_shed * Ab_susc * sum(vS[,,]) + v_reduced_shed * v_Ab_susc * sum(vS2[,,]))/ Ntot)

##################################################################################################################################
# initial conditions
# equilibrium solution approximated to speed up balancing of demography
# no equilibrium solution for infection at this point
##################################################################################################################################

## initial states

S_ini_p[ , , ] <- user()

initial(Sm[ , , ]) <- 0
initial(S[ , , ]) <- S_ini_p[i, j, k]
initial(I[ , , ]) <- 0
initial(R[ , , ]) <- 0
initial(S2[ , , ]) <- 0
initial(I2[ , , ]) <- 0
initial(vSm[ , , ]) <- 0
initial(vS[ , , ]) <- 0
initial(vI[ , , ]) <- 0
initial(vR[ , , ]) <- 0
initial(vS2[ , , ]) <- 0
initial(vI2[ , , ]) <- 0

initial(tt) <- 1
initial(Itot) <- 0
initial(Ntot) <- sum(S_ini_p)
initial(Itot_patch[,]) <- 0
initial(incidence) <- 0
initial(w_incidence) <- 0
initial(new_inf) <- 0
initial(new_rinf) <- 0
initial(new_vinf) <- 0
initial(new_vrinf) <- 0
initial(seropoz_A4[ , ]) <- 0
initial(Re) <- beta / gamma

##################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(vax) <- N_age
dim(vaxp) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- c(N_age, nr, nc)
dim(p_S) <- c(N_age, nr, nc)
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- c(N_age, nr, nc)
dim(p_I2) <- N_age
dim(p_vSm) <- c(N_age, nr, nc)
dim(p_vS) <- c(N_age, nr, nc)
dim(p_vI) <- N_age
dim(p_vR) <- N_age
dim(p_vS2) <- c(N_age, nr, nc)
dim(p_vI2) <- N_age


dim(Sm) <- c(N_age, nr, nc)
dim(S) <- c(N_age, nr, nc)
dim(I) <- c(N_age, nr, nc)
dim(R) <- c(N_age, nr, nc)
dim(S2) <- c(N_age, nr, nc)
dim(I2) <- c(N_age, nr, nc)
dim(vSm) <- c(N_age, nr, nc)
dim(vS) <- c(N_age, nr, nc)
dim(vI) <- c(N_age, nr, nc)
dim(vR) <- c(N_age, nr, nc)
dim(vS2) <- c(N_age, nr, nc)
dim(vI2) <- c(N_age, nr, nc)

dim(S_ini_p) <- c(N_age, nr, nc)
dim(norm_p_infection) <- c(N_age, nr, nc)
dim(norm_p_infection_mAb) <- c(N_age, nr, nc)
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- c(N_age, nr, nc)
dim(norm_p_reinfection) <- c(N_age, nr, nc)

dim(norm_p_v_infection) <- c(N_age, nr, nc)
dim(norm_p_v_infection_mAb) <- c(N_age, nr, nc)
dim(norm_p_v_gamma) <- N_age
dim(norm_p_v_sigma) <- N_age
dim(norm_p_v_sigma_m) <- c(N_age, nr, nc)
dim(norm_p_v_reinfection) <- c(N_age, nr, nc)
dim(norm_p_rho_vM) <- c(N_age, nr, nc)
dim(norm_p_rho_vS) <- c(N_age, nr, nc)
dim(norm_p_rho_vI) <- N_age
dim(norm_p_rho_vR) <- N_age
dim(norm_p_rho_vS2) <- c(N_age, nr, nc)
dim(norm_p_rho_vI2) <- N_age

dim(new_waned_mAb) <- c(N_age, nr, nc)
dim(new_infections) <- c(N_age, nr, nc)
dim(new_infections_mAb) <- c(N_age, nr, nc)
dim(new_recoveries) <- c(N_age, nr, nc)
dim(new_waned) <- c(N_age, nr, nc)
dim(new_reinfections) <- c(N_age, nr, nc)
dim(new_recoveries_two) <- c(N_age, nr, nc)

dim(v_new_waned_mAb) <- c(N_age, nr, nc)
dim(v_new_infections) <- c(N_age, nr, nc)
dim(v_new_infections_mAb) <- c(N_age, nr, nc)
dim(v_new_recoveries) <- c(N_age, nr, nc)
dim(v_new_waned) <- c(N_age, nr, nc)
dim(v_new_reinfections) <- c(N_age, nr, nc)
dim(v_new_recoveries_two) <- c(N_age, nr, nc)
dim(new_waned_vax_M) <- c(N_age, nr, nc)
dim(new_waned_vax_S) <- c(N_age, nr, nc)
dim(new_waned_vax_I) <- c(N_age, nr, nc)
dim(new_waned_vax_R) <- c(N_age, nr, nc)
dim(new_waned_vax_S2) <- c(N_age, nr, nc)
dim(new_waned_vax_I2) <- c(N_age, nr, nc)

dim(outflow_Sm) <- c(N_age, nr, nc)
dim(outflow_S) <- c(N_age, nr, nc)
dim(outflow_I) <- c(N_age, nr, nc)
dim(outflow_R) <- c(N_age, nr, nc)
dim(outflow_S2) <- c(N_age, nr, nc)
dim(outflow_I2) <- c(N_age, nr, nc)

dim(outflow_vSm) <- c(N_age, nr, nc)
dim(outflow_vS) <- c(N_age, nr, nc)
dim(outflow_vI) <- c(N_age, nr, nc)
dim(outflow_vR) <- c(N_age, nr, nc)
dim(outflow_vS2) <- c(N_age, nr, nc)
dim(outflow_vI2) <- c(N_age, nr, nc)

dim(new_Sm) <- c(N_age, nr, nc)
dim(new_S) <- c(N_age, nr, nc)
dim(new_I) <- c(N_age, nr, nc)
dim(new_R) <- c(N_age, nr, nc)
dim(new_S2) <- c(N_age, nr, nc)
dim(new_I2) <- c(N_age, nr, nc)

dim(new_vSm) <- c(N_age, nr, nc)
dim(new_vS) <- c(N_age, nr, nc)
dim(new_vI) <- c(N_age, nr, nc)
dim(new_vR) <- c(N_age, nr, nc)
dim(new_vS2) <- c(N_age, nr, nc)
dim(new_vI2) <- c(N_age, nr, nc)

dim(vaxd_Sm) <- c(N_age, nr, nc)
dim(vaxd_S) <- c(N_age, nr, nc)
dim(vaxd_I) <- c(N_age, nr, nc)
dim(vaxd_R) <- c(N_age, nr, nc)
dim(vaxd_S2) <- c(N_age, nr, nc)
dim(vaxd_I2) <- c(N_age, nr, nc)

dim(N) <- c(N_age, nr, nc)
dim(IN_patch) <- c(nr, nc)
dim(I2N_patch) <- c(nr, nc)
dim(vIN_patch) <- c(nr, nc)
dim(vI2N_patch) <- c(nr, nc)
dim(N_patch) <- c(nr, nc)
dim(correction_ex) <- c(nr, nc)

dim(new_births) <- c(nr, nc)
dim(births_protected) <- c(nr, nc)
dim(births_not_protected) <- c(nr, nc)
dim(seropoz_A4) <- c(nr, nc)
dim(rate_infection) <- c(nr, nc)
dim(rate_infection_mAb) <- c(nr, nc)
dim(rate_internal_infection) <- c(nr, nc)
dim(rate_external_infection) <- c(nr, nc)
dim(rate_reinfection) <- c(nr, nc)
dim(rate_infection_vaccinated) <- c(nr, nc)
dim(rate_reinfection_vaccinated) <- c(nr, nc)
dim(rate_infection_mAb_vaccinated) <- c(nr, nc)
dim(p_infection) <- c(nr, nc)
dim(p_infection_mAb) <- c(nr, nc)
dim(p_reinfection) <- c(nr, nc)
dim(p_v_infection) <- c(nr, nc)
dim(p_v_infection_mAb) <- c(nr, nc)
dim(p_v_reinfection) <- c(nr, nc)
dim(external_I1) <- c(nr, nc)
dim(external_I2) <- c(nr, nc)
dim(external_vI1) <- c(nr, nc)
dim(external_vI2) <- c(nr, nc)


dim(Itot_patch) <- c(nr, nc)

