####################################################
## Meta-population Model 2024 - odin.dust version ##
####################################################

# This stochastic model simulates MERS-CoV transmission in a structured population,
# where the total population is split into a grid of 25 equally sized sub populations.
# sub populations can contribute to the foi in their neighbouring patches
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is approximated as 360 days 
# It is written using odin.dust

N_age <- user(49) #number of age classes
nr <- user(5, min = 3) # number of rows in grid of sub-pops
nc <- user(5, min = 3) # number of cols in grid of sub-pops

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
Ab_susc <- user() # proportion of susceptibility experienced if previously infected
reduced_shed <- user() # proportion of shedding/infectiousness seen in reinfections. default = no difference

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

# proportion infectious naive per patch
IN_patch[ , ] <- sum(I[ , i, j])/sum(N[ , i, j]) # where i = rows and j = cols

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
I2N_patch[ , ] <- sum(I2[ , i, j])/sum(N[ , i, j]) # where i = rows and j = cols

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

# background foi for introduction
foi_bg_usr <- user()
foi_bg <- if(tt < 3600) foi_bg_usr else 0

# balancing external and internal (within patch) foi
correction_ex[,] <- user()

# frequency dependent rate of infection
rate_internal_infection[ , ] <- beta *  (sum(I[ , i, j]) / sum(N[, i, j]))  + 
  beta * reduced_shed * (sum(I2[ , i, j]) / sum(N[ , i, j]))

rate_external_infection[ , ] <- beta * external_I1[i, j] + 
  beta * reduced_shed * external_I2[i, j]

rate_infection[ , ] <- (1 - correction_ex[i, j] * connectivity) * rate_internal_infection[i, j] + 
  connectivity * rate_external_infection[i, j] + foi_bg
  
rate_infection_mAb[ , ] <- mAb_susc * rate_infection[i, j]
rate_reinfection[ , ] <- Ab_susc * rate_infection[i, j]

#####################
## mortality rates ##
#####################
mu[] <- user() # user-defined age-dependent mortality rate

###########################
## recovery rate I --> R ## where R is non-infectious and completely immune to further infection
###########################
gamma <- user(1/14) # (/day) 

####################################################
## rate at which complete immunity wanes R --> S2 ## 
####################################################
sigma <- user() # (/day) to be taken from catalytic work eventually

###############################################################
## rate at which maternally-acquired immunity wanes Sm --> S ##
###############################################################
sigma_m <- user() # (/day) to be taken from catalytic work eventually

##############################################################################################################
# CONVERTING THESE RATES --> PROBABILITIES
# the above boil down to 7 rates which are converted to probabilities below
##############################################################################################################
p_infection[ , ] <- 1 - exp(-rate_infection[i, j])
p_infection_mAb[ , ] <- 1 - exp(-rate_infection_mAb[i, j])
p_reinfection[ , ] <- 1 - exp(-rate_reinfection[i, j])
p_mu[] <- 1 - exp(-mu[i]) # prob death
p_gamma <- 1 - exp(-gamma) # prob recovery
p_sigma <- 1 - exp(-sigma) # prob waned
p_sigma_m <- 1 - exp(-sigma_m) #prob mAbs waned

##############################################################################################################
# OUTFLOWS
# 1. infection, recovery, waning immunity and death
# 2. ageing
##############################################################################################################

#compartments are:
# Sm = susceptible but protected by mAbs
# S = susceptible
# I = infected and infectious
# R = recovered and completely immune
# S2 = immunity has waned to some degree
# I2 = infected and infectious for the 2nd+ time

# probability of leaving each compartment for any reason 
# (other than ageing which is dealt with later)
p_Sm[ , , ] <- 1 - exp(- (sigma_m + rate_infection_mAb[j, k] + mu[i])) 
p_S[ , , ] <- 1 - exp(- (rate_infection[j, k] + mu[i])) 
p_I[] <- 1 - exp(- (gamma + mu[i])) 
p_R[] <- 1 - exp(- (sigma + mu[i])) 
p_S2[ , , ] <- 1 - exp(- (rate_reinfection[j, k] + mu[i])) 
p_I2[] <- 1 - exp(- (gamma + mu[i])) 

# outflows due to infection, recovery or death
outflow_Sm[ , , ] <- rbinom(Sm[i, j, k], p_Sm[i, j, k])
outflow_S[ , , ] <- rbinom(S[i, j, k], p_S[i, j, k])
outflow_I[ , , ] <- rbinom(I[i, j, k], p_I[i])
outflow_R[ , , ] <- rbinom(R[i, j, k], p_R[i])
outflow_S2[ , , ] <- rbinom(S2[i, j, k], p_S2[i, j, k])
outflow_I2[ , , ] <- rbinom(I2[i, j, k], p_I2[i])


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
norm_p_sigma_m[ , , ] <- p_sigma_m/(p_sigma_m + p_infection_mAb[j, k] + p_mu[i])
norm_p_infection_mAb[ , , ] <- p_infection_mAb[j, k]/(p_sigma_m + p_infection_mAb[j, k] + p_mu[i])
norm_p_infection[ , , ] <- p_infection[j, k]/(p_infection[j, k] + p_mu[i])
norm_p_reinfection[ , , ] <- p_reinfection[j, k]/(p_reinfection[j, k] + p_mu[i])
norm_p_gamma[] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[] <- p_sigma/(p_sigma + p_mu[i])

# number of new infections, recoveries and newly susceptible
new_waned_mAb[ , , ] <- rbinom(outflow_Sm[i, j, k], norm_p_sigma_m[i, j, k])
new_infections_mAb[ , , ] <- rbinom(outflow_Sm[i, j, k], norm_p_infection_mAb[i, j, k])
new_infections[ , , ] <- rbinom(outflow_S[i, j, k], norm_p_infection[i, j, k])
new_recoveries[ , , ] <- rbinom(outflow_I[i, j, k], norm_p_gamma[i])
new_waned[ , , ] <- rbinom(outflow_R[i, j, k], norm_p_sigma[i])
new_reinfections[ , , ] <- rbinom(outflow_S2[i, j, k], norm_p_reinfection[i, j, k])
new_recoveries_two[ , , ] <- rbinom(outflow_I2[i, j, k], norm_p_gamma[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)
N_0 <- user(1000) # user-defined initial population size per patch
birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial patch population size
new_births[ , ] <- rpois(birth_rate) #per day per patch
births_protected[ , ] <- rbinom(new_births[i, j], seropoz_A4[i, j]) # protected by mAbs, patch specific
births_not_protected[ , ] <- new_births[i, j] - births_protected[i, j] #  NOT protected by mAbs, patch specific

#########################
## importation process ##
#########################

importation_rate <- user()
imported_cases <- rpois(importation_rate) #per day
imp_t <- user() # a user defined time at which cases are imported

###############################################################################################################
# EQUATIONS for movement of individuals between disease states and age classes
# time-step = 1 day
###############################################################################################################

## STEP 1 - disease state changes and births
## all importations (whether using rate or pulse) occur into age class 25 (~2 years old)

# no need for tt %% here as no inflows
new_Sm[1, , ] <- Sm[1, j, k] - outflow_Sm[1, j, k] + births_protected[j, k]
new_Sm[2:N_age, , ] <- Sm[i, j, k] - outflow_Sm[i, j, k]

new_S[1, , ] <- S[1, j, k] - outflow_S[1, j, k] + new_waned_mAb[1, j, k] + births_not_protected[j, k] 
new_S[2:N_age, , ] <- S[i, j, k] - outflow_S[i, j, k] + new_waned_mAb[i, j, k]

new_I[ , , ] <- I[i, j, k] - outflow_I[i, j, k] + new_infections[i, j, k] + new_infections_mAb[i, j, k]
new_I[25, nr -2, nc - 2] <- I[i, j, k] - outflow_I[i, j, k] + new_infections[i, j, k] + new_infections_mAb[i, j, k] + 
  imported_cases # overwrite for central patch

new_R[ , , ] <- R[i, j, k] - outflow_R[i, j, k] + new_recoveries[i, j, k] + new_recoveries_two[i, j, k]
new_S2[ , , ] <- S2[i, j, k] - outflow_S2[i, j, k] + new_waned[i, j, k]
new_I2[ , , ] <- I2[i, j, k] - outflow_I2[i, j, k] + new_reinfections[i, j, k]


## STEP 2 update with ageing

update(Sm[1, , ]) <- if(tt %% 30 == 0) 0 else new_Sm[1, j, k]
update(Sm[2:48, , ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j, k] else new_Sm[i, j, k]
update(Sm[N_age, , ]) <- if(tt %% 30 == 0) new_Sm[i - 1, j, k] + new_Sm[i, j, k] else new_Sm[i, j, k]

update(S[1, , ]) <- if(tt %% 30 == 0) 0 else new_S[1, j, k]
update(S[2:48, , ]) <- if(tt %% 30 == 0) new_S[i - 1, j, k] else new_S[i, j, k]
update(S[N_age, , ]) <- if(tt %% 30 == 0) new_S[i - 1, j, k] + new_S[i, j, k] else new_S[i, j, k]

update(I[1 , , ]) <-  if(tt %% 30 == 0) 0 else new_I[1, j, k]
update(I[2:48, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] else new_I[i, j, k]
update(I[25, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] else if(tt == imp_t) 1 + new_I[i, j, k] else new_I[i, j, k]
update(I[N_age, , ]) <- if(tt %% 30 == 0) new_I[i - 1, j, k] + new_I[i, j, k] else new_I[i, j, k]

update(R[1, , ]) <- if(tt %% 30 == 0) 0 else new_R[1, j, k]
update(R[2:48, , ]) <- if(tt %% 30 == 0) new_R[i - 1, j, k] else new_R[i, j, k]
update(R[N_age, , ]) <- if(tt %% 30 == 0) new_R[i - 1, j, k] + new_R[i, j, k] else new_R[i, j, k]

update(S2[1, , ]) <- if(tt %% 30 == 0) 0 else new_S2[1, j, k]
update(S2[2:48, , ]) <- if(tt %% 30 == 0) new_S2[i - 1, j, k] else new_S2[i, j, k]
update(S2[N_age, , ]) <- if(tt %% 30 == 0) new_S2[i - 1, j, k] + new_S2[i, j, k] else new_S2[i, j, k]

update(I2[1, , ]) <- if(tt %% 30 == 0) 0 else new_I2[1, j, k]
update(I2[2:48, , ]) <- if(tt %% 30 == 0) new_I2[i - 1, j, k] else new_I2[i, j, k]
update(I2[N_age, , ]) <- if(tt %% 30 == 0) new_I2[i - 1, j, k] + new_I2[i, j, k] else new_I2[i, j, k]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

update(Itot) <- sum(I[ , , ]) + sum(I2[ , , ])
update(Ntot) <- sum(Sm[, ,]) + sum(S[ ,  ,  ]) + sum(I[ ,  ,  ])+
  sum(R[ ,  ,  ]) + sum(S2[ ,  , ]) + sum(I2[ , , ])

## record total population size for use in FOI
N[ , , ] <- Sm[i, j, k] + S[i, j, k] + I[i, j, k] + R[i, j, k] + S2[i, j, k] + I2[i, j, k]
## record adult seroprevalence for use in births_protected
seropoz_A4[ , ] <- (S2[N_age, i, j] + I[N_age, i, j] + I2[N_age, i, j] + R[N_age, i, j])/ 
  (Sm[N_age, i , j] + S[N_age, i, j] + S2[N_age, i, j] + I[N_age, i, j] + I2[N_age, i, j] + R[N_age, i, j])

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

initial(tt) <- 1
initial(Itot) <- 0
initial(Ntot) <- sum(S_ini_p)

##################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- c(N_age, nr, nc)
dim(p_S) <- c(N_age, nr, nc)
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- c(N_age, nr, nc)
dim(p_I2) <- N_age

dim(Sm) <- c(N_age, nr, nc)
dim(S) <- c(N_age, nr, nc)
dim(I) <- c(N_age, nr, nc)
dim(R) <- c(N_age, nr, nc)
dim(S2) <- c(N_age, nr, nc)
dim(I2) <- c(N_age, nr, nc)

dim(S_ini_p) <- c(N_age, nr, nc)
dim(norm_p_infection) <- c(N_age, nr, nc)
dim(norm_p_infection_mAb) <- c(N_age, nr, nc)
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- c(N_age, nr, nc)
dim(norm_p_reinfection) <- c(N_age, nr, nc)

dim(new_waned_mAb) <- c(N_age, nr, nc)
dim(new_infections) <- c(N_age, nr, nc)
dim(new_infections_mAb) <- c(N_age, nr, nc)
dim(new_recoveries) <- c(N_age, nr, nc)
dim(new_waned) <- c(N_age, nr, nc)
dim(new_reinfections) <- c(N_age, nr, nc)
dim(new_recoveries_two) <- c(N_age, nr, nc)

dim(outflow_Sm) <- c(N_age, nr, nc)
dim(outflow_S) <- c(N_age, nr, nc)
dim(outflow_I) <- c(N_age, nr, nc)
dim(outflow_R) <- c(N_age, nr, nc)
dim(outflow_S2) <- c(N_age, nr, nc)
dim(outflow_I2) <- c(N_age, nr, nc)

dim(new_Sm) <- c(N_age, nr, nc)
dim(new_S) <- c(N_age, nr, nc)
dim(new_I) <- c(N_age, nr, nc)
dim(new_R) <- c(N_age, nr, nc)
dim(new_S2) <- c(N_age, nr, nc)
dim(new_I2) <- c(N_age, nr, nc)

dim(N) <- c(N_age, nr, nc)
dim(IN_patch) <- c(nr, nc)
dim(I2N_patch) <- c(nr, nc)
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
dim(p_infection) <- c(nr, nc)
dim(p_infection_mAb) <- c(nr, nc)
dim(p_reinfection) <- c(nr, nc)
dim(external_I1) <- c(nr, nc)
dim(external_I2) <- c(nr, nc)

