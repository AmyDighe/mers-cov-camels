#############################
## Single Patch Model 2024 ##
#############################

# This simple stochastic model simulates MERS-CoV transmission in a single, homogenously mixed, 
# age-stratified population of dromedary camels (with ageing)
# In this model 1 year is approximated as 360 days 
# It is written using Odin
# It can then be run from a user-edited R script (e.g. your own or ~R/6_estimate_R0.R or ~R/7_simulate_CCS.R)

N_age <- 49 #number of age classes

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

# background foi for introduction
foi_bg_usr <- user()
foi_bg <- if(tt < 3600) foi_bg_usr else 0 # corresponds to 5 new infections per 1000 susceptible individuals per year

# frequency dependent rate of infection
# when reduced_shed = 1 I and I2 are essentially the same compartment
rate_infection <- beta * (sum(I[1:N_age]) / sum(N[1:N_age]))  + 
  beta * reduced_shed * (sum(I2[1:N_age]) / sum(N[1:N_age])) +
  foi_bg
rate_infection_mAb <- mAb_susc * rate_infection
rate_reinfection <- Ab_susc * rate_infection

#####################
## mortality rates ##
#####################
# user-defined age-dependent mortality rate
mu[] <- user()

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
#p_alpha <- 1 - exp(-alpha) # prob birth
p_infection <- 1 - exp(-rate_infection)
p_infection_mAb <- 1 - exp(-rate_infection_mAb)
p_reinfection <- 1 - exp(-rate_reinfection)
p_mu[1:N_age] <- 1 - exp(-mu[i]) # prob death
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
p_Sm[1:N_age] <- 1 - exp(- (sigma_m + rate_infection_mAb + mu[i])) 
p_S[1:N_age] <- 1 - exp(- (rate_infection + mu[i])) 
p_I[1:N_age] <- 1 - exp(- (gamma + mu[i])) 
p_R[1:N_age] <- 1 - exp(- (sigma + mu[i])) 

p_S2[1:N_age] <- 1 - exp(- (rate_reinfection + mu[i])) 
p_I2[1:N_age] <- 1 - exp(- (gamma + mu[i])) 

# outflows due to infection, recovery or death
outflow_Sm[1:N_age] <- rbinom(Sm[i], p_Sm[i])
outflow_S[1:N_age] <- rbinom(S[i], p_S[i])
outflow_I[1:N_age] <- rbinom(I[i],  p_I[i])
outflow_R[1:N_age] <- rbinom(R[i],  p_R[i])
outflow_S2[1:N_age] <- rbinom(S2[i],  p_S2[i])
outflow_I2[1:N_age] <- rbinom(I2[i],  p_I2[i])


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
norm_p_sigma_m[1:N_age] <- p_sigma_m/(p_sigma_m + p_infection_mAb + p_mu[i])
norm_p_infection_mAb[1:N_age] <- p_infection_mAb/(p_sigma_m + p_infection_mAb + p_mu[i])
norm_p_infection[1:N_age] <- p_infection/(p_infection + p_mu[i])
norm_p_reinfection[1:N_age] <- p_reinfection/(p_reinfection + p_mu[i])
norm_p_gamma[1:N_age] <- p_gamma/(p_gamma + p_mu[i])
norm_p_sigma[1:N_age] <- p_sigma/(p_sigma + p_mu[i])

# number of new infections, recoveries and newly susceptible
new_waned_mAb[1:N_age] <- rbinom(outflow_Sm[i],  norm_p_sigma_m[i])
new_infections_mAb[1:N_age] <- rbinom(outflow_Sm[i],  norm_p_infection_mAb[i])
new_infections[1:N_age] <- rbinom(outflow_S[i],  norm_p_infection[i])
new_recoveries[1:N_age] <- rbinom(outflow_I[i],  norm_p_gamma[i])
new_waned[1:N_age] <- rbinom(outflow_R[i],  norm_p_sigma[i])
new_reinfections[1:N_age] <- rbinom(outflow_S2[i],  norm_p_reinfection[i])
new_recoveries_2[1:N_age] <- rbinom(outflow_I2[i],  norm_p_gamma[i])

###################
## birth process ##
###################

delta <- user() # modulates the seasonality of births (1 being strongly seasonal, 0 being not at all seasonal)
pi <- 3.14159 # odin doesn't have pi

# calculating a seasonal birthrate, with a one year periodicity, not too sharp peak. 
# can use alpha rather than p_alpha here (bc not coming from a finite pop)

birth_rate <- N_0 * alpha * (1 + (delta * (cos(2 * pi * tt / 360)))) # N_0 is the initial population size
new_births <- rpois(birth_rate) #per day

births_protected <- rbinom(new_births,  seropoz_A4) # protected by mAbs
births_not_protected <- new_births - births_protected #  NOT protected by mAbs

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
new_Sm[1] <- Sm[1] - outflow_Sm[1] + births_protected
new_Sm[2:N_age] <- Sm[i] - outflow_Sm[i]

new_S[1] <- S[1] - outflow_S[1] + new_waned_mAb[1] + births_not_protected 
new_S[2:N_age] <- S[i] - outflow_S[i] + new_waned_mAb[i]

new_I[1:24] <-  I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i]
new_I[25] <- I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i] + imported_cases
new_I[26:N_age] <-  I[i] - outflow_I[i] + new_infections[i] + new_infections_mAb[i]

new_R[1:N_age] <- R[i] - outflow_R[i] + new_recoveries[i] + new_recoveries_2[i]
new_S2[1:N_age] <- S2[i] - outflow_S2[i] + new_waned[i]
new_I2[1:N_age] <- I2[i] - outflow_I2[i] + new_reinfections[i]



## STEP 2 update with ageing

update(Sm[1]) <- if(tt %% 30 == 0) 0 else new_Sm[1]
update(Sm[2:48]) <- if(tt %% 30 == 0) new_Sm[i - 1] else new_Sm[i]
update(Sm[N_age]) <- if(tt %% 30 == 0) new_Sm[i - 1] + new_Sm[i] else new_Sm[i]

update(S[1]) <- if(tt %% 30 == 0) 0 else new_S[1]
update(S[2:48]) <- if(tt %% 30 == 0) new_S[i - 1] else new_S[i]
update(S[N_age]) <- if(tt %% 30 == 0) new_S[i - 1] + new_S[i] else new_S[i]

update(I[1]) <-  if(tt %% 30 == 0) 0 else new_I[1]
update(I[2:24]) <- if(tt %% 30 == 0) new_I[i - 1] else new_I[i]
update(I[25]) <- if(tt %% 30 == 0) new_I[i - 1] else if(tt == imp_t) 25 + new_I[i] else new_I[i]
update(I[26:48]) <- if(tt %% 30 == 0) new_I[i - 1] else new_I[i]
update(I[N_age]) <- if(tt %% 30 == 0) new_I[i - 1] + new_I[i] else new_I[i]

update(R[1]) <- if(tt %% 30 == 0) 0 else new_R[1]
update(R[2:48]) <- if(tt %% 30 == 0) new_R[i - 1] else new_R[i]
update(R[N_age]) <- if(tt %% 30 == 0) new_R[i - 1] + new_R[i] else new_R[i]

update(S2[1]) <- if(tt %% 30 == 0) 0 else new_S2[1]
update(S2[2:48]) <- if(tt %% 30 == 0) new_S2[i - 1] else new_S2[i]
update(S2[N_age]) <- if(tt %% 30 == 0) new_S2[i - 1] + new_S2[i] else new_S2[i]

update(I2[1]) <- if(tt %% 30 == 0) 0 else new_I2[1]
update(I2[2:48]) <- if(tt %% 30 == 0) new_I2[i - 1] else new_I2[i]
update(I2[N_age]) <- if(tt %% 30 == 0) new_I2[i - 1] + new_I2[i] else new_I2[i]

update(tt) <- tt + 1 # used to count time, must start at one for %% conditioning to work

## record total population size and seroprevalence in each age group, and in adults >4

N[] <- Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i]

update(seroprevalence[1:N_age]) <- (I[i] + R[i] + S2[i] + I2[i]) / (Sm[i] + S[i] + I[i] + R[i] + S2[i] + I2[i])

update(seropoz_A4) <- (sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]))/ (sum(Sm[N_age]) + sum(S[N_age]) + sum(S2[N_age]) + sum(I[N_age]) + sum(I2[N_age]) + sum(R[N_age]))

##################################################################################################################################
# initial conditions
# equilibrium solution approximated to speed up balancing of demography
# no equilibrium solution for infection at this point
##################################################################################################################################

## initial states

initial(Sm[1:N_age]) <- 0
initial(S[1:N_age]) <- S_ini_p[i]
initial(I[1:N_age]) <- 0
initial(R[1:N_age]) <- 0
initial(S2[1:N_age]) <- 0
initial(I2[1:N_age]) <- 0

initial(tt) <- 1
#initial(N[1:N_age]) <- S_ini_p[i]
initial(seroprevalence[1:N_age]) <- 0
initial(seropoz_A4) <- 0

## initial population size for use in birthrate

N_0 <- user(1000) # user-defined

## setting initial conditions using the equilibrium solution for age distribution

births_detr[1:360] <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * i / 360)))) # change to fixed N_0 = huge to avoid NaNs
births_det[1:360] <- rpois(births_detr[i]) 


## if we start the model with the equilibrium amount in each of the first 48 month-wide compartments,
## birthrate will be set to balance summed death rate of the equilibrium age distribution

a_max <- 32 ## estimated max number of years spent in the last age class before death (only 1% of the population above after 32 yrs)
## (calculated in "calculating_a_max.R using exp decay mod)
ind1[] <- user() # indexes of birth influx for demographic equilibrium solution
ind2[] <- user()

# setting initial number of camels at demographic equilibrium
S_ini[1] <- 0
S_ini[2:48] <- sum(births_det[ind1[i]:ind2[i]]) * exp(- (30 * sum(mu[1:(i - 1)])))

# special treatment for the final open-ended compartment 49
yearly_influx <- sum(births_det[1:360]) * exp(-(30 * (sum(mu[1:48])))) #annual inflow into the final age compartment
yr[1:a_max] <- i  # number of years of influx before max age expectancy reached
# camels remaining from each cohort to enter in the last 32 years influx that remain
cohort_remaining[1:a_max] <- yearly_influx * exp(- (360 * yr[i] * mu[N_age])) 
S_ini[N_age] <- sum(cohort_remaining[1:a_max])

# getting proportion of animals in each age-class, multiplying by N_0 and rounding to whole animals
S_ini_p[1:N_age] <- round((S_ini[i] / sum(S_ini[1:N_age])) * N_0)


################################################################################################################################

# OUTPUTS

################################################################################################################################

#########################################
## number of individuals in each state ##
#########################################

# total number of individuals still protected by maternal immunity
output(M) <- sum(Sm[1:N_age])

output(S_1) <- sum(S[1:N_age]) # susceptible individuals never infected
output(I_1) <- sum(I[1:N_age]) # individuals infectious for the 1st time
output(R_1) <- sum(R[1:N_age]) # individuals recovered from a 1st infection
output(S_2) <- sum(S2[1:N_age]) # susceptible individuals whose immunity has waned
output(I_2) <- sum(I2[1:N_age]) # individuals infectious for the 2nd+ time

output(N_pop[1:N_age]) <- N[i]

output(Stot) <- sum(S[1:N_age]) + sum(S2[1:N_age]) + sum(Sm[1:N_age]) # total number of susceptible individuals
output(Itot) <- sum(I[1:N_age]) + sum(I2[1:N_age]) # total number of infectious individuals
output(Rtot) <- sum(R[1:N_age]) 

output(Ntot) <- sum(N[1:N_age]) # total number of individuals
output(N_C) <- sum(Sm[1:12]) + sum(S[1:12]) + sum(S2[1:12]) + sum(I[1:12]) + sum(I2[1:12]) + sum(R[1:12])
output(N_J) <- sum(Sm[13:24]) + sum(S[13:24]) + sum(S2[13:24]) + sum(I[13:24]) + sum(I2[13:24]) + sum(R[13:24])
output(N_A) <- sum(Sm[25:N_age]) + sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age])

####################
## seroprevalence ##
####################
output(seropoz_A) <- 100 * (sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]))/(sum(Sm[25:N_age]) + sum(S[25:N_age]) + sum(S2[25:N_age]) + sum(I[25:N_age]) + sum(I2[25:N_age]) + sum(R[25:N_age]))
output(seropoz_J) <- 100 * (sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24])) / (sum(Sm[1:24]) + sum(S[1:24]) + sum(I[1:24]) + sum(R[1:24]) + sum(S2[1:24]) + sum(I2[1:24]))
output(seropz_tot) <- 100 * (sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]))/(sum(Sm[1:N_age]) + sum(S[1:N_age]) + sum(I[1:N_age]) + sum(R[1:N_age]) + sum(S2[1:N_age]) + sum(I2[1:N_age]))

############################
## age at first infection ##
############################
output(reinf_1) <- new_reinfections[1]
output(reinf_2) <- new_reinfections[2]
output(incidence_new_inf) <- sum(new_infections[1:N_age])
output(incidence_indig_inf) <- sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age])
output(total_incidence) <- (sum(new_infections[1:N_age]) + sum(new_reinfections[1:N_age]))

###########
## other ##
###########
output(inf) <- rate_infection
output(birthrate) <- birth_rate
output(births) <- new_births
output(importations) <- imported_cases
output(yy) <- yr[12]
output(Se) <- (sum(S[]) + Ab_susc * reduced_shed * sum(S2[]))/ sum(N[])
#output(foi) <- beta * sum(I[1:N_age]) / sum(N[1:N_age]) # foi applying to first infections only

################################################################################################################################

# assigning DIMENSIONS needed for arrays

################################################################################################################################## 

dim(mu) <- N_age
dim(p_mu) <- N_age
dim(p_Sm) <- N_age
dim(p_S) <- N_age
dim(p_I) <- N_age
dim(p_R) <- N_age
dim(p_S2) <- N_age
dim(p_I2) <- N_age

dim(Sm) <- N_age
dim(S) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(S2) <- N_age
dim(I2) <- N_age

dim(S_ini) <- N_age
dim(S_ini_p) <- N_age
dim(yr) <- a_max
dim(cohort_remaining) <- a_max
dim(norm_p_infection) <- N_age
dim(norm_p_infection_mAb) <- N_age
dim(norm_p_gamma) <- N_age
dim(norm_p_sigma) <- N_age
dim(norm_p_sigma_m) <- N_age
dim(norm_p_reinfection) <- N_age

dim(new_waned_mAb) <- N_age
dim(new_infections) <- N_age
dim(new_infections_mAb) <- N_age
dim(new_recoveries) <- N_age
dim(new_waned) <- N_age
dim(new_reinfections) <- N_age
dim(new_recoveries_2) <- N_age

dim(outflow_Sm) <- N_age
dim(outflow_S) <- N_age
dim(outflow_I) <- N_age
dim(outflow_R) <- N_age
dim(outflow_S2) <- N_age
dim(outflow_I2) <- N_age

dim(new_Sm) <- N_age
dim(new_S) <- N_age
dim(new_I) <- N_age
dim(new_R) <- N_age
dim(new_S2) <- N_age
dim(new_I2) <- N_age

dim(births_det) <- 360
dim(births_detr) <- 360
dim(seroprevalence) <- N_age
dim(N) <- N_age
dim(ind1) <- 48
dim(ind2) <- 48
dim(N_pop) <- N_age