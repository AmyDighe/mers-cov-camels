# This script provides the starting conditions needed to initiate the dynamic model
# close to the zoographic equilibrium

# it does not depend on other scripts

# it generates:
#  1. the functions stoch_init (for the meta-population models in ~dynamic-odin-models/) 
#     and stoch_init_sp (for the single patch models in ~dynamic-odin-models/) used
#     to initiate the model close to equilibrium 

###############################################################
## stochastic initialisation at rough zoographic equilibrium ##
###############################################################
stoch_init <- function(alpha, delta, N_0, mu, N_age, n_r, n_c){
  
  set.seed(4)

births_detr <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * seq(1:360) / 360))))
births_det <- vector(length = 360) 

for(i in 1:360){
  births_det[i] <- rpois(n = 1, lambda = births_detr[i])
}

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

# setting initial number of camels at demographic equilibrium
S_ini <- vector(length = N_age)
S_ini_p <- vector(length = N_age)
S_ini_tot <- vector(length = N_age)

#### for the first 48 month-wide age classes
for(i in 1:(N_age-1)){
  S_ini[i] <- sum(births_det[ind1[i]:ind2[i]]) * exp(- (30 * sum(mu[1:(i - 1)])))
}

S_ini[1] <- 0

#### special treatment for the final open-ended compartment 49
## (a_max calculated in "calculating_a_max.R using exp decay mod)

## estimated max number of years spent in the last age class before death
a_max <- 32 

#annual inflow into the final age compartment
yearly_influx <- sum(births_det[]) * exp(-(30 * (sum(mu[1:48])))) 

# number of years of influx before max age expectancy reached
yr <- c(1:a_max)  

# camels remaining from each cohort to enter in the last 32 years that remain
cohort_remaining <- vector(length = a_max)
for(i in 1:a_max){
  cohort_remaining[i] <- yearly_influx * exp(- (360 * yr[i] * mu[N_age]))  
}
S_ini[N_age] <- sum(cohort_remaining[1:a_max])

# getting the mean expected number of animals in each age-class
S_ini_p <- vector(length = N_age)
for(i in 1:N_age){
  S_ini_p[i] <- N_0 * S_ini[i]/sum(S_ini)
}
S_ini_p <- rep(S_ini_p, n_c*n_r)

# sampling using this mean to get whole numbers of animals
S_ini_n <- vector(length = N_age*n_r*n_c)
for(i in 1:(n_r * n_c * N_age)){
  S_ini_n[i] <- rpois(n = 1, lambda = S_ini_p[i])
}
S_ini_n <- array(S_ini_n, dim = c(N_age, n_r, n_c))

return(S_ini_n)
}



###############################################################
## stochastic initialisation at rough zoographic equilibrium ##
###############################################################
stoch_init_sp <- function(alpha, delta, N_0, mu, N_age){
  
  set.seed(4)
  
  births_detr <- 10000000 * alpha * (1 + (delta * (cos(2 * pi * seq(1:360) / 360))))
  births_det <- vector(length = 360) 
  
  for(i in 1:360){
    births_det[i] <- rpois(n = 1, lambda = births_detr[i])
  }
  
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
  
  # setting initial number of camels at demographic equilibrium
  S_ini <- vector(length = N_age)
  S_ini_p <- vector(length = N_age)
  S_ini_tot <- vector(length = N_age)
  
  #### for the first 48 month-wide age classes
  for(i in 1:(N_age-1)){
    S_ini[i] <- sum(births_det[ind1[i]:ind2[i]]) * exp(- (30 * sum(mu[1:(i - 1)])))
  }
  
  S_ini[1] <- 0
  
  #### special treatment for the final open-ended compartment 49
  ## (a_max calculated in "calculating_a_max.R using exp decay mod)
  
  ## estimated max number of years spent in the last age class before death
  a_max <- 32 
  
  #annual inflow into the final age compartment
  yearly_influx <- sum(births_det[]) * exp(-(30 * (sum(mu[1:48])))) 
  
  # number of years of influx before max age expectancy reached
  yr <- c(1:a_max)  
  
  # camels remaining from each cohort to enter in the last 32 years that remain
  cohort_remaining <- vector(length = a_max)
  for(i in 1:a_max){
    cohort_remaining[i] <- yearly_influx * exp(- (360 * yr[i] * mu[N_age]))  
  }
  S_ini[N_age] <- sum(cohort_remaining[1:a_max])
  
  # getting the mean expected number of animals in each age-class
  S_ini_p <- vector(length = N_age)
  for(i in 1:N_age){
    S_ini_p[i] <- N_0 * S_ini[i]/sum(S_ini)
  }
  
  # sampling using this mean to get whole numbers of animals
  S_ini_n <- vector(length = N_age)
  for(i in 1:N_age){
    S_ini_n[i] <- rpois(n = 1, lambda = S_ini_p[i])
  }
  return(S_ini_n)
}
