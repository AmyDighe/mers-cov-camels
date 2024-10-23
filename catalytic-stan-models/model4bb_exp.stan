#include likelihoods.stan 
data{
  int S; //number of studies
  int A; //number of age classes
  int N[S,A]; //number of camels tested per age class per study
  int pos[S,A]; //number of seropositive camels per age class per study
  matrix[S,A] age1; //lower bound per age class per study
  matrix[S,A] age2; //upper bound per age class per study
  real sens[S]; // serotest sensitivity
  real spec[S]; // serotest specificity 
  real mabs; // mAbs assumed present or not
  real mu_0; // mortality rate in < 2 year olds
  real mu; // mortality rate in >=2 year olds
  // int B; // number of ages for seroprevalence (fine)
  // vector[B] age1_fine; // fine resolution ages for seroprevalence visualisation
  // vector[B] age2_fine; // fine resolution ages for seroprevalence visualisation
  
}

parameters{
  vector<lower = 0.0001, upper = 10>[S] foi; // force of infection parameter per study
  real<lower = 0.0001, upper = 1> sigma_r; // rate of waning Abs following infection
  real<lower = 0.0001, upper = 10> sigma_m; // rate of waning maternal Abs
  real<lower = 0.000001> k; // overdispersion
}

model{
  k ~ normal(0.000001,0.5); //prior for overdispersion
  for(s in 1:S){
    for(a in 1:A){
      if(N[s,a] != 0){
        target+= model4exp_bb_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, k, sens[s], spec[s], mabs, mu_0, mu);
      } else {
        target+= 0;
      }
    }
  }
}

generated quantities{
  matrix <lower = 0, upper = 100> [S,A] seroprevalence;
  matrix[S,A] log_lik;
 // matrix <lower = 0, upper = 100> [S,B] fine_seroprev;
  
  for(s in 1:S){
    for(a in 1:A){
      seroprevalence[s,a] = seroprev_exp(foi[s], sigma_r, sigma_m, mabs, age1[s,a], age2[s,a], sens[s], spec[s], mu_0, mu);
      if(N[s,a] != 0){
        log_lik[s,a] = model4exp_bb_lpmf(pos[s,a]| N[s,a], foi[s], age1[s,a], age2[s,a], sigma_r, sigma_m, k, sens[s], spec[s], mabs, mu_0, mu);
      } else {
        log_lik[s,a] = 0;
      }
    }

    // for(b in 1:B){
    //   fine_seroprev[s,b] = seroprev(foi[s], sigma_r, sigma_m, mabs, age1_fine[b], age2_fine[b], sens[s], spec[s]);
    // }
  }
}
