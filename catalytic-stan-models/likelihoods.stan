functions{
// assuming uniform age distribution within each age class

// predicted seroprevalence due to infection
   real prev4_int(real foi, 
                  real sigma_r, 
                  real sigma_m,
                  real M, 
                  real age1, 
                  real age2){
                    
    real pp;
  
    pp = 1 / (age2 - age1) * (
           foi/(foi + sigma_r) * (
              (age2 - age1) + 1 / (sigma_r + foi) * (
                exp(-(foi + sigma_r)*age2) - exp(-(foi + sigma_r) * age1)
              )
            ) -
        M * foi / (foi + sigma_r - sigma_m) * (
               1 / (foi + sigma_r) * (
                exp(-(sigma_r + foi) * age2) - exp(-(sigma_r + foi) * age1)
              ) -
               1 / sigma_m * (exp(-sigma_m * age2)- exp(-sigma_m * age1))
           )
        );
    
    return(pp);
}

//predicted total seroprevalence including mAbs and considering test sens/spec
    real seroprev(real foi, 
                  real sigma_r, 
                  real sigma_m,
                  real mabs, 
                  real age1, 
                  real age2,
                  real sens,
                  real spec){
                    
                    real pred_prev;
                    real pred_mab;
                    real pred_prev_tot;
                    real obs_pred_prev;
                    real m_initial_true;
                    real m_initial;
                    
              if(mabs > 0) {
       m_initial_true = prev4_int(foi, sigma_r, sigma_m, 0, 4.0000, 10.0000);
       m_initial = sens * m_initial_true + (1 - spec) * (1 - m_initial_true);
    } else {
       m_initial = 0;
    }
          pred_prev = prev4_int(foi, sigma_r, sigma_m, m_initial, age1, age2);
          pred_mab = 1/(age2 - age1)*((-m_initial/sigma_m)*(exp(-sigma_m*age2) - exp(-sigma_m*age1)));
          pred_prev_tot = pred_prev + pred_mab;  
          obs_pred_prev = sens * pred_prev_tot + (1 - spec) * (1 - pred_prev_tot);
    
    return(obs_pred_prev);
    }

// likeihood of model given data
    real model4av_bb_lpmf(int seropos,
                    int N,
                    real foi,
                    real age1,
                    real age2,
                    real sigma_r,
                    real sigma_m,
                    real k,
                    real sens, 
                    real spec,
                    real mabs
                    ){
                      
            real obs_pred_prev;
            real alpha;
            real beta;    
            real loglik;
    
    obs_pred_prev = seroprev(foi, sigma_r, sigma_m, mabs, age1, age2, sens, spec);
    alpha = (((N-1)/k) - 1)* obs_pred_prev;
    beta = (((N-1)/k) - 1)* (1 - obs_pred_prev);
    
    loglik = beta_binomial_lpmf(seropos|N, alpha, beta);
    
    return loglik;
    }
    
// assuming exponential age distribution within each age class

// predicted seroprevalence due to infection  
    real prev4_int_exp(real foi,
                        real sigma_r, 
                        real sigma_m,
                        real M, 
                        real age1, 
                        real age2,
                        real mu_0,
                        real mu){
                    
    real pp;
    
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
        );
          
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
        );
      
    } else if(age1<2 && age2>=2){
      
      real F1;
      real F2;
      real pp1;
      real pp2;
      
      F1 = mu * (exp(-mu_0 * age1) - exp(-2 * mu_0));
      
      F2 = mu_0 * exp(-2 * mu_0) * (1 - exp(-mu_0 * (age2 - 2)));
      
      pp1 = mu / (exp(-mu * age1) - exp(-mu * 2)) * (
        foi/(foi + sigma_r) * (
          (exp(-mu * age1) - exp(-mu * 2)) / mu - 
          (exp(-(foi + sigma_r + mu)*age1) - exp(-(foi + sigma_r + mu)*2))/(foi + sigma_r + mu)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu) * age1) - exp(-(sigma_m + mu) * 2)) / (sigma_m + mu) - 
              (exp(-(foi + sigma_r + mu) * age1) - exp(-(foi + sigma_r + mu) * 2)) / (foi + sigma_r + mu)
            )
          ) 
        );
      
      pp2 = mu / (exp(-mu * 2) - exp(-mu * age2)) * (
        foi/(foi + sigma_r) * (
          (exp(-mu * 2) - exp(-mu * age2)) / mu - 
          (exp(-(foi + sigma_r + mu)*2) - exp(-(foi + sigma_r + mu)*age2))/(foi + sigma_r + mu)) - (
            M * foi / (foi + sigma_r - sigma_m) * (
              (exp(-(sigma_m + mu) * 2) - exp(-(sigma_m + mu) * age2)) / (sigma_m + mu) - 
              (exp(-(foi + sigma_r + mu) * 2) - exp(-(foi + sigma_r + mu) * age2)) / (foi + sigma_r + mu)
            )
          ) 
        );
      
      pp = (F1 * pp1 + F2 * pp2) / (F1 + F2);
      
    }
    
    return(pp);
}

// predicted seroprevalence due to mAbs
    real prev_mAb_exp(real foi,
                        real sigma_r, 
                        real sigma_m,
                        real M, 
                        real age1, 
                        real age2,
                        real mu_0,
                        real mu){
                    
    real mp;
    
    if(age1<2 && age2<2){
      
      mp = (((mu_0 * M) / (sigma_m + mu_0)) * (
        exp(-(mu_0 + sigma_m)* age1) - exp(-(mu_0 + sigma_m) * age2)))
      / (exp(-mu_0 * age1) - exp(-mu_0 * age2));
          
    } else if(age1>=2 && age2>=2){
      
      mp = (((mu * M) / (sigma_m + mu)) * (
        exp(-(mu + sigma_m)* age1) - exp(-(mu + sigma_m) * age2)))
      / (exp(-mu * age1) - exp(-mu * age2));
      
    } else if(age1<2 && age2>=2){
      
      real F1;
      real F2;
      real mp1;
      real mp2;
      
      F1 = mu * (exp(-mu_0 * age1) - exp(-2 * mu_0));
      
      F2 = mu_0 * exp(-2 * mu_0) * (1 - exp(-mu_0 * (age2 - 2)));
      
      mp1 = (((mu_0 * M) / (sigma_m + mu_0)) * (
        exp(-(mu_0 + sigma_m)* age1) - exp(-(mu_0 + sigma_m) * 2)))
      / (exp(-mu_0 * age1) - exp(-mu_0 * 2));
      
      mp2 = (((mu * M) / (sigma_m + mu)) * (
        exp(-(mu + sigma_m)* 2) - exp(-(mu + sigma_m) * age2)))
      / (exp(-mu * 2) - exp(-mu * age2));
      
      mp = (F1 * mp1 + F2 * mp2) / (F1 + F2);
      
    }
    
    return(mp);
}

//predicted total seroprevalence including mAbs and considering test sens/spec    
    real seroprev_exp(real foi, 
                  real sigma_r, 
                  real sigma_m,
                  real mabs, 
                  real age1, 
                  real age2,
                  real sens,
                  real spec,
                  real mu_0,
                  real mu){
                    
                    real pred_prev;
                    real pred_mab;
                    real pred_prev_tot;
                    real obs_pred_prev;
                    real m_initial;
                    
              if(mabs > 0) {
       m_initial = prev4_int_exp(foi, sigma_r, sigma_m, 0, 4.0000, 10.0000, mu_0, mu);
    } else {
       m_initial = 0;
    }
          pred_prev = prev4_int_exp(foi, sigma_r, sigma_m, m_initial, age1, age2, mu_0, mu);
          pred_mab = prev_mAb_exp(foi, sigma_r, sigma_m, m_initial, age1, age2, mu_0, mu);
          pred_prev_tot = pred_prev + pred_mab;  
          obs_pred_prev = sens * pred_prev_tot + (1 - spec) * (1 - pred_prev_tot);
    
    return(obs_pred_prev);
    }

// likeihood of model given data 
    real model4exp_bb_lpmf(int seropos,
                    int N,
                    real foi,
                    real age1,
                    real age2,
                    real sigma_r,
                    real sigma_m,
                    real k,
                    real sens, 
                    real spec,
                    real mabs,
                    real mu_0,
                    real mu
                    ){
                      
            real obs_pred_prev;
            real alpha;
            real beta;    
            real loglik;
    
    obs_pred_prev = seroprev_exp(foi, sigma_r, sigma_m, mabs, age1, age2, sens, spec, mu_0, mu);
    alpha = (((N-1)/k) - 1)* obs_pred_prev;
    beta = (((N-1)/k) - 1)* (1 - obs_pred_prev);
    loglik = beta_binomial_lpmf(seropos|N, alpha, beta);
    
    return loglik;
    }
}

