# mers-cov-camels
![mers-camels-logo](https://github.com/user-attachments/assets/5a65014a-8a19-4f5b-be58-8644fc8205c3)

This repository contains the code used to model MERS-CoV transmission and vaccination in dromedary camels, accompanying the manuscript ***"Modelling transmission of MERS-CoV in dromedary camel populations and the potential impact of animal vaccination"*** (link pending).

# Summary of the research
Outbreaks of Middle East Respiratory syndrome coronavirus (MERS-CoV) in humans are driven by recurring zoonotic spillover from camels, leading to demand for camel vaccination. With two vaccine candidates shown to reduce infectiousness, there is a need to better understand transmission of MERS-CoV in camels and assess the potential impact of vaccination. To help address this, we used age-stratified seroprevalence data and a combination of modelling methodologies to estimate key epidemiological quantities including MERS-CoV transmissibility in camels and to estimate vaccine impact on infection incidence. Transmissibility was higher in the Middle East (R0 range 3-34) compared to Africa (2-15) and South Asia (2-4), highlighting the need for setting-specific vaccination strategies. Modelling suggested that even if the vaccine only reduced infectiousness rather than susceptibility to infection, vaccinating calves could achieve large reductions in incidence in moderate and high transmission settings, and interrupt transmission in low transmission settings, provided coverage was high (70-90%).

# Structure of the repository
## Data
- ***seroprevalence data*** used to fit catalytic models and estimate transmissibility of MERS-CoV in camels were collated by the systematic review [here](https://doi.org/10.1016/j.epidem.2019.100350) and supplemented by publicly available data from this later study by [Kandeil et al., 2019](https://doi.org/10.3390/v11080717) and are available in `data/real/sero_rna_syst.csv`.
- ***simulated data*** for catalytic model validation are available in `data/simulated/`.
- ***viral shedding data***: data extracted from [Haagmans et al., 2015](https://doi.org/10.1126/science.aad1283) and [Alharbi et al., 2019](https://doi.org/10.1038/s41598-019-52730-4) collected during vaccine trials in dromedaries and used to parameterise the infectiousness of reinfected and vaccinated animals is available in `data/shedding/`.

## Dependencies
- `dependencies.R` loads all the R libraries used within the code as well as sourcing `utils.R` which contains the functions used.

## Code
The analysis code consists of the catalytic models `catalytic-models-stan` (written in [stan](hamiltonian monte carlo algorithm](https://mc-stan.org/docs/2_24/cmdstan-guide/mcmc-config.html)), tranmsission models `dynamic-odin-models` (written using [odin](https://github.com/mrc-ide/odin/tree/a27f172ad11505c58353833f2e49905f34a0eec5)), and a series of analysis scripts used to run analysis tasks in R ordered from 1-11 in the order that you would want to use them taking into account their dependencies on one and other. At the top of each analysis script in R is a description of the scripts dependencies and its outputs.

### Catalytic models
- The catalytic models used to estimate the force of infection can be found in `catalytic-stan-models/`. The `likelihoods.stan` file includes the model functions used to predict seroprevalence and estimate log likelihood. The models are written and fitted in `rstan` which employs a [hamiltonian monte carlo algorithm](https://mc-stan.org/docs/2_24/cmdstan-guide/mcmc-config.html) to fit the models within a Bayesian framework.
  
### Dynamic transmission models
- The dynamic models of MERS-CoV transmission amongst camels used to estimate the reproduction number, R0, the critical community size, CCS, the periodicity of epidemics, and finally the potantial impact of vaccination, can be found in `dynamic-odin-models`. These include ths single patch and meta-population models in there original form and each extended to simulate vaccination (four models in total). The models are written using [odin](https://github.com/mrc-ide/odin/tree/a27f172ad11505c58353833f2e49905f34a0eec5) which provides a domain specific language that looks like R but is compiled in C allowing for rapid development and solving of models centred around ordinary differential equations. Due to the run time of the meta population models. Due to the nature of the meta-population models, [dust](https://mrc-ide.github.io/dust/articles/dust.html) was used to compile the odin models (interfaced in [odin-dust](https://mrc-ide.github.io/odin.dust/)) to allow them to run in parallel greatly shortening the run time (rather than relying in running them in a single thread as is the case if using R).

### Analysis scripts
- `1_process_raw_data.R`: this script processes the raw seroprevalence data from `data/sero_rna_syst.csv` and saves reshaped data used to fit the catalytic models and estimate force of infection
- `2_fit_real_data.R`: this script takes the processed data and fits the catalytic models to it.
- `3_process_fits.R`: this script takes the fits from the catalytic model and saves in a tidy way for using in tables and figures.Note that the raw fits (the HMC chains) are too large to be uploaded to github so you will need to rerun `1_process_raw_data.R` and `2_fit_real_data.R` on your local machine to access these, but the processed fits summarising the mean and credible intervals are available in `fits/processed_real/exp/`.
- `4_compare_models_tableS1S2.R`: this script takes the raw model fits from the catalytic model and compares model fit using DIC, generating Tables S1 and S2.
- `5_shedding_AUC.R`: this script takes the viral shedidng data from `data/shedding` and estimates relative infectiousness of reinfections and vaccinated animals to inform parameterisation of the dynamic transmission models.
- `6_estimate_R0_table1_figS4_tableS4.R`: this script takes our estimated force of infection, FoI, and estimates the reproduction number, R0 by calibrating the transmission rate beta to reflect the FoI estimated in the catalytic modelling. It generates Table 1 in the main text as well as the sensitivity analysis results in Table S4, and Figure S4 showing how we chose our beta transmission rate values to take forward into the remaining analysis in the supplementary material. 
- `7_simulate_CCS_sp.R`: this script estimates the critical community size (CCS) - the population size above which the extinction of MERS-CoV transmission by chance becomes unlikely, which we define as the population size needed for transmission to persist in <50% of stochastic based on the observation that the drop off from persistence to extinction across all model runs
happens very quickly over a relatively small change in population size.
- `8_simulate_CCS_mp.R`: this script estimates the CCS using the metapopulation model instead. This has a run time of approximately 24 hours.
- `9_simulate_periodicity.R`: this script estimates the periodicity of MERS-CoV infections in camels based on our assumption that births follow the same seasonality as observed in KSA for a range of R0 values reflecting our range of FoI estimates across the different camel populations studied. 
- `10_estimate_optimal_age_vax.R`
- `11_estimate_persistence_with_vax.R`
- `S1_create_sim_data.R`: this script simulates age stratified seroprevalence data used during model validation and problem solving.
- `S2_fit_sim_data.R`: this script takes the simulated data produced in `S1_creat_sim_data.R` and fits the catalytic models 1-4 to it.
- `S3_fit_real_data_LOO`: this script looks at the impact each dataset has on the overall estimates of force of infection and antibody waning, by leaving each out in turn. The output is used to produce FigS2 and S3 in `generate_figS2S3.R`.
- There are then 4 scripts that use outputs of the analysis scripts to generate the more complex figures for the main manuscript, and X scripts that generate the figures and tables for the supplementary material.

## Outputs
The summarised model fits and other forms of data presented in the manuscript are available in `fits/proccessed_real/`, `tables/` or `figs/`  but very large model fits could not be uploaded to github so the user will need to regenerate these themselves by rerunnning the models to populate the `fits/real_data` folders if they wish to use the raw fits. Please contact me (amydighe@gmail.com) if you wish to do this or you run into any issues with dependencies or errors. The figures and tables generated by the analysis scripts above are output into `figs/` and `tables/`.
