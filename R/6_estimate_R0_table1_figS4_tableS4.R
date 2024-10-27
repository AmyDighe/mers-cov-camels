# This script takes our estimated force of infection, FoI, and estimates the reproduction number, R0 

# it depends on:
#  1. the single-patch transmission model ~dynamic-odin-models/single_patch_model.R
#  2. the estimates for FoI from our best fitting model ~fits/processed_real/exp/sens_spec_1/fit4bb.rds 
#     and as a sensitivity analysis, our second best fitting model ~fits/processed_real/exp/sens_spec_1/fit3bb.rds 
#  3. the function "R0_table", "foi_to_R0_simple" and "R0_col" in the source file ~/utils.R
#  4. packages: odin, dplyr & lemon, 

# it generates:
#  1. For each scenario[j], a vector of mean FoIs ~generated_data/mean_foi_scenario[j].rds 
#     estimates from the dynamic model "sir_model" to match up with those estimated in the 
#     catalytic modelling enabling calibration  of the transmission parameter beta. 
#     This vector is only used in this script.
#  1. Table 1 for the main text: ~tables/table1_R0.csv
#  2. Table S4 for the supplementary material: ~tables/tableS2_R0.csv
#  3. Figure S4 for the supplementary material: ~figs/FigS4.png

source("dependencies.R")

# load the transmission model
sir_model <- odin::odin("dynamic-odin-models/single_patch_model.R", verbose = FALSE, skip_cache = TRUE)

##########################
## customise parameters ##
##########################

# input an initial population size
N_0 <- 1e+7

# set a level of seasonality for births 
# (1 being strongly seasonal, 0 being not at all seasonal)
delta <-  1 

# input the time period that you wish to run the model for (in days)
time_period <- 50*360 + last_imp 
t <- seq(0:time_period)

# read in the FoI and rate of waning of mAbs from catalytic modelling best fitting model
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
omega_1 <- par_esti["sigma_m", "mode_ms"] / 360 # annual --> daily

# read in the FoI and rate of waning of mAbs from catalytic modelling SECOND BEST fitting model
par_esti_2 <- readRDS("fits/processed_real/exp/sens_spec_1/fit3bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
omega_2 <- par_esti_2["sigma_m", "mode_ms"] / 360 # annual --> daily


# convert annual FoI -> daily FoI to match model time step
foi_df <- par_esti[1:23,]/360
foi_range <- range(matrix(foi_df))
foi_df_2 <- par_esti_2[1:23,]/360
foi_range_2 <- range(matrix(foi_df_2))
foi_df$study <- unique(data_sero$STUDY_COUNTRY)
foi_df_2$study <- unique(data_sero$STUDY_COUNTRY)


# grid of pars for different immunity scenarios
par_grid_R0 <- data.frame(Ab_susc = c(rep(0.75, times = 3), rep(1, times = 2), rep(0.25, times = 2), rep(0.75, times = 6)), # susceptibility of previously infected animals
                          sigma = c(rep(1/30, times = 7), rep(1/90, times = 2), rep(1, times = 2), rep(1/30, times = 2)), # rate of waning of immunity
                          red_shed = c(1/90, 1/2, 1/4, rep(c(1/90, 1/2), times = 5)),# relative infectiousness of reinfections
                          omega = c(rep(omega_1, times = 11), rep(omega_2, times = 2))) # rate of waning of mAbs


# scenario names
scenario <- vector(length = (dim(par_grid_R0)[1]))
for(i in 1:(dim(par_grid_R0)[1])){
    scenario[i] <- paste("AbsSusc", par_grid_R0$Ab_susc[i], 
                         "DurCompImm=", 1/par_grid_R0$sigma[i], 
                        "RedShed=",  round(100*par_grid_R0$red_shed[i], 0),
                        "DurmAbs=", round(12/(360*par_grid_R0$omega[i]), 1),
                         sep = "_")
}

# Range of beta to try (transmission rate parameter) - calibrated to cover span of FoIs from catalytic modelling
beta_list <- list(seq(0.08, 3.6, by = 0.08),
                  seq(0.08, 0.4, by = 0.04),
                  seq(0.08, 0.8, by = 0.08),
                  seq(0.08, 3.2, by = 0.08),
                  seq(0.08, 0.4, by = 0.04),
                  seq(0.08, 3.6, by = 0.08),
                  seq(0.08, 0.8, by = 0.08),
                  seq(0.08, 3.6, by = 0.08),
                  seq(0.08, 0.8, by = 0.08),
                  seq(0.08, 3.6, by = 0.08),
                  seq(0.08, 0.4, by = 0.04),
                  seq(0.08, 3.6, by = 0.08),
                  seq(0.08, 0.4, by = 0.04))
names(beta_list) <- scenario



###############
## run model ##
###############

for(j in 1:(dim(par_grid_R0)[1])){
  
  Ab_susc <- par_grid_R0$Ab_susc[j]
  sigma <- par_grid_R0$sigma[j]
  reduced_shed <- par_grid_R0$red_shed[j]
  beta_vector <- beta_list[[j]]
  mean_foi <- vector(length = length(beta_vector))
  omega <- par_grid_R0$omega[j]
  
  for(i in 1:(length(beta_vector))){
    
   beta <- beta_vector[i]
    x <- sir_model$new(alpha = alpha, beta = beta, gamma = gamma, sigma = sigma, sigma_m = omega, Ab_susc = Ab_susc, 
                       mAb_susc = mAb_susc, reduced_shed = reduced_shed, mu = mu, N_0 = N_0,
                       importation_rate = importation_rate, imp_t = 1, delta = delta, ind1 = ind1, ind2 = ind2,
                       foi_bg_usr = foi_bg_usr)
    if(beta < 0.15){
      nruns <- 50
    } else{
      nruns <- 10
    }
    
    out <- as.data.frame(replicate(nruns, x$run(t)[, c(349, 352, 354, 356)]))
    out_I1 <- out[,grep("I_1", colnames(out))]
    out_I2 <- out[,grep("I_2", colnames(out))]
    out_Itot <- out[,grep("Itot", colnames(out))]
    out_N <- out[,grep("Ntot", colnames(out))]
    idx_persist <- which(out_Itot[last_imp + (50*360),] > 0)
    out_I1_persist <- out_I1[,idx_persist]
    out_I2_persist <- out_I2[,idx_persist]
    out_N_persist <- out_N[,idx_persist]
    
    #matplot(out_Itot, type = "l")
    #abline(v = c(last_imp + 25*360, last_imp + 35*360))
    
    if(length(idx_persist) == 0){
      mean_foi[i] <- NA
    } else if(length(idx_persist) > 1){
      I1_N <- colMeans(out_I1_persist[(last_imp + 40*360):(last_imp + 50*360),] / out_N_persist[(last_imp + 40*360):(last_imp + 50*360),])  
      I2_N <- colMeans(out_I2_persist[(last_imp + 40*360):(last_imp + 50*360),] / out_N_persist[(last_imp + 40*360):(last_imp + 50*360),])
      foi <- beta * I1_N + beta * reduced_shed * I2_N
      mean_foi[i] <- mean(foi)
    } else if(length(idx_persist) == 1) {
      I1_N <- mean(out_I1_persist[(last_imp + 40*360):(last_imp + 50*360)] / out_N_persist[(last_imp + 40*360):(last_imp + 50*360)])
      I2_N <- mean(out_I2_persist[(last_imp + 40*360):(last_imp + 50*360)] / out_N_persist[(last_imp + 40*360):(last_imp + 50*360)])
      foi <- beta * I1_N + beta * reduced_shed * I2_N
      mean_foi[i] <- mean(foi)
    }
    
    print(paste("i =", i, "j = ", j, "foi= ", mean_foi[i], sep = " ")) }
  
  saveRDS(file = paste("generated_data/mean_foi_", scenario[j], ".rds", sep = ""), object = mean_foi)
  
}

#########################
### pulling R0 values ###
#########################

################
### R0 table ###
################
R0_tab_model4 <- vector(mode = "list", length = (dim(par_grid_R0)[1]) - 2)
names(R0_tab_model4) <- scenario[1:((dim(par_grid_R0)[1]) - 2)]
for(j in 1:((dim(par_grid_R0)[1]) - 2)){
  R0_tab_model4[[j]] <- R0_table(beta_vector = beta_list[[j]],
                          foi_vector = mean_foi[[j]],
                          foi_df = foi_df,
                          duration_infection = 14) 
}

R0_full_model4 <- bind_rows(R0_tab_model4, .id = "column_label")

R0_tab_model3 <- vector(mode = "list", length = 2)
names(R0_tab_model3) <- scenario[12:13]
for(j in 12:13){
  R0_tab_model3[[j-11]] <- R0_table(beta_vector = beta_list[[j]],
                          foi_vector = mean_foi[[j]],
                          foi_df = foi_df_2,
                          duration_infection = 14) 
}
R0_full_model3 <- bind_rows(R0_tab_model3, .id = "column_label")

# bind together model 3 and model 4 estimates
R0_full <- bind_rows(R0_full_model4, R0_full_model3)

R0_appendix <- R0_full %>%
  mutate(R0_comp = paste(round(R0,1), " (", round(R0_lower,1), ", ",
                         round(R0_upper, 1), ")", sep = ""))
R0_appendix_wide <- R0_appendix %>% 
  dplyr::select(-R0, -R0_lower, -R0_upper)%>%
  pivot_wider(names_from = column_label, values_from = R0_comp)%>%
  mutate(region = unique(data_sero$REGION_COUNTRY_STUDY))

R0_appendix_order <- R0_appendix_wide[order(R0_appendix_wide$region),]

## subset core assumptions for main text
R0_mainTxt <- R0_appendix_order[, c(1:4)]

write.csv(file = "tables/table1_R0.csv", R0_mainTxt)

## for the full table we actually want to transpose it for ease of reading 
R0_app <- t(R0_appendix_order[,c(-1, -15)])
colnames(R0_app) <- R0_appendix_order$study
R0_app_p <- cbind(R0_app, format(round(par_grid_R0, 3), nsmalls = 2))
row.names(R0_app_p) <- NULL

## save full table for sensitivity analyses in the Supplementary Material
write.csv(file = "tables/tableS4_R0.csv", R0_app_p)

###################################################################################################
# R versus the transmission rate beta - plots to show selection of beta values for core scenarios #
###################################################################################################

beta_esti <- matrix(nrow = 23, ncol = (dim(par_grid_R0)[1]))
beta_esti_lower <- beta_esti
beta_esti_upper <- beta_esti
mean_foi <- vector(mode = "list", length = (dim(par_grid_R0)[1]))

for(i in 1:(dim(par_grid_R0)[1])){
  mean_foi[[i]] <- readRDS(file = paste("generated_data/mean_foi_", scenario[i], ".rds", sep = ""))
}

for(j in 1:(dim(par_grid_R0)[1])){
  
  beta_esti[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                     foi_vector = mean_foi[[j]],
                                     foi_cat = foi_df$mode_ms,
                                     duration_infection = 14)$beta
  beta_esti_lower[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                           foi_vector = mean_foi[[j]],
                                           foi_cat = foi_df$`2.5%`, 14)$beta
  beta_esti_upper[, j] <- foi_to_R0_simple(beta_vector = beta_list[[j]],
                                           foi_vector = mean_foi[[j]],
                                           foi_cat = foi_df$`97.5%`, 14)$beta
}

rownames(beta_esti)<- foi_df$study
rownames(beta_esti_lower)<- foi_df$study
rownames(beta_esti_upper)<- foi_df$study

dat <- data.frame(foi_df, R0_1 = 14 * beta_esti[,1],
                  R0_25 = 14 * beta_esti[,3],
                  R0_50 = 14 * beta_esti[,2])

cuts <- data.frame(R0_1 = c(3.5, 7, 14),
                   R0_25 = c(2.0, 3.5, 5.0), 
                   R0_50 = c(1.75, 2.3, 3.0),
                   y = rep(0.01, 3))
dat_long <- pivot_longer(dat, 5:7, names_to = "shedding", values_to = "R0")
dat_long$region <- rep(c("South Asia", "Africa", "Africa", "Africa",
                         "Middle East", "Middle East", "Middle East", "Middle East",
                         "Africa", "Africa", "Africa",
                         "Middle East", "Middle East", "Middle East", "Middle East", "Middle East",
                         "South Asia", "South Asia", "Africa", "Africa", "Africa", "Middle East", 
                         "Africa"), each = 3)

cuts_long <- pivot_longer(cuts, 1:3, names_to = "shedding", values_to = "R0")
cuts_long <- cuts_long %>%
  group_by(shedding)


shedding.labs <- c("1%", "50%")
names(shedding.labs) <- unique(cuts_long%>% filter(shedding != "R0_25")%>% pull(shedding))

p1 <- ggplot(dat_long%>% filter(shedding != "R0_25"))+
  geom_point(aes(x = R0, y = mode_ms, col = region), size = 3)+
  theme_minimal()+
  ylab("FOI")+
  geom_vline(aes(xintercept = R0), lty = 2, lwd = 1, col = "firebrick", 
             data = cuts_long%>% filter(shedding != "R0_25"))+
  geom_text(data = cuts_long%>%filter(shedding == "R0_1"), aes(x = R0 + 1.4, y = y, label = R0), 
            col = "firebrick", size = 3)+
  geom_text(data = cuts_long%>%filter(shedding == "R0_50"), aes(x = R0 + 0.2, y = y, label = R0), 
            col = "firebrick", size = 3)+
  lemon::facet_rep_grid(~shedding, scale = "free_x", labeller = labeller(shedding = shedding.labs))+
  scale_color_viridis_d()+
  theme(text = element_text(size = 21))+
  ggtitle("relative infectiousness of reinfections")+
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "bottom")+
  theme(axis.line = element_line(colour = "black"))


ggsave(filename= "figs/FigS4.png", p1, width = 10, height = 5, unit = "in")




