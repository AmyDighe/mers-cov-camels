# This script estimates the optimal age for vaccination to result in the greatest 
# reduction of incidence of MERS-CoV infection in camels

# it depends on:
#  1. the transmission model "dynamic-odin-models/single_patch_model_vax.R"
#  2. the estimated rate of waning of mAbs from the best fitting catalytic model
#     "fits/processed_real/exp/sens_spec_1/fit4bb.rds

# it generates:
#  1. Very large model output files summarising changes in incidence 
#     for different target ages under different efficacy scenarios and 
#     transmission settings and saves them to ~generated_data/
#  1. Figure 4 panel 1 in the main text (~figs/Fig4_panel1.png)
#  2. Figure S6 in the SM (~figs/FigS6.png)
#  3. FIgure S7 in the SM (~figs/FigS7.png)



##############################
# optimal age of vaccination #
##############################

# load single patch model 
sir_model_vax <- odin::odin("dynamic-odin-models/single_patch_model_vax.R", 
                            verbose = FALSE, skip_cache = TRUE)

# input the time period that you wish to run the model for (in days)
time_period <- 30*360 
t <- seq(0:time_period)

# specify parameters ##########################################################

# demography
seasonality <- 1 # relative strength of seasonality of births (1 = full strength, 0 = not seasonal)
N_0 <- 1000000

# user-input immunity parameters
waning <- 1/30 # average rate at which complete immunity wanes following infection (/day)
Ab_susc <- 0.75 # relative susceptibility following complete immunity waning
par_esti <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds") %>%
  dplyr::select(mode_ms, `2.5%`, `97.5%`)
sigma_m <- par_esti["sigma_m", "mode_ms"] / 360 # read in mAbs waning from catalytic modelling

# vaccination
age_targ_idx <- c(1:12, seq(14, 24, by = 2), seq(28, 48, by = 4))
v_gamma <- gamma # rate of recovery from infection in vaccinated animals
v_sigma <- waning # rate of waning complete infection induced immunity in vaccinated animals
v_sigma_m <- sigma_m # rate of waning of mAbs in vaccinated animals
v_mAb_susc <- mAb_susc # relative susceptibility of vaccinated animals with mAbs
coverage <- 0.8
rho <- c(1/360, 1/(3*360), 1/(10*360))
beta <- list("shedding_0.01" = c(3.5, 7, 14)/14,
             "shedding_0.25" = c(2.0, 3.5, 5.0)/14,
             "shedding_0.50" = c(1.75, 2.3, 3.0)/14)

# scenario 1: vaccine reduces infectiousness of naive and previously infected animals
# as well as reducing susceptibility to infections
pars_1 <- expand.grid(beta = beta[["shedding_0.01"]], rho = rho, 
                      v_shed = 0.01, v_reduced_shed = 0.0015, reduced_shed = 0.01,
                      v_susc = 0.75, v_Ab_susc = 0.50)
pars_50 <- expand.grid(beta = beta[["shedding_0.50"]], rho = rho, 
                       v_shed = 0.5, v_reduced_shed = 0.33, reduced_shed = 0.5,
                       v_susc = 0.75, v_Ab_susc = 0.50)
# scenario 2: core scenario - vaccine reduces infectiousness of naive and previously infected animals
# but does not reduce susceptibility to infections
pars_1_transblock <- expand.grid(beta = beta[["shedding_0.01"]], rho = rho, 
                                 v_shed = 0.01, v_reduced_shed = 0.0015, reduced_shed = 0.01,
                                 v_susc = 1, v_Ab_susc = 0.75)
pars_50_transblock <- expand.grid(beta = beta[["shedding_0.50"]], rho = rho, 
                                  v_shed = 0.5, v_reduced_shed = 0.33, reduced_shed = 0.5,
                                  v_susc = 1, v_Ab_susc = 0.75)
# scenario 3: vaccine only serves to boost naturally acquired immunity
# providing no effect in previously naive animals
pars_1_boostnat <- expand.grid(beta = beta[["shedding_0.01"]], rho = rho, 
                               v_shed = 1, v_reduced_shed = 0.0015, reduced_shed = 0.01,
                               v_susc = 1, v_Ab_susc = 0.50)
pars_50_boostnat <- expand.grid(beta = beta[["shedding_0.50"]], rho = rho, 
                                v_shed = 1, v_reduced_shed = 0.33, reduced_shed = 0.5,
                                v_susc = 1, v_Ab_susc = 0.50)

pars_vax <- rbind(pars_1, pars_50, pars_1_transblock, pars_1_boostnat,
                  pars_50_transblock, pars_50_boostnat)

# save so we can call again for persistence analysis 
saveRDS(pars_vax, "R/pars_vax.rds")

## run model ###################################################################

# set up vectors to store incidence
incidence <- matrix(nrow = dim(pars_vax)[1], ncol = 25)
incidence_l <- incidence
incidence_u <- incidence
# set up array to store incidence by state --. infectiousness-weighted incidence
#dim(age targeted, age, scenario)
age_w_inc <- array(dim = c(25, 49, dim(pars_vax)[1]))
age_w_inc_l <- age_w_inc
age_w_inc_u <- age_w_inc

for(i in c(1:(dim(pars_vax)[1]))){
  beta <- pars_vax$beta[i]
  
  for(j in 1:25){
    vaxp <- rep(0, 49)
    if(j <25){
      vaxp[age_targ_idx[j]] <- coverage
    }
    print(vaxp)
    # include any user-defined parameters as arguments here
    x <- sir_model_vax$new(alpha = alpha, beta = pars_vax$beta[i], gamma = gamma, sigma = waning, sigma_m = sigma_m, Ab_susc = Ab_susc, 
                           mAb_susc = mAb_susc, reduced_shed = pars_vax$reduced_shed[i], mu = mu, N_0 = N_0,
                           importation_rate = importation_rate, imp_t = imp_t, delta = seasonality, ind1 = ind1, ind2 = ind2,
                           v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_susc  = pars_vax$v_susc[i], v_mAb_susc = v_mAb_susc, 
                           v_Ab_susc = pars_vax$v_Ab_susc[i], v_shed = pars_vax$v_shed[i], v_reduced_shed = pars_vax$v_reduced_shed[i], 
                           vaxp = vaxp, rho = pars_vax$rho[i],
                           foi_bg_usr = foi_bg_usr)
    
    n_runs <- 200
    
    out <- as.data.frame(replicate(n_runs, x$run(t)[, c(613, 668:716)])) #save just weighted inc by age and total inc
    # output total incidence
    out_inc <- out[7200 : 10800, grep("total_incidence", colnames(out))]
    incidence[i, j] <- mean(colSums(out_inc[ , ]))
    incidence_l[i, j] <- quantile(colSums(out_inc[ , ]), 0.05)
    incidence_u[i, j] <- quantile(colSums(out_inc[ , ]), 0.95)
    # output the weighted incidence reflecting infectiousness, by age
    out_weighted_inc <- matrix(colSums(out[7200 : 10800, grep("weighted", colnames(out))]), 
                               nrow = n_runs, ncol = 49, byrow = T)
    
    age_w_inc[j, , i] <- colMeans(out_weighted_inc[,])
    age_w_inc_l[j, , i] <- apply(out_weighted_inc, 2, quantile, probs = 0.05, na.rm = TRUE)
    age_w_inc_u[j, , i] <- apply(out_weighted_inc, 2, quantile, probs = 0.95, na.rm = TRUE)
    
    print(paste(i,j, incidence[i,j], sep = " "))
  } #end of j loop
  saveRDS(incidence[i, ], file = paste("generated_data/vax_age_opti/inc/", "incidence_", i, ".rds", sep = ""))
  saveRDS(age_w_inc[, ,i], file = paste("generated_data/vax_age_opti/w_inc/", "w_incidence_", i, ".rds", sep = ""))
  saveRDS(incidence_l[i, ], file = paste("generated_data/vax_age_opti/inc/", "incidence_l_", i, ".rds", sep = ""))
  saveRDS(age_w_inc_l[, ,i], file = paste("generated_data/vax_age_opti/w_inc/", "w_incidence_l_", i, ".rds", sep = ""))
  saveRDS(incidence_u[i, ], file = paste("generated_data/vax_age_opti/inc/", "incidence_u_", i, ".rds", sep = ""))
  saveRDS(age_w_inc_u[, ,i], file = paste("generated_data/vax_age_opti/w_inc/", "w_incidence_u_", i, ".rds", sep = ""))
} #end of i loop
saveRDS(incidence[,], file = "generated_data/vax_age_opti/inc/incidence_s2_shed1.rds")
saveRDS(age_w_inc[,], file = "generated_data/vax_age_opti/output/weighted_incidence_s2_shed1.rds")
saveRDS(incidence_l[,], file = "generated_data/vax_age_opti/inc/incidence_s2_shed1_l.rds")
saveRDS(age_w_inc_l[,], file = "generated_data/vax_age_opti/output/weighted_incidence_s2_shed1_l.rds")
saveRDS(incidence_u[,], file = "generated_data/vax_age_opti/inc/incidence_s2_shed1_u.rds")
saveRDS(age_w_inc_u[,], file = "generated_data/vax_age_opti/output/weighted_incidence_s2_shed1_u.rds")


# process results

#######################
## extract incidence ##
#######################
n_scenario <- dim(pars_vax)[1]
inc_mat <- matrix(NA, nrow = n_scenario, ncol = 25)
inc_mat_l <- inc_mat
inc_mat_u <- inc_mat

for(i in 1:54){
  inc_mat[i,] <- readRDS(paste("generated_data/vax_age_opti/inc/", 
                               "incidence_", 
                               i,
                               ".rds", sep = ""))
  inc_mat_l[i,] <- readRDS(paste("generated_data/vax_age_opti/inc/",
                                 "incidence_l_",
                                 i,
                                 ".rds", sep = ""))
  inc_mat_u[i,] <- readRDS(paste("generated_data/vax_age_opti/inc/",
                                 "incidence_u_",
                                 i,
                                 ".rds", sep = ""))
}

scenario <- rep(c( "base", "base", 
                   "transblock","boostnat", 
                   "transblock", "boostnat"), each = 9)
shed <- rep(c("1%", "50%", "1%", "1%", "50%", "50%"), each = 9)
beta_rho <- vector(length = n_scenario)
for(i in 1:n_scenario){
  beta_rho[i] <- paste("R0", " = ", format(round(14 * pars_vax$beta[i], 1), nsmall = 1), 
                       ", \U03C1", " = ", 1/(360*pars_vax$rho[i]), " yrs", sep = "")}
beta_rho_qual <- paste0("R0 = ", rep(c("low", "medium", "high"), 18),
                        ", ", "\U03C1", " = ", 1/(360*pars_vax$rho))
# central values
inc_df <- as.data.frame(inc_mat)
names(inc_df) <- c(age_targ_idx, 49)
inc_df <- cbind(inc_df, pars_vax)
inc_df$scenario <- scenario
inc_df$beta_rho <- beta_rho
inc_df$beta_rho_qual <- beta_rho_qual
inc_df$shed <- shed
inc_df_long <- pivot_longer(inc_df, cols = 1:25, names_to = "age_targeted",
                            values_to = "incidence")
inc_df_long$age_targeted <- factor(inc_df_long$age_targeted, levels = as.character(c(1:49)))

# lower quantile (5%)
inc_l_df <- as.data.frame(inc_mat_l)
names(inc_l_df) <- c(age_targ_idx, 49)
inc_l_df <- cbind(inc_l_df, pars_vax)
inc_l_df$scenario <- scenario
inc_l_df$beta_rho <- beta_rho
inc_l_df$beta_rho_qual <- beta_rho_qual
inc_l_df$shed <- shed
inc_l_df_long <- pivot_longer(inc_l_df, cols = 1:25, names_to = "age_targeted",
                              values_to = "incidence_l")
inc_l_df_long$age_targeted <- factor(inc_l_df_long$age_targeted, levels = as.character(c(1:49)))

# upper quantile (95%)
inc_u_df <- as.data.frame(inc_mat_u)
names(inc_u_df) <- c(age_targ_idx, 49)
inc_u_df <- cbind(inc_u_df, pars_vax)
inc_u_df$scenario <- scenario
inc_u_df$beta_rho <- beta_rho
inc_u_df$beta_rho_qual <- beta_rho_qual
inc_u_df$shed <- shed
inc_u_df_long <- pivot_longer(inc_u_df, cols = 1:25, names_to = "age_targeted", names_prefix = "V",
                              values_to = "incidence_u")
inc_u_df_long$age_targeted <- factor(inc_u_df_long$age_targeted, levels = as.character(c(1:49)))

# merge central, lower and upper
inc_all <- merge(inc_df_long, inc_l_df_long) %>%
  merge(inc_u_df_long)

inc_df_long$scenario <- factor(inc_df_long$scenario, levels = c("base", 
                                                                "transblock",
                                                                "boostnat"))
inc_all$scenario <- factor(inc_all$scenario, levels = c("base", 
                                                        "transblock",
                                                        "boostnat"))
inc_all$beta_rho_qual <- factor(inc_all$beta_rho_qual, 
                                levels  = c("R0 = high, \U03C1 = 1",
                                            "R0 = high, \U03C1 = 3",
                                            "R0 = high, \U03C1 = 10",
                                            "R0 = medium, \U03C1 = 1",
                                            "R0 = medium, \U03C1 = 3",
                                            "R0 = medium, \U03C1 = 10",
                                            "R0 = low, \U03C1 = 1",
                                            "R0 = low, \U03C1 = 3",
                                            "R0 = low, \U03C1 = 10"))



##############################
# plot just the core results #
##############################

# counterfactual with no vaccination --> % change
cfact <- inc_all %>% filter(age_targeted == 49) %>%
  mutate(cf_inc = incidence,
         cf_inc_l = incidence_l,
         cf_inc_u = incidence_u)%>%
  select(-age_targeted, -incidence, -incidence_l, -incidence_u)

# merge with baseline counterfactual
df <- merge(inc_all, cfact, by = c("scenario", "beta_rho", "beta_rho_qual",
                                   "shed", "rho", "beta", "v_shed", "v_reduced_shed",
                                   "reduced_shed")) %>%
  mutate(inc_change = incidence/cf_inc,
         inc_change_l = incidence_u/cf_inc,
         inc_change_u = incidence_l/cf_inc)

# re-level factors for plotting
inc_core <- df%>%filter(age_targeted !=49, scenario == "transblock", shed == "1%") %>%
  mutate(beta_rho = factor(beta_rho,
                            levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                       "R0 = 14.0, \U03C1 = 3 yrs",
                                       "R0 = 14.0, \U03C1 = 10 yrs",
                                       "R0 = 7.0, \U03C1 = 1 yrs",
                                       "R0 = 7.0, \U03C1 = 3 yrs",
                                       "R0 = 7.0, \U03C1 = 10 yrs",
                                       "R0 = 3.5, \U03C1 = 1 yrs",
                                       "R0 = 3.5, \U03C1 = 3 yrs",
                                       "R0 = 3.5, \U03C1 = 10 yrs")),
         rho = round(1/(rho*360))) %>%
           mutate(rho = factor(rho),
                  beta = factor(beta)) 

# colour by transmission intensity, lty by rho

# plot by magnitude of reduction in incidence
# core results

beta.labs <- c(" ", " ", " ")
names(beta.labs) <- c(1, 0.5, 0.25)

q1 <- ggplot(data = inc_core)+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave(q1, filename = "figs/Fig4_panel1.png", height = 11, width = 11, unit = "in")
saveRDS(q1, file = "figs/q1.rds")

#############################
# plot all scenarios for SM #
#############################

scenario.labs <- c("S1: reduces susc. & transmission", 
                   "S2: reduces transmission only", 
                   "S3: ineffective in naive animals")
names(scenario.labs) <- c("base", "transblock", "boostnat")
shed.labs <- c(" ", " ")
names(shed.labs) <- c("1%", "50%")


q2 <- ggplot(data = df%>%filter(age_targeted !=49, scenario == "base", shed == "1%")%>%
               mutate(beta_rho = factor(beta_rho,
                                        levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                                    "R0 = 14.0, \U03C1 = 3 yrs",
                                                    "R0 = 14.0, \U03C1 = 10 yrs",
                                                    "R0 = 7.0, \U03C1 = 1 yrs",
                                                    "R0 = 7.0, \U03C1 = 3 yrs",
                                                    "R0 = 7.0, \U03C1 = 10 yrs",
                                                    "R0 = 3.5, \U03C1 = 1 yrs",
                                                    "R0 = 3.5, \U03C1 = 3 yrs",
                                                    "R0 = 3.5, \U03C1 = 10 yrs")),
                      rho = round(1/(rho*360))) %>%
               mutate(rho = factor(rho),
                      beta = factor(beta)) )+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1.2)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q3 <- ggplot(data = df%>%filter(age_targeted !=49, scenario == "boostnat", shed == "1%")%>%
               mutate(beta_rho = factor(beta_rho,
                                        levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                                    "R0 = 14.0, \U03C1 = 3 yrs",
                                                    "R0 = 14.0, \U03C1 = 10 yrs",
                                                    "R0 = 7.0, \U03C1 = 1 yrs",
                                                    "R0 = 7.0, \U03C1 = 3 yrs",
                                                    "R0 = 7.0, \U03C1 = 10 yrs",
                                                    "R0 = 3.5, \U03C1 = 1 yrs",
                                                    "R0 = 3.5, \U03C1 = 3 yrs",
                                                    "R0 = 3.5, \U03C1 = 10 yrs")),
                      rho = round(1/(rho*360))) %>%
               mutate(rho = factor(rho),
                      beta = factor(beta)) )+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1.2)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")


beta.labs <- c(" ", " ", " ")
names(beta.labs) <- c(0.125, 2.3/14, 3/14)

q4 <- ggplot(data = df%>%filter(age_targeted !=49, scenario == "transblock", shed == "50%")%>%
               mutate(beta_rho = factor(beta_rho,
                                        levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                                    "R0 = 14.0, \U03C1 = 3 yrs",
                                                    "R0 = 14.0, \U03C1 = 10 yrs",
                                                    "R0 = 7.0, \U03C1 = 1 yrs",
                                                    "R0 = 7.0, \U03C1 = 3 yrs",
                                                    "R0 = 7.0, \U03C1 = 10 yrs",
                                                    "R0 = 3.5, \U03C1 = 1 yrs",
                                                    "R0 = 3.5, \U03C1 = 3 yrs",
                                                    "R0 = 3.5, \U03C1 = 10 yrs")),
                      rho = round(1/(rho*360))) %>%
               mutate(rho = factor(rho),
                      beta = factor(beta)) )+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1.2)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q5 <- ggplot(data = df%>%filter(age_targeted !=49, scenario == "base", shed == "50%")%>%
               mutate(beta_rho = factor(beta_rho,
                                        levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                                    "R0 = 14.0, \U03C1 = 3 yrs",
                                                    "R0 = 14.0, \U03C1 = 10 yrs",
                                                    "R0 = 7.0, \U03C1 = 1 yrs",
                                                    "R0 = 7.0, \U03C1 = 3 yrs",
                                                    "R0 = 7.0, \U03C1 = 10 yrs",
                                                    "R0 = 3.5, \U03C1 = 1 yrs",
                                                    "R0 = 3.5, \U03C1 = 3 yrs",
                                                    "R0 = 3.5, \U03C1 = 10 yrs")),
                      rho = round(1/(rho*360))) %>%
               mutate(rho = factor(rho),
                      beta = factor(beta)) )+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1.2)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q6 <- ggplot(data = df%>%filter(age_targeted !=49, scenario == "boostnat", shed == "50%")%>%
               mutate(beta_rho = factor(beta_rho,
                                        levels  = c("R0 = 14.0, \U03C1 = 1 yrs",
                                                    "R0 = 14.0, \U03C1 = 3 yrs",
                                                    "R0 = 14.0, \U03C1 = 10 yrs",
                                                    "R0 = 7.0, \U03C1 = 1 yrs",
                                                    "R0 = 7.0, \U03C1 = 3 yrs",
                                                    "R0 = 7.0, \U03C1 = 10 yrs",
                                                    "R0 = 3.5, \U03C1 = 1 yrs",
                                                    "R0 = 3.5, \U03C1 = 3 yrs",
                                                    "R0 = 3.5, \U03C1 = 10 yrs")),
                      rho = round(1/(rho*360))) %>%
               mutate(rho = factor(rho),
                      beta = factor(beta)) )+
  geom_line(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                col = rho), lwd = 1.2)+
  geom_point(aes(x = as.numeric(age_targeted), y = 100 - 100*inc_change,
                 col = rho), size = 4, pch = "+")+
  geom_ribbon(aes(x = as.numeric(age_targeted), ymin = 100 - 100*inc_change_l, ymax = 100 - 100*inc_change_u,
                  fill = rho))+
  theme_minimal()+
  facet_wrap(~beta,labeller = labeller(beta = beta.labs), nrow = 3)+
  scale_x_continuous(breaks = c(0, 6, 12, 18, 24, 30, 36, 42, 48))+
  ylab("% reduction in incidence") +
  xlab("age targeted (months)") +
  ylim(0,100)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", direction = -1)+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect (yrs)", alpha = 0.15, direction = -1)+
  geom_vline(aes(xintercept = 6), col = "firebrick", lty = 2, lwd = 1)+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# arrange the plots in 2 rows
prow <- plot_grid(
  q1 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q2 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q3 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  q4 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q5 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q6 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  align = 'vh',
  hjust = -1,
  nrow = 2,
  rel_widths = c(1, -0.15, 1, -0.15, 1,
                 1, -0.15, 1, -0.15, 1)
)


legend <- get_legend(
  # create some space to the left of the legend
  q1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)


#create common x and y labels
y.grob <- textGrob("reduction in incidence (%)", 
                   gp=gpar(col="black", fontsize=20), rot=90)

x.grob <- textGrob("age targeted (months)", 
                   gp=gpar(col="black", fontsize=20))

#add to plot
qtot <- grid.arrange(arrangeGrob(prow, left = y.grob, bottom = x.grob))

# add the legend
qtot2 <- plot_grid(qtot, legend, nrow = 2, rel_heights = c(3, .3))

ggsave(filename = "figs/FigS7.png", height = 13, width = 16, unit = "in")


############################################
# plot change in weighted incidence by age #
# ############################################
output_list <- vector(mode = "list", length = dim(pars_vax)[1])
output_l_list <- vector(mode = "list", length = dim(pars_vax)[1])
output_comp <- vector(mode = "list", length = dim(pars_vax)[1])
names(output_comp) <- 1:(dim(pars_vax)[1])
output_comp_change <- output_comp
names(output_comp_change) <- 1:(dim(pars_vax)[1])

# extract the difference in incidence compared to without vaccination
for(i in 19:27){
  output_list[[i]] <- readRDS(paste("generated_data/vax_age_opti/w_inc/",
                                    "w_incidence_", i,
                                    ".rds", sep = ""))
  output_comp[[i]] <- rowSums(output_list[[i]][,19:49])
  output_comp_change[[i]] <- output_comp[[i]][1:24] - output_comp[[i]][rep(25,24)]
}
output_comp_change <- do.call(rbind.data.frame, output_comp_change)
names(output_comp_change) <- age_targ_idx[1:24]

# extract the average number of animals per age class in order to calculate incidence/1000
i = 1
vaxp <- rep(0, 49)
x <- sir_model_vax$new(alpha = alpha, beta = pars_vax$beta[i], gamma = gamma, sigma = waning, sigma_m = sigma_m, Ab_susc = Ab_susc,
                       mAb_susc = mAb_susc, reduced_shed = pars_vax$reduced_shed[i], mu = mu, N_0 = N_0,
                       importation_rate = importation_rate, imp_t = imp_t, delta = seasonality, ind1 = ind1, ind2 = ind2,
                       v_gamma = v_gamma, v_sigma = v_sigma, v_sigma_m = v_sigma_m, v_susc  = pars_vax$v_susc[i], v_mAb_susc = v_mAb_susc,
                       v_Ab_susc = pars_vax$v_Ab_susc[i], v_shed = pars_vax$v_shed[i], v_reduced_shed = pars_vax$v_reduced_shed[i],
                       vaxp = vaxp, rho = pars_vax$rho[i],
                       foi_bg_usr = foi_bg_usr)
n_runs <- 100
out <- as.data.frame(replicate(n_runs, x$run(t)[, 619:667]))

N_age_temp <- matrix(colMeans(out[7200:10800, ]), nrow = 100, ncol = 49, byrow = T)
N_age <- colMeans(N_age_temp)
N_age_mat <- matrix(replicate(24, N_age), nrow=24, byrow = T)
saveRDS(N_age, file = "generated_data/N_age_mat.rds")

N_age <- readRDS("generated_data/N_age_mat.rds")
N_18 <- sum(N_age[19:49]) # population over 18 months old

# merge with parameters from model runs
op <- cbind(output_comp_change, pars_vax[19:27,]) %>%
  mutate(beta_rho = factor(paste0("R0 = ", round(op$beta * 14, 1),
                                  ", 1/rho = ", round(1/(op$rho*360),0), " yrs"),
                           levels = c("R0 = 3.5, 1/rho = 1 yrs",
                                      "R0 = 3.5, 1/rho = 3 yrs",
                                      "R0 = 3.5, 1/rho = 10 yrs",
                                      "R0 = 7, 1/rho = 1 yrs",
                                      "R0 = 7, 1/rho = 3 yrs",
                                      "R0 = 7, 1/rho = 10 yrs",
                                      "R0 = 14, 1/rho = 1 yrs",
                                      "R0 = 14, 1/rho = 3 yrs",
                                      "R0 = 14, 1/rho = 10 yrs")))
# make long for plot
op_long <- pivot_longer(op, 1:24,names_to = "age_targeted",values_to = "change")
op_long$age_targeted <- factor(op_long$age_targeted,
                                             levels = age_targ_idx[1:24])
# plot core scenario
data <- op_long
max_change <-  max(1000*data$change/(10*N_18)) # divide by number of years * number of camels
min_change <-  min(1000*data$change/(10*N_18))
ggplot(data)+
  geom_tile(aes(x = beta_rho, y = age_targeted, fill = 1000*change/(10*N_18)))+
  theme_minimal()+
  scale_fill_gradientn(name = "Change in mean\nannual weighted\nincidence per 1000\n>18-month-olds",
                       colours = c("blue","white","red"),
                       values = scales::rescale(c(min_change, 0, max_change)),
                       guide = "colorbar", limits=c(min_change, 11),
                       breaks = c(+10, 0, -15, -30, - 45, -60))+
  ylab(" ")+
  xlab(" ")+
  scale_x_discrete(labels = rep(c("             ", "              ", "               "),
                                times = 4))+
  theme(text = element_text(size = 19),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5), col = "white", lwd = 1)

ggsave(filename = "figs/FigS6.png", height = 8, width = 10, unit = "in")
