############################################
# look at the overall periodicity patterns #
############################################
sp <- readRDS(file = "generated_data/period.rds")
par_grid_period <- sp$par_grid

period_sp <- sp$period

cuts <- apply(period_sp, 2, cut, c(0, 340, 380, 700, 740, 1060, 1100, 1420, 1460, Inf), 
              labels= c("other", "1 yr", "other", "2 yrs", "other", "3 yrs", "other", "4 yrs", "other"))
cuts_df <- cbind(period = t(cuts), par_grid_period) %>%
  pivot_longer(cols = 1:100, names_to = "run", values_to = "period")%>%
  mutate(period = factor(period, levels = c(NA, "other", "4 yrs",
                                            "3 yrs", "2 yrs", "1 yr")),
         pop = factor(pop),
         R0 = factor(case_when((beta == (3.5/14) | beta == (2/14) | beta == (1.75/14)) ~ "low",
                               (beta == (7/14) | beta == (5/14) | beta == (2.3/14)) ~ "medium",
                               (beta == (14/14) | beta == (5/14) | beta == (3/14)) ~ "high",
                               TRUE~as.character(NA)), levels = c("low", "medium", "high")))

# with seasonality of births half-strength
g1 <- ggplot(data = cuts_df%>%filter(seasonality == 0.5, shedding != 0.25, !is.na(period)), 
             aes(x = pop, fill = period))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c("lightgrey", '#fde725', '#35b779', '#31688e', '#472d7b'))+
  theme_minimal()+
  theme(text= element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(shedding~R0, scales = "free")+
  ylab("proportion of endemic simulations")+
  xlab("population size")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")+
  guides(fill = guide_legend(reverse=TRUE))

ggsave(g1, filename = "figs/barchart_period.png")

# with seasonality of births full-strength
g2 <- ggplot(data = cuts_df%>%filter(seasonality == 1, shedding != 0.25, !is.na(period)), 
             aes(x = pop, fill = period))+
  geom_bar(position = "fill")+
  scale_fill_manual(values = c("lightgrey", '#fde725', '#35b779', '#31688e', '#472d7b'))+
  theme_minimal()+
  theme(text= element_text(size = 25))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(shedding~R0, scales = "free")+
  ylab("proportion of endemic simulations")+
  xlab("population size")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom")+
  guides(fill = guide_legend(reverse=TRUE))
ggsave(g2, filename = "figs/barchart_period2.png")

#######################################################################################
## example with annual periodicities - large population for clarity of overall trend ##
#######################################################################################

# extract proportional break down of periods
cut_sum <- cuts_df %>% 
  group_by(pop, shedding, beta, waning, susc, seasonality, persist, persist2, R0) %>%
  summarise(p1 = 100* mean(period == "1 yr", na.rm = T),
            p2 = 100* mean(period == "2 yrs", na.rm = T),
            p3 = 100* mean(period == "3 yrs", na.rm = T),
            p4 = 100* mean(period == "4 yrs", na.rm = T),
            other = 100*mean(period == "other", na.rm = T))%>%
  mutate(pop = as.numeric(as.character(pop)))

# indexes of good examples of annual transmission in par_grid_period
idx_select <- c(337:339, 343:348)

# run these parameter sets
# for(i in idx_select){
#   x <- sir_model$new(alpha = alpha, beta = cut_sum$beta[i], gamma = gamma, 
#                      sigma = cut_sum$waning[i], sigma_m = omega, 
#                      Ab_susc = cut_sum$susc[i], mAb_susc = mAb_susc, 
#                      reduced_shed = cut_sum$shedding[i], mu = mu,
#                      N_0 = 1000000000, importation_rate = importation_rate, 
#                      imp_t = 1, delta = cut_sum$seasonality[i], 
#                      ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)
#   
#   out <- as.data.frame(replicate(100, x$run(t)[, c(354, 369)])) # extract Itot and birthrate
#   saveRDS(out, file = paste("generated_data/time_series_examples_period/out", i, ".rds", sep = ""))
# }


# Combine all outputs
annual_examples <- list(length = length(idx_select))
for(i in 1:(length(idx_select))){
  annual_examples[[i]] <- format_out(filepath = paste0("generated_data/time_series_examples_period/out", idx_select[i], ".rds"),
                                     time_period = 12600, par_grid  = cut_sum, n = 2) 
}

annual_examples_df <- bind_rows(annual_examples)

annual_examples_df <- annual_examples_df %>%
  mutate(seasonality = factor(seasonality, levels = c(0, 0.5, 1)))

# subset data when shedding = 0.01 and R0 is high
s1_r0h <- annual_examples_df %>% filter(shedding == 0.01, R0 == "high",
                                        (t > 360*25)& (t<360*28))

p1 <- ggplot()+
  geom_line(data = s1_r0h %>% filter(var == "Itot"),
            aes(x = t, y = values, group = interaction(run_no, seasonality), colour = seasonality), 
            alpha = 0.5, linewidth = 1)+
  scale_color_manual("strength of\nseasonality\nof calving", 
                     values = c("mediumorchid1", "mediumpurple1","mediumpurple4"),
                     labels = c("0%", "50%", "100%"))+
  scale_x_continuous(name = " ", breaks = c(seq(360*25, 360*28, by = 360)),
                     labels = c(1:4))+
  annotate("rect", xmin = c((360*25), (360*25) + 270, (360*25) + 630, (360*25) + 990),
           xmax = c((360*25) + 90, (360*25) + 450, (360*25) + 810, 360*28), ymin = min(s1_r0h %>% filter(var == "Itot") %>% pull(values)) -1000, 
           ymax = max(s1_r0h %>% filter(var == "Itot") %>% pull(values)) + 1000,
           alpha = .05,fill = "blue")+
  theme_classic()+
  theme(text = element_text(size = 18),
        legend.key.width = unit(1, 'cm'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14))+
  ylab(" ")+
  ggtitle("High transmission intensity setting, low relative infectiousness of reinfections (1%)")

# subset data when shedding = 0.01 and R0 is medium
s1_r0m <- annual_examples_df %>% filter(shedding == 0.01, R0 == "medium",
                                        (t > 360*25)& (t<360*28))

p2 <- ggplot(data = s1_r0m %>% filter(var == "Itot"))+
  geom_line(aes(x = t, y = values, group = interaction(run_no, seasonality), colour = seasonality), 
            alpha = 0.5, linewidth = 1)+
  scale_color_manual("strength of\nseasonality\nof calving",
                     values = c("mediumorchid1", "mediumpurple1","mediumpurple4"),
                     labels = c("0%", "50%", "100%"))+
  theme_classic()+
  scale_x_continuous(name = " ", breaks = c(seq(360*25, 360*28, by = 360)),
                     labels = c(1:4))+
  annotate("rect", xmin = c((360*25), (360*25) + 270, (360*25) + 630, (360*25) + 990),
           xmax = c((360*25) + 90, (360*25) + 450, (360*25) + 810, 360*28), ymin = min(s1_r0m %>% filter(var == "Itot") %>% pull(values)) -1000, 
           ymax = max(s1_r0m %>% filter(var == "Itot") %>% pull(values)) + 1000,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 18),
        legend.key.width = unit(1, 'cm'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14))+
  ylab("Infectious individuals")+
  ggtitle("Moderate transmission intensity, low relative infectiousness of reinfections (1%)")

# subset data when shedding = 0.50 and R0 is high
s50_r0h <- annual_examples_df %>% filter(shedding == 0.50, R0 == "high",
                                         (t > 360*25)& (t<360*28))
p3 <- ggplot()+
  geom_line(data = s50_r0h %>% filter(var == "Itot"),
            aes(x = t, y = values, group = interaction(run_no, seasonality), colour = seasonality), 
            alpha = 0.5, linewidth = 1)+
  scale_color_manual("strength of\nseasonality\nof calving",
                     values = c("mediumorchid1", "mediumpurple1","mediumpurple4"),
                     labels = c("0%", "50%", "100%"))+
  theme_classic()+
  scale_x_continuous(name = "time (years)", breaks = c(seq(360*25, 360*28, by = 360)),
                     labels = c(1:4))+
  annotate("rect", xmin = c((360*25), (360*25) + 270, (360*25) + 630, (360*25) + 990),
           xmax = c((360*25) + 90, (360*25) + 450, (360*25) + 810, 360*28), ymin = min(s50_r0h %>% filter(var == "Itot") %>% pull(values)) -1000, 
           ymax = max(s50_r0h %>% filter(var == "Itot") %>% pull(values)) + 1000,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 18),
        legend.key.width = unit(1, 'cm'),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 14))+
  ylab(" ")+
  ggtitle("High transmission intensity setting, high relative infectiousness of reinfections (50%)")


p4 <- ggplot()+
  geom_line(data = s50_r0h %>% filter(var == "Itot", t > 360*27),
            aes(x = t, y = values, group = interaction(run_no, seasonality, var), colour = seasonality, lty = var), 
            alpha = 0.5, linewidth = 1)+
  geom_line(data = s50_r0h %>% filter(var == "birthrate", t > 360*27, seasonality != "0"),
            aes(x = t, y = (values*18)+53000000, group = interaction(run_no, seasonality, var), colour = seasonality, lty = var), 
            alpha = 0.5, linewidth = 1)+
  scale_linetype_manual("", values = c(2, 1), labels = c("Births", "Infections"))+
  scale_color_manual("strength of\nseasonality\nof calving",
                     values = c("mediumorchid1", "mediumpurple1","mediumpurple4"),
                     labels = c("0%", "50%", "100%"))+
  theme_classic()+
  scale_y_continuous(
    name = "Infectious individuals",
    sec.axis = sec_axis(~ (.-53000000)/18, name= "Births")
  )+
  scale_x_continuous(name = "time of year", breaks = c(seq((360*28) - 345, (360*28)-15, by = 30), 
                                                       seq(((360*28) - 360), 360*28, by = 30)),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "June",
                                "July", "Aug", "Sept", "Oct","Nov", "Dec",
                                rep(c(""), 13)))+
  annotate("rect", xmin = c((360*27), (360*27) + 270),
           xmax = c((360*27) + 90, (360*27) + 360), ymin = min(s50_r0h %>% filter(var == "birthrate") %>% pull(values)* 18) +53000000 -1000, 
           ymax = max(s50_r0h %>% filter(var == "birthrate") %>% pull(values)* 18) + 53000000 + 1000,
           alpha = .05,fill = "blue")+
  theme(text = element_text(size = 18),
        axis.title.x = element_text(vjust=-0.5),
        axis.title.y = element_text(vjust= 2),
        axis.ticks.x = element_line(color = c(rep(NA, 13), rep("black", 13))),
        axis.line.x.bottom=element_line(color="black"),
        panel.grid.major.x = element_line(color = c(rep(NA, 13), rep("gray95", 13))),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        legend.box.background = element_rect(color = "black"),
        plot.title = element_text(size = 16))+
  guides(colour = guide_legend(order = 2),
         linetype = guide_legend(order = 1))+
  ggtitle(" ")

pcol1 <- plot_grid(p1+theme(legend.position = "none"),
                   p2+theme(legend.position = "none"),
                   p3+theme(legend.position = "none"),
                   p4, ncol = 1, rel_heights = c(1,1,1,2), align = "v")

ggsave(pcol1, height = 12, width = 9, unit = "in", file = "figs/Fig3_panel1.png")


############################################
## examples of extra-annual periodicities ##
############################################

# indexes of good examples of annual transmission in par_grid_period
idx_select_extra_an <- c(342, 237, 306)

# run these parameter sets
persist <- vector(length = (dim(cut_sum)[1]))
period <- matrix(ncol = dim(cut_sum)[1], nrow = 100)

for(i in idx_select_extra_an){
  i = 306
  x <- sir_model$new(alpha = alpha, beta = cut_sum$beta[i], gamma = gamma,
                     sigma =cut_sum$waning[i], sigma_m = omega,
                     Ab_susc = cut_sum$susc[i], mAb_susc = mAb_susc,
                     reduced_shed = cut_sum$shedding[i], mu = mu,
                     N_0 = cut_sum$pop[i], importation_rate = importation_rate,
                     imp_t = 1, delta = cut_sum$seasonality[i],
                     ind1 = ind1, ind2 = ind2, foi_bg_usr = foi_bg_usr)
  
  out <- as.data.frame(replicate(100, x$run(t)[, c(354)])) # extract Itot
  
  persist[i] <- sum(out[360*35, ] > 0) # 35yrs after model initiation (25yrs after background FOI switched off)
  
  if(persist[i] > 50){
    
    idx_persist <- which(out[360*35, ] > 0)
    out_persist <- out[,idx_persist]
    
    for(j in 1:(persist[i])){
      ACF <- acf(out_persist[(360*25):(360*35), j], lag.max = 360*5, plot = FALSE)
      acf_df <- data.frame(acf = ACF$acf[180:(360*5)], lag = ACF$lag[180:(360*5)])
      # write something in here about making sure acf is above the significance level
      period[j,i] <- acf_df[which.max(acf_df$acf), ]$lag 
    }
  } else {
    period[,i] <- NA
  }
  
  print(paste0(i))
  
  saveRDS(out, file = paste("generated_data/time_series_examples_period/out_extra_an", i, ".rds", sep = ""))
}

saveRDS(period[,idx_select_extra_an], file = "generated_data/time_series_examples_period/period_extra_an.rds")

outB <- format_out_extra(filepath = "generated_data/time_series_examples_period/out_extra_an342.rds",
                         time_period = 12600, par_grid  = cut_sum, n = 1)

period_recovered <- as.data.frame(readRDS("generated_data/time_series_examples_period/period_extra_an.rds"))
period_recovered$run_no <- paste0("V", 1:100)

# merge period with outB to identify biennial examples
outB <- merge(outB, period_recovered[,c(1,4)], by = "run_no") %>%
  rename(period = V1)

# plot the biennial example

pb <- ggplot(data = outB %>% filter(t>(360*25), (period > 710) & (period < 730)))+
  geom_line(aes(x = (t - (360*25))/360, y =values, group = run_no), col = '#31688e', alpha = 0.1)+
  geom_line(data = outB %>% filter(t>(360*25), run_no == "V64"), 
            aes(x = (t - (360*25))/360, y =values), col = '#31688e', linewidth = 1.5)+
  xlab(" ")+
  ylab(" ")+
  theme_classic()+  
  ggtitle("2-year periodicity")+
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5))


outT <- format_out_extra(filepath = "generated_data/time_series_examples_period/out_extra_an237.rds",
                         time_period = 12600, par_grid  = cut_sum, n = 1)
# merge period with outB to identify biennial examples
outT <- merge(outT, period_recovered[,c(2,4)], by = "run_no") %>%
  rename(period = V2)

# plot the triennial example

pt <- ggplot(data = outT %>% filter(t>(360*25), (period > 1060) & (period < 1100)))+
  geom_line(aes(x = (t - (360*25))/360, y =values, group = run_no), col = '#35b779', alpha = 0.1)+
  geom_line(data = outT %>% filter(t>(360*25), run_no == "V72"), 
            aes(x = (t - (360*25))/360, y =values), col = '#35b779', linewidth = 1.5)+
  xlab(" ")+
  ylab("infectious individuals")+
  theme_classic()+  
  ggtitle("3-year periodicity")+
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5))


outQ <- format_out_extra(filepath = "generated_data/time_series_examples_period/out_extra_an306.rds",
                         time_period = 12600, par_grid  = cut_sum, n = 1)
# merge period with outB to identify triennial examples
outQ <- merge(outQ, period_recovered[,c(3,4)], by = "run_no") %>%
  rename(period = V3)



ggsave(pcol2, file = "figs/Fig3_panel2.png", height = 12, width = 9, unit = "in")

ggsave(ptot, file = "figs/Fig3_whole.png", height = 12, width = 18, unit = "in")