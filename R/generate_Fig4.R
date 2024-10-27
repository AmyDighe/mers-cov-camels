##################################
# read in vaccine impact results #
##################################

# infectiousness and susceptibility parameters under 3 vaccine efficacy scenarios, 
# 2 differing assumed relationships between RNA shedding and infectiousness 
# and 3 different transmission intensities 
pars_vax <- readRDS("R/pars_vax.rds")
# year index lower and upper
yr_idx_l <- (seq(20*360, 34*360, by = 360) - 7200)/10
yr_idx_u <- yr_idx_l + 35
# vaccination coverage vector
coverage <- c(0, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# patch structure
n_r <- 5 # number of rows in grid of sub-pops
n_c <- 5 # number of cols in grid of sub-pops
N_patch <- n_r*n_c

# define arrays
patch_years <- array(NA, dim = c(54, N_patch, length(yr_idx_l), length(coverage)))
an_Re <- array(NA, dim = c(54, length(yr_idx_l), length(coverage), 100))
an_inc <- an_Re
an_winc <- an_Re
persist <- matrix(NA, nrow = 54, ncol = length(coverage))
inc <- persist
winc <- persist
Re <- persist

# populate with saved outputs for population size 1
for(j in 1:54){# (dim(pars_vax)[1])){
  persist[j,] <- readRDS(paste0("generated_data/persistence/av_persist_2mil", j, "_1.rds"))
  inc[j,] <- readRDS(paste0("generated_data/persistence/av_inc_2mil", j, "_1.rds"))
  #winc[j,] <- readRDS(paste0("generated_data/persistence/av_winc_2mil", j, "_1.rds"))
  an_inc[j,,,] <- readRDS(paste0("generated_data/persistence/an_inc_2mil", j, "_1.rds"))
  #an_winc[j,,,] <- readRDS(paste0("generated_data/persistence/an_winc", j, "_1.rds"))
  #an_Re[j,,,] <- readRDS(paste0("generated_data/persistence/an_Re_2mil", j, "_1.rds"))
  patch_years[j,,,] <- readRDS(paste0("generated_data/persistence/patch_years_2mil", j, "_1.rds"))
}

patch_years2 <- array(NA, dim = c(dim(pars_vax)[1], N_patch, length(yr_idx_l), length(coverage)))
an_Re2 <- array(NA, dim = c(dim(pars_vax)[1], length(yr_idx_l), length(coverage), 100))
an_inc2 <- an_Re2
an_winc2 <- an_Re2
persist2 <- matrix(NA, nrow = dim(pars_vax)[1], ncol = length(coverage))
inc2 <- persist2
winc2 <- persist2
Re2 <- persist2

# populate with saved outputs for population size 2
for(j in 1:(dim(pars_vax)[1])){
  persist2[j,] <- readRDS(paste0("generated_data/persistence/av_persist", j, "_2.rds"))
  inc2[j,] <- readRDS(paste0("generated_data/persistence/av_inc", j, "_2.rds"))
  winc2[j,] <- readRDS(paste0("generated_data/persistence/av_winc", j, "_2.rds"))
  an_inc2[j,,,] <- readRDS(paste0("generated_data/persistence/an_inc", j, "_2.rds"))
  #an_winc2[j,,,] <- readRDS(paste0("generated_data/persistence/an_winc", j, "_2.rds"))
  an_Re2[j,,,] <- readRDS(paste0("generated_data/persistence/an_Re", j, "_2.rds"))
  patch_years2[j,,,] <- readRDS(paste0("generated_data/persistence/patch_years", j, "_2.rds"))
}

###########################################################
## produce Figure XB: impact of vaccination of incidence ##
###########################################################

# reshape annual incidence data --> total with upper and lower bound
## take just first ten years following vaccination
df <- an_inc[,1:10,,]
## sum over these 10 years
df <- rowSums(aperm(df, c(1,3,4,2)), dims = 3)
## make long with coverage and iteration as variables
library(R.utils)
df_test <- as.data.frame(wrap(df, map=list(1, NA)))
df_test <- cbind(df_test, pars_vax[1:54,])
df_test$set <- 1:54
df_tmp <- df_test %>% pivot_longer(cols = 1:800, names_to = "cov_iter", values_to = "incidence")
df_tmp$coverage <- rep(coverage, 100*54)
df_tmp$iteration <- rep(rep(1:100, each = 8), 54)

## calculate % reduction using coverage = 0 as baseline
baseline <- df_tmp %>% filter(coverage == 0)%>% rename(baseline = incidence)
df_tmp2 <- merge(baseline %>% select(set, iteration, baseline), 
                 df_tmp, by = c("set", "iteration"))%>%
  mutate(reduction = 100*(baseline - incidence)/baseline)

## df with means and quantiles across 100 runs
df_plot <- df_tmp2 %>% 
  group_by(set, coverage, beta, rho) %>%
  summarise(mean = mean(reduction),
            lower = quantile(reduction, probs = 0.025),
            upper = quantile(reduction, probs = 0.975))

## plot

beta.labs <- c("R0 = 14", "R0 = 7", "R0 = 3.5")
names(beta.labs) <- c(1, 0.5, 0.25)

p0 <- ggplot(df_plot %>% filter(set >=19, set <= 27, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                      labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_rep_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3, strip.position = "right")+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")


p1 <- readRDS("figs/q1.rds")
ggsave(p1, file = "figs/Fig4_update.png", height = 10, width = 9, unit = "in")

legend <- get_legend(p0)

plots <- plot_grid(p1+theme(legend.position = "none"), NULL,
          p0+theme(legend.position = "none"), rel_widths = c(1, -0.15, 1), nrow = 1)

plot_tot <- plot_grid(plots, legend, nrow = 2, rel_heights = c(1, 0.1))

ggsave(plot_tot, file = "figs/Fig4.png", height = 10, width = 18, unit = "in")

# add in 75,000 pop equivalent

# populate with saved outputs for population size 2
for(j in 1:(dim(pars_vax)[1])){
  persist2[j,] <- readRDS(paste0("generated_data/persistence/av_persist", j, "_2.rds"))
  inc2[j,] <- readRDS(paste0("generated_data/persistence/av_inc", j, "_2.rds"))
  winc2[j,] <- readRDS(paste0("generated_data/persistence/av_winc", j, "_2.rds"))
  an_inc2[j,,,] <- readRDS(paste0("generated_data/persistence/an_inc", j, "_2.rds"))
  #an_winc2[j,,,] <- readRDS(paste0("generated_data/persistence/an_winc", j, "_2.rds"))
  an_Re2[j,,,] <- readRDS(paste0("generated_data/persistence/an_Re", j, "_2.rds"))
  patch_years2[j,,,] <- readRDS(paste0("generated_data/persistence/patch_years", j, "_2.rds"))
}

# reshape annual incidence data --> total with upper and lower bound
## take just first ten years following vaccination
df <- an_inc2[,1:10,,]
## sum over these 10 years
df <- rowSums(aperm(df, c(1,3,4,2)), dims = 3)
## make long with coverage and iteration as variables
df_test <- as.data.frame(wrap(df, map=list(1, NA)))
df_test <- cbind(df_test, pars_vax[1:54,])
df_test$set <- 1:54
df_tmp <- df_test %>% pivot_longer(cols = 1:800, names_to = "cov_iter", values_to = "incidence")
df_tmp$coverage <- rep(coverage, 100*54)
df_tmp$iteration <- rep(rep(1:100, each = 8), 54)

## calculate % reduction using coverage = 0 as baseline
baseline <- df_tmp %>% filter(coverage == 0)%>% rename(baseline = incidence)
df_tmp2 <- merge(baseline %>% select(set, iteration, baseline), 
                 df_tmp, by = c("set", "iteration"))%>%
  mutate(reduction = ifelse(test = baseline > 10000, yes = (100*(baseline - incidence)/baseline), no = NA))

## df with means and quantiles across 100 runs
df_plot2 <- df_tmp2 %>% filter(set >=19, set <= 27, baseline >0)%>% 
  group_by(set, coverage, beta, rho) %>%
  drop_na(reduction)%>%
  summarise(mean = mean(reduction),
            lower = quantile(reduction, probs = 0.025),
            upper = quantile(reduction, probs = 0.975))%>%
  mutate(lower = if_else(lower<0, 0, lower))

p2 <- ggplot(df_plot2 %>% filter(set >=19, set <= 27, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  facet_rep_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3, strip.position = "right")+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

plots <- plot_grid(p0+theme(legend.position = "none"), NULL,
                   p2+theme(legend.position = "none"), rel_widths = c(1, -0.15, 1), nrow = 1)

plot_tot <- plot_grid(plots, legend, nrow = 2, rel_heights = c(1, 0.1))

ggsave(plot_tot, file = "figs/Fig5.png", height = 10, width = 18, unit = "in")
ggsave(p2, file = "figs/Fig4B_smallpop.png", height = 10, width = 9, unit = "in")
saveRDS(df_plot2, file = "fig4B_small_pop.rds")

##################################################
## produce Table X: persistence and vaccination ##
##################################################

# reshape persistence data
pf <- cbind(as.data.frame(persist),pars_vax[1:50,])
pf$set <- 1:50
names(pf)[1:8] <- coverage
pf <- pf %>% pivot_longer(cols = 1:8, names_to = "coverage", 
                          values_to = "persistence")%>%
  filter(set >=19, set <= 27)

# save table
write.csv(pf, file = "tables/persistence_vax.csv")

# same for second population
pf2 <- cbind(as.data.frame(persist2),pars_vax)
pf2$set <- 1:54
names(pf2)[1:8] <- coverage
pf2 <- pf2 %>% pivot_longer(cols = 1:8, names_to = "coverage", 
                          values_to = "persistence")%>%
  filter(set >=19, set <= 27)

# save table
write.csv(pf2, file = "tables/persistence_vax_pop2.csv")

###################################################################
## difference in reduction in incidence between scenario 1 and 2 ##
###################################################################

# first get the df_plot equiv for scenario 2 for each pop size
# pop 1
df_plot3 <- df_plot %>% filter(set >=1, set <= 9)
# pop 2
df_plot4 <- df_tmp2 %>% filter(set >=1, set <= 9, baseline >0)%>% 
  group_by(set, coverage, beta, rho) %>%
  drop_na(reduction)%>%
  summarise(mean = mean(reduction),
            lower = quantile(reduction, probs = 0.025),
            upper = quantile(reduction, probs = 0.975))%>%
  mutate(lower = if_else(lower<0, 0, lower))

# merge and calculate the difference
dp1 <- merge(df_plot %>% filter(set >=19, set <= 27, coverage >0) %>%
               mutate(set = set - 18), df_plot3 %>% filter(coverage > 0), 
             by = c("set", "coverage", "beta", "rho"))%>%
  mutate(dif = mean.y - mean.x)

dp2 <- merge(df_plot2 %>% filter(set >=19, set <= 27, coverage >0) %>%
               mutate(set = set - 18), df_plot4 %>% filter(coverage > 0), 
             by = c("set", "coverage", "beta", "rho"))%>%
  mutate(dif = mean.y - mean.x)

#######################################################################
## reduction in incidence by coverage for all 6 scenarios for the SM ##
#######################################################################

## for the population of 2 million

beta.labs <- c(" ", " ", " ")
names(beta.labs) <- c(1, 0.5, 0.25)

q0 <- ggplot(df_plot %>% filter(set >=19, set <= 27, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q1 <- ggplot(df_plot %>% filter(set >=1, set <= 9, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q2 <- ggplot(df_plot %>% filter(set >=28, set <= 36, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
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

q3 <- ggplot(df_plot %>% filter(set >=37, set <= 45, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q4 <- ggplot(df_plot %>% filter(set >=10, set <= 18, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
  theme(text = element_text(size = 22),
        axis.title.x = element_text(vjust= -1),
        axis.title.y = element_text(vjust = +3),
        plot.margin = unit(c(0,1,1,1), "cm"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

q5 <- ggplot(df_plot %>% filter(set >=46, set <= 54, coverage > 0))+
  geom_point(aes(x = 100*coverage, y = mean, col = as.factor(rho)), size = 4, pch = "+")+
  geom_line(aes(x = 100*coverage, y = mean, col = as.factor(rho), group = as.factor(rho)), lwd = 1)+
  geom_ribbon(aes(x = 100*coverage, ymin = lower, ymax = upper,
                  group = as.factor(rho), fill = as.factor(rho)), alpha = 0.3)+
  scale_color_viridis_d(option = "C", name = "duration of vaccine effect",
                        labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_fill_viridis_d(option = "C", name = "duration of vaccine effect",
                       labels = c("10 yrs", "3 yrs", "1 yr"))+
  scale_x_continuous(breaks = c(40, 50, 60, 70,  80, 90, 100))+
  facet_wrap(~beta, labeller = labeller(beta = beta.labs), nrow = 3)+
  theme_minimal()+
  xlab("% vaccine coverage")+
  ylim(0,100)+
  ylab(" ")+
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
  q0 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q1 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q2 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  q3 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q4 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
  NULL,
  q5 + theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
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

x.grob <- textGrob("coverage (%)", 
                   gp=gpar(col="black", fontsize=20))

#add to plot
qtot <- grid.arrange(arrangeGrob(prow, left = y.grob, bottom = x.grob))

# add the legend
qtot2 <- plot_grid(qtot, legend, nrow = 2, rel_heights = c(3, .3))



ggsave(qtot2, file = "figs/4B_SM.png", height = 13, width = 16, unit = "in")




###################################################
## produce Table SX: persistence and vaccination ##
###################################################

# reshape persistence data
pf_all <- cbind(as.data.frame(persist),pars_vax[1:54,])
pf_all$set <- 1:54
names(pf_all)[1:8] <- coverage
pf_all <- pf_all %>% pivot_longer(cols = 1:8, names_to = "coverage", 
                          values_to = "persistence")

# save table
write.csv(pf_all, file = "tables/persistence_vax_SM1.csv")

# same for second population
pf2_all <- cbind(as.data.frame(persist2),pars_vax)
pf2_all$set <- 1:54
names(pf2_all)[1:8] <- coverage
pf2_all <- pf2_all %>% pivot_longer(cols = 1:8, names_to = "coverage", 
                            values_to = "persistence")

# save table
write.csv(pf2_all, file = "tables/persistence_vax_pop2_SM.csv")
