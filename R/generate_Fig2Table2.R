# plot the table and graph showing CCS under different model assumptions

# read in persistence estimates from single patch model
results_sp <- readRDS(here("generated_data", "ccs", "persistence.rds")) %>%
  mutate(seasonality = recode_factor(seasonality, `1` = 1.0, `2` = 0.5))

# read in persistence estimates from metapopulation model
results_mp <- readRDS(here("generated_data", "ccs", "persistence_mp.rds"))%>%
  mutate(seasonality = recode_factor(seasonality, `1` = 1.0, `2` = 0.5))


## Table of CCS estimates ######################################################

# interpolate to get the CCS (pop at which size persistence = 50% chance)

## single patch
minp <- results_sp %>% group_by(beta, shedding, seasonality) %>% 
  filter(persist<=50) %>%
  summarise(population = max(pop),
            persist = persist[which.max(pop)])

maxp <- results_sp %>% group_by(beta, shedding, seasonality) %>% 
  filter(persist>50) %>%
  summarise(population = min(pop),
            persist = persist[which.min(pop)])

ccs_sp <- rbind(minp, maxp) %>% 
  group_by(beta, shedding, seasonality) %>%
  summarise(ccs = round((approx(x = persist, 
                               y = population, 
                               xout = 51, 
                               method = "linear")$y), -2))

## meta-pop
minp_mp <- results_mp %>% group_by(beta, shedding, seasonality) %>% 
  filter(persist<=50) %>%
  summarise(population = max(pop),
            persist = persist[which.max(pop)])

maxp_mp <- results_mp %>% group_by(beta, shedding, seasonality) %>% 
  filter(persist>50) %>%
  summarise(population = min(pop),
            persist = persist[which.min(pop)])

ccs_mp <- rbind(minp_mp, maxp_mp) %>% 
  group_by(beta, shedding, seasonality) %>%
  summarise(ccs = round((approx(x = persist, 
                          y = 25*population, 
                          xout = 51, 
                          method = "linear")$y), -2))

# merge
ccs <- merge(ccs_sp, ccs_mp, by = c("shedding", "beta", "seasonality"))

# arrange table
names(ccs)[4:5] <- c("ccs_sp", "ccs_mp")
ccs_d1 <- ccs %>% filter(seasonality == 1.0)%>%
  rename(ccs_sp_1 = ccs_sp,
         ccs_mp_1 = ccs_mp)
ccs_d05 <- ccs %>% filter(seasonality == 0.5) %>% 
  select(-beta, -shedding, -seasonality) %>%
  rename(ccs_sp_05 = ccs_sp,
         ccs_mp_05 = ccs_mp)
ccs_tab <- cbind(ccs_d1, ccs_d05) %>%
  mutate(R0 = 14*beta)

write.csv(ccs_tab, file = here("tables", "table3.csv"))



## Plot persistence by population size #########################################
  
# labels for multipanel plot
## recode seasonality as a factor
shedding.labs <- c("rinf = 1%", "rinf = 25%", "rinf = 50%")
names(shedding.labs) <- unique(results_sp$shedding)

seasonality.labs <- c("\U03B4 = 1.0", "\U03B4 = 0.5")
names(seasonality.labs) <- rev(unique(results_sp$seasonality))

null.labs <- c(" ", " ")
names(null.labs) <- rev(unique(results_sp$seasonality))

# create plots
shed_1_p <- ggplot(results_sp %>% filter(shedding == 0.01))+
  geom_line(aes(x = pop, y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = seasonality.labs,
                                                           shedding = shedding.labs),
                 scale = "free")+
  scale_x_continuous(limits=c(1e+02, 1e+06), trans = "log10") +
  scale_y_continuous(limits=c(0,100))+
  geom_hline(yintercept = 50, lwd = 0.5, lty = 2)+
  theme_minimal()+
  scale_color_manual(name = "R0", labels = c("low", "moderate", "high"),
                     values = c("hotpink2", "navy", "olivedrab1"))+
  theme(text = element_text(size = 20))+
  xlab("population size (log10 scale)")+
  ylab(" ")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  theme(axis.line =element_line(color="black"))+
  geom_line(data = results_mp %>% filter(shedding == 0.01),
            aes(x = 25*pop, y = persist,  col = as.factor(beta*14)), size = 1, lty = 2)+
  theme(legend.position = "bottom")

shed_25_p <- ggplot(results_sp %>% filter(shedding == 0.25))+
  geom_line(aes(x = pop, y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = null.labs,
                                                           shedding = shedding.labs),
                 scale = "free")+
  scale_x_continuous(limits=c(1e+02, 1e+06), trans = "log10") +
  scale_y_continuous(limits=c(0,100))+
  geom_hline(yintercept = 50, lwd = 0.5, lty = 2)+
  theme_minimal()+
  scale_color_manual(name = "R0", labels = c("3.50", "7.00", "14.00"),
                     values = c("salmon1", "cornflowerblue", "goldenrod"))+
  theme(text = element_text(size = 20))+
  xlab("population size (log10 scale)")+
  ylab("persistence (% of stochastic model runs)")+
  theme(axis.title.x = element_blank())+
  theme(axis.line =element_line(color="black"))+
  geom_line(data = results_mp %>% filter(shedding == 0.25),
            aes(x = 25*pop, y = persist,  col = as.factor(beta*14)), size = 1, lty = 2)+
  theme(legend.position = "right")

shed_50_p <- ggplot(results_sp %>% filter(shedding == 0.50))+
  geom_line(aes(x = pop, y = persist, col = as.factor(beta*14)),
            size = 1)+
  facet_rep_grid(shedding~seasonality, labeller = labeller(seasonality = null.labs,
                                                           shedding = shedding.labs),
                 scale = "free")+
  scale_x_continuous(limits=c(1e+02, 1e+06), trans = "log10") +
  scale_y_continuous(limits=c(0,100))+
  geom_hline(yintercept = 50, lwd = 0.5, lty = 2)+
  theme_minimal()+
  scale_color_manual(name = "R0", labels = c("3.50", "7.00", "14.00"),
                     values = c("hotpink2", "navy", "olivedrab1"))+
  theme(text = element_text(size = 20))+
  xlab("population size (log10 scale)")+
  ylab(" ")+
  theme(axis.line =element_line(color="black"))+
  geom_line(data = results_mp %>% filter(shedding == 0.50),
            aes(x = 25*pop, y = persist,  col = as.factor(beta*14)), size = 1, lty = 2)+
  theme(legend.position = "right")

dummy_dat <- data.frame(x = c(1, 2, 1, 2), y = c(1, 2, 1, 3), model = c(rep("homogeneous", 2),
                                                                        rep("structured", 2)))
shed_dummy_p <- ggplot(dummy_dat)+
  geom_line(aes(x = x, y = y, lty = model),
            size = 1.5)+
  scale_linetype_manual(values = c(1, 2))+
  theme_minimal()+
  theme(text = element_text(size = 24),
        legend.position = "bottom",
        legend.key.size = unit(3.5, "cm"),
        legend.title = element_blank())

prow <- plot_grid(
  shed_1_p+theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none"),
  NULL,
  shed_50_p+theme(plot.margin = unit(c(0,0,0,0), "cm"), legend.position = "none"),
  align = 'vh',
  hjust = -1,
  nrow = 3,
  rel_heights = c(1,-0.2,1)
)

# add a shared legend
legend_b <- get_legend(
  shed_dummy_p + 
    guides(linetype = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add a shared legend for R0
legend_a <- get_legend(
  shed_1_p + 
    guides(col = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
figure_2 <- plot_grid(prow, legend_b, legend_a, ncol = 1, rel_heights = c(1, .1, .1))

ggsave(width = 10, height = 8, filename = "figs/figure_2.png", figure_2)
