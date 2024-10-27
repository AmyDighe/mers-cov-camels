# generates figure S2 & S3

###########################################################
# Figure S2 - robustness of model ranking to LOO analysis #
###########################################################

# comparing DIC
looDIC4bb <- readRDS("fits/real_data/sens_spec_1/exp/LOO/4bb_fit.rds")
DIC_4bb <- rep(NA, 23)
for(i in 1:23){
  DIC_4bb[i] <- DIC(fit = looDIC4bb[[i]])$DIC_mode
}

looDIC3bb <- readRDS("fits/real_data/sens_spec_1/exp/LOO/3bb_fit.rds")
DIC_3bb <- rep(NA, 23)
for(i in 1:23){
  DIC_3bb[i] <- DIC(fit = looDIC3bb[[i]])$DIC_mode
}

looDIC2bb <- readRDS("fits/real_data/sens_spec_1/exp/LOO/2bb_fit.rds")
DIC_2bb <- rep(NA, 23)
for(i in 1:23){
  DIC_2bb[i] <- DIC(fit = looDIC2bb[[i]])$DIC_mode
}

looDIC1bb <- readRDS("fits/real_data/sens_spec_1/exp/LOO/1bb_fit.rds")
DIC_1bb <- rep(NA, 23)
for(i in 1:23){
  DIC_1bb[i] <- DIC(fit = looDIC1bb[[i]])$DIC_mode
}

# pull together into a single dataframe
DIC_df <- data.frame("1bb" = DIC_1bb, "2bb" = DIC_2bb, "3bb" = DIC_3bb, "4bb" = DIC_4bb)
saveRDS(file = "fits/real_data/sens_spec_1/exp/LOO/DIC_exp.rds", DIC_df)

DIC_df <- readRDS("fits/real_data/sens_spec_1/exp/LOO/DIC_exp.rds")
names(DIC_df)<- c("1", "2", "3", "4")
DIC_df$LOO <- unique(data_sero$STUDY_COUNTRY)
DIC_long <- pivot_longer(DIC_df, 1:4, names_to = "model", 
                         values_to = "DIC")

# plot DIC by model for each set of datasets
ggplot(data = DIC_long)+
  geom_point(aes(x = LOO, y = DIC, col = model), size = 3)+
  scale_color_manual(values= c("cornflowerblue", "#355C7D", "#C06C84", "#F67280", "#F8B185"))+
  scale_x_discrete(labels = pretty_names)+
  theme_minimal()+
  xlab("data-set left out")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 25))#+
ggsave(filename = "figs/figS2.png")

#################################################################
# Figure S3 - robustness of parameter estimates to LOO analysis #
#################################################################

cv4bb <- readRDS("fits/processed_real/exp/sens_spec_1/fit4_loo.rds")
cv4bb_df <- do.call(rbind.data.frame, cv4bb) 
cv4bb_df <- cv4bb_df[grepl("sigma|k" , rownames(cv4bb_df)) & !grepl("log_lik", rownames(cv4bb_df)),
                     c("mode_ms", "2.5%", "97.5%")] %>%
  dplyr::rename(low = "2.5%", upp = "97.5%")
cv4bb_df$set <- rep(unique(data_sero$STUDY_COUNTRY), each = 3)
cv4bb_df$par <- rep(c("sigma_r", "sigma_m", "k"), times = 23)

## plot 4bb
par.labs4 <- c("k", "omega", "sigma")
names(par.labs4) <- c("k", "sigma_m", "sigma_r")
ggplot(cv4bb_df)+
  geom_point(aes(x = set, y = mode_ms))+
  geom_errorbar(aes(ymin = low, ymax = upp, x = set))+
  facet_wrap(~par, scale = "free_y", nrow = 3, labeller = labeller(par = par.labs4))+
  theme_minimal()+
  scale_x_discrete(labels = pretty_names)+
  xlab("dataset left out")+
  ylab("value")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 20))+
  geom_hline(data = main4bb_df, aes(yintercept = value), lty = 2, col = "hotpink", lwd = 1)

ggsave(filename = "figs/figS3.png")
