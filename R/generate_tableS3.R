# save table of FoI and global par estimates for manuscript
# Table 3.3 comparison table across the 12 models
fit1_exp1 <- readRDS("fits/processed_real/exp/sens_spec_1/fit1bb.rds")
fit2_exp1 <- readRDS("fits/processed_real/exp/sens_spec_1/fit2bb.rds")
fit3_exp1 <- readRDS("fits/processed_real/exp/sens_spec_1/fit3bb.rds")
fit4_exp1 <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds")
fit1_exp2 <- readRDS("fits/processed_real/exp/sens_spec_2/fit1bb.rds")
fit2_exp2 <- readRDS("fits/processed_real/exp/sens_spec_2/fit2bb.rds")
fit3_exp2 <- readRDS("fits/processed_real/exp/sens_spec_2/fit3bb.rds")
fit4_exp2 <- readRDS("fits/processed_real/exp/sens_spec_2/fit4bb.rds")


foi_tab <- matrix(NA, ncol = 8, nrow = 26)
fits <- list(fit1_exp1 = fit1_exp1,
             fit2_exp1 = fit2_exp1,
             fit3_exp1 = fit3_exp1,
             fit4_exp1 = fit4_exp1,
             fit1_exp2 = fit1_exp2,
             fit2_exp2 = fit2_exp2,
             fit3_exp2 = fit3_exp2,
             fit4_exp2 = fit4_exp2)

for(i in 1:8){
  fit <- fits[[i]]
  foi_tab[,i] <- extract_foi(fit)
}

colnames(foi_tab) <- names(fits)
foi_tab_df <- as.data.frame(foi_tab)
foi_tab_df$region <- c(unique(data_sero$REGION_COUNTRY_STUDY), rep("z_global", 3))
foi_tab_df$par <- c(rep("foi", 23), "sigma_r", "omega", "k")
foi_tab_df <- foi_tab_df[order(foi_tab_df$region),]

write.csv(file = "tables/tableS3_FoI.csv", foi_tab_df)
