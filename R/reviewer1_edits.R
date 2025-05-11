################################################################################################
## addressing question 2 from reviewer 1 re. resolution of age class strata and effect on FoI ##
################################################################################################

# We will compare the FoI estimates from Pakistan (with 5 age classes) with 
# the estimates we would have got if there had been only 2 age classes available

source("dependencies.R")

## read in the seroprevalence data

data_sero <- readRDS("data/real/data_sero.rds") # dataframe STUDY * VARIABLES

## aggregate Pakistan, KSA, data with 4 age classes --> 2 age classes

agg_subset <- data_sero %>%
  group_by(STUDY_COUNTRY)%>%
  mutate(row_number = dplyr::row_number(),
         target = if_else((STUDY == "hemida" | STUDY == "saqib"), "Y", "N"))%>%
  filter(target == "Y")%>%
  mutate(class = if_else(row_number<3, "1", "2"))%>%
  group_by(
    REGION_COUNTRY_STUDY, STUDY_COUNTRY, STUDY, COUNTRY, REGION, TEST_TYPE, class)%>%
  summarise(SERO_POS = sum(SERO_POS),
            SERO_N = sum(SERO_N),
            LOW_AGE = min(LOW_AGE),
            UPP_AGE = max(UPP_AGE))%>%
  mutate(AGE_MID = LOW_AGE + (UPP_AGE-LOW_AGE)/2,
         seroprevalence = SERO_POS/SERO_N)

# replace the relavent datasets with their aggregated versions
data_sero <- bind_rows(data_sero %>% filter((STUDY != "hemida" & STUDY != "saqib")),
                       agg_subset)

dat_full_agg <- reshape_data(data_sero = data_sero)

## pull test sensitivity and specificity

STUDY_TEST_TYPE <- unique(data_sero %>% dplyr::select(STUDY_COUNTRY, TEST_TYPE))
STUDY_TEST_TYPE$id <- 1:nrow(STUDY_TEST_TYPE)
STUDY_TEST_TYPE_1 <- merge(STUDY_TEST_TYPE, TEST_SPEC_SENS_1)
STUDY_TEST_TYPE_1 <- STUDY_TEST_TYPE_1[order(STUDY_TEST_TYPE_1$id),]
STUDY_TEST_TYPE_2 <- merge(STUDY_TEST_TYPE, TEST_SPEC_SENS_2)
STUDY_TEST_TYPE_2 <- STUDY_TEST_TYPE_2[order(STUDY_TEST_TYPE_2$id),]

# seroconversion + mAbs + sero-reversion

fit4bb_exp_real <- stan(
  file = here::here("catalytic-stan-models/model4bb_exp.stan"),
  data = list(
    S = nrow(dat_full_agg$SEROPOS),
    A =  ncol(dat_full_agg$SEROPOS),
    N = dat_full_agg$N_CAMELS,
    pos = dat_full_agg$SEROPOS,
    age1 = dat_full_agg$AGE_L,
    age2 = dat_full_agg$AGE_U,
    sens = STUDY_TEST_TYPE_1$SENS,
    spec = STUDY_TEST_TYPE_1$SPEC,
    mabs = 1,
    mu_0 = mu_0_cat,
    mu = mu_cat
  ),
  chains = 4,
  iter = 10000,
  warmup = 5000,
  verbose = TRUE
)

saveRDS(fit4bb_exp_real, "fits/real_data/sens_spec_1/exp/fit4bb_reviewer1.rds")
fit4_summary <- rstan::summary(fit4bb_exp_real)

# pair results with studies
results <- cbind(as.data.frame(fit4_summary$summary[1:23,]), 
                 as.data.frame(dat_full_agg$SEROPOS), 
                 as.data.frame(row.names(dat_full_agg$SEROPOS)))
names(results) <- c(names(results)[1:13], "study")

