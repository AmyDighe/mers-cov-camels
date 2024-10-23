# This script processes the raw data collected in https://pubmed.ncbi.nlm.nih.gov/31201040/ 

# it depends on:
  #  1. the raw data file in ~/data/real/sero_rna_syst.csv 
  #  2. the dplyr package
  #  3. the function "reshape_data" in the source file ~/utils.R

# it generates:
  #  1. the seroprevalence data w/ one row per age-study-group "~/data/real/data_sero.rds"
  #  2. the seroprevalence data as an age class x study padded matrix for rstan "~/data/real/SERO_POS.rds"
  #  3. the lower age bound for each element of the age class x study matrix "~/data/real/AGE_L.rds"
  #  4. the upper age bound for each element of the age class x study matrix "~/data/real/AGE_U.rds"
  #  5. the number of camels for each element of the age class x study matrix "~/data/real/N_CAMELS.rds"


# read in the data from systematic review

data_raw <- read.csv("data/real/sero_rna_syst.csv")


# select columns of interest

data_min <- data_raw%>%
  dplyr::select(STUDY, COUNTRY, REGION, LOW_AGE, UPP_AGE, N_AGE_2, RNA_POS,
                RNA_N, SERO_POS, SERO_N, TEST_TYPE,
                SAMPLE)



# filter only those studies which have serology data
# remove studies which did targeted sampling linked to human cases

data_sero <- data_min %>%
  dplyr::filter(!is.na(SERO_POS))%>%
  dplyr::filter(SAMPLE != "targetted_epilink")%>%
  dplyr::filter(STUDY != "wernery_phylo")%>% #responded to camel infection
  dplyr::filter(STUDY != "ali_cross")%>% #overlap with ali_syst
  mutate(STUDY_COUNTRY = paste(COUNTRY, STUDY, sep ="_"),
         REGION_COUNTRY_STUDY = paste(REGION, COUNTRY, STUDY, sep = "_"))


# replace any upper age limits of exactly 2 years with 1.99999
# to avoid issues with age distribution stemming from 
# assuming mortality rates change at 2 years

data_sero$UPP_AGE[which(data_sero$UPP_AGE == 2)] <- 1.99999


# reshape data for stan (cannot have ragged arrays so must pad)

dat_full <- reshape_data(data_sero = data_sero)


# save  processed data

saveRDS(dat_full$data_sero, file = "data/real/data_sero.rds")
saveRDS(dat_full$SEROPOS, file = "data/real/SEROPOS.rds")
saveRDS(dat_full$AGE_L, file = "data/real/AGE_L.rds")
saveRDS(dat_full$AGE_U, file = "data/real/AGE_U.rds")
saveRDS(dat_full$N_CAMELS, file = "data/real/N_CAMELS.rds") 

