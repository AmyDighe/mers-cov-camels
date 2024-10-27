# This script estimates the relative infectiousness of reinfected camels based on viral shedding data 

# Background:
#   1. there is a big difference between the amount of viral RNA shed by naive camels 
# when infected for the first time cf. reinfected camels:
# chadox vaccine field study - https://www.nature.com/articles/s41598-019-52730-4#Sec20

#   2. there are also differences in shedding between vaccinated and unvaccinated 
#      subsequently infected animals
# chadox vaccine field study - https://www.nature.com/articles/s41598-019-52730-4#Sec20
# mva vaccine field study - https://www.science.org/doi/10.1126/science.aad1283

# data- were extracted from the papers using https://plotdigitizer.com 


# this script depends on:
#  1. the data extracted from the chadOx vaccine study ~/data/shedding/chadox.csv 
#  2. the dplyr and bayestestR packages

# it generates:
#  1. measures of the relative reduction in infectiousness in reinfected and vaccinated animals
#     which are then used to inform parameters all models in the folder: ~dynamic-odin-models/

source("dependencies.R")


##############################
## first the chadox vaccine ##
##############################

vl_data <- read.csv("data/shedding/chadox.csv")

ggplot(vl_data)+
  geom_line(aes(x = x, y = y, col = group))+
  theme_bw()

AUC_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative"))$x,
                                y = (vl_data%>%filter(group == "seronegative"))$y,
                                method = "trapezoid")
AUC_seropos <- area_under_curve(x = (vl_data%>%filter(group == "seropositive"))$x,
                                y = (vl_data%>%filter(group == "seropositive"))$y,
                                method = "trapezoid")

AUC_vax_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative_vax"))$x,
                                    y = (vl_data%>%filter(group == "seronegative_vax"))$y,
                                    method = "trapezoid")
AUC_vax_seropos<- area_under_curve(x = (vl_data%>%filter(group == "seropositive_vax"))$x,
                                   y = (vl_data%>%filter(group == "seropositive_vax"))$y,
                                   method = "trapezoid")

redshed_inf <- AUC_seropos/AUC_seroneg # seropositivity reduces shedding by 92x 
v_redshed_if_seropos <- AUC_vax_seropos/AUC_seropos # vax seropositive animals reduces their shedding by around 8x
redshed_v_pos_cf_v_neg <- AUC_vax_seropos/AUC_vax_seroneg # vax reduces shedding 1000x more in animals that are already seropositive cf seronegative animals
v_redshed_if_seroneg <- AUC_vax_seroneg/AUC_seroneg #  vaccine does not reduce shedding in seronegative animals
redshed_if_vax_and_inf <- AUC_vax_seropos/AUC_seroneg # vaccination and seropositivity reduces shedding by ~724x compared to unvaccinated naive animals
redshed_vax_neg_cf_unvax_pos <- AUC_vax_seroneg/AUC_seropos #infection alone is 126x better at reducing shed than vaccination alone

## what about log AUC values?

logAUC_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative"))$x,
                                   y = log((vl_data%>%filter(group == "seronegative"))$y),
                                   method = "trapezoid")
logAUC_seropos <- area_under_curve(x = (vl_data%>%filter(group == "seropositive"))$x,
                                   y = log((vl_data%>%filter(group == "seropositive"))$y),
                                   method = "trapezoid")

logAUC_vax_seroneg <- area_under_curve(x = (vl_data%>%filter(group == "seronegative_vax"))$x,
                                       y = log((vl_data%>%filter(group == "seronegative_vax"))$y),
                                       method = "trapezoid")
logAUC_vax_seropos<- area_under_curve(x = (vl_data%>%filter(group == "seropositive_vax"))$x,
                                      y = log((vl_data%>%filter(group == "seropositive_vax"))$y),
                                      method = "trapezoid")
logredshed_inf <- logAUC_seropos/logAUC_seroneg # seropositivity reduces shedding by 2x 
logv_redshed_if_seropos <- logAUC_vax_seropos/logAUC_seropos # vax seropositive animals reduces their shedding by around 1.5x
logredshed_v_pos_cf_v_neg <- logAUC_vax_seropos/logAUC_vax_seroneg # vax reduces shedding 3x more in animals that are already seropositive cf seronegative animals
logv_redshed_if_seroneg <- logAUC_vax_seroneg/logAUC_seroneg #  vaccine does not reduce shedding in seronegative animals
logredshed_if_vax_and_inf <- logAUC_vax_seropos/logAUC_seroneg # vaccination and seropositivity reduces shedding by ~3x compared to unvaccinated naive animals
logredshed_vax_neg_cf_unvax_pos <- logAUC_vax_seroneg/logAUC_seropos #infection alone is 2x better at reducing shed than vaccination alone


##############################
## second - the MVA vaccine ##
##############################

# this study did not include seronegative animals so we cannot directly use this 
# to parameterise the infectiousness of reinfected animals but it provides a comparison
# for the vaccination parameters

mva_rna_vax <- read.csv("data/shedding/rna_vax.csv")%>%
  group_by(dpi, group)%>%
  summarise(y = mean(y))
mva_rna_cntrl <- read.csv("data/shedding/rna_control.csv")%>%
  group_by(dpi, group)%>%
  summarise(y = mean(y))
mva_inf_vax <- read.csv("data/shedding/inf_vax.csv")%>%
  group_by(dpi, group)%>%
  summarise(y = mean(y))
mva_inf_cntrl <- read.csv("data/shedding/inf_control.csv")%>%
  group_by(dpi, group)%>%
  summarise(y = mean(y))

mva_rna <- rbind(mva_rna_vax, mva_rna_cntrl)
mva_inf <- rbind(mva_inf_vax, mva_inf_cntrl)

ggplot(mva_rna, aes(x = dpi, y = log10(y), col = group))+
  geom_line(lwd = 1)+
  theme_minimal()

ggplot(mva_inf, aes(x = dpi, y = log10(y), col = group))+
  geom_line(lwd = 1)+
  theme_minimal()

AUC_MVAcntrl_seroneg_rna <- area_under_curve(x = (mva_rna%>%filter(group == "control"))$dpi,
                                         y = (mva_rna%>%filter(group == "control"))$y,
                                         method = "trapezoid")

AUC_MVAvax_seroneg_rna <- area_under_curve(x = (mva_rna%>%filter(group == "vax"))$dpi,
                                y = (mva_rna%>%filter(group == "vax"))$y,
                                method = "trapezoid")

AUC_MVAcntrl_seroneg_inf <- area_under_curve(x = (mva_inf%>%filter(group == "control"))$dpi,
                                             y = (mva_inf%>%filter(group == "control"))$y,
                                             method = "trapezoid")

AUC_MVAvax_seroneg_inf <- area_under_curve(x = (mva_inf%>%filter(group == "vax"))$dpi,
                                           y = (mva_inf%>%filter(group == "vax"))$y,
                                           method = "trapezoid")

# relative infectiousness of vaccinated animals
100*AUC_MVAvax_seroneg_inf/AUC_MVAcntrl_seroneg_inf # ~2%
100*AUC_MVAvax_seroneg_rna/AUC_MVAcntrl_seroneg_rna # ~7%


