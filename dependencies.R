# rstan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# odin.dust
library(odin.dust)

# other R packages
library(tidyr)
library(dplyr)
library(readr)
library(bayestestR)
library(here)
library(emdbook)
library(cowplot)
library(ggcorrplot)
library(bayesplot)
library(modeest)
library(lemon)
library(grid)
library(gridExtra)
# functions and fixed parameters
source("utils.R")
