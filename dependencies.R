# rstan
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# odin.dust
library(odin.dust)

# other R packages and version used in model development
library(tidyverse) # version 2.0.0
library(readr) # version 1.4.0
library(bayestestR) # version 0.10.0
library(here) # version 1.0.1
library(emdbook) # version 1.3.12
library(binom) # version 1.1.1
# for figures
library(cowplot) # version 1.1.1
library(ggcorrplot) # version 0.1.3
library(bayesplot) # version 1.8.0
library(modeest) # version 2.4.0
library(lemon) # version 0.4.5
library(grid) # version 4.1.2
library(gridExtra) # version 2.3

# functions and fixed parameters
source("utils.R")
