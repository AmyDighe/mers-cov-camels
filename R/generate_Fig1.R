# read in fits
fit1bbexp <- readRDS("fits/processed_real/exp/sens_spec_1/fit1bb.rds")
fit2bbexp <- readRDS("fits/processed_real/exp/sens_spec_1/fit2bb.rds")
fit3bbexp <- readRDS("fits/processed_real/exp/sens_spec_1/fit3bb.rds")
fit4bbexp <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb.rds")

fit4bbexpf <- readRDS("fits/processed_real/exp/sens_spec_1/fit4bb_f.rds")
SEROPOS <- readRDS("data/real/SEROPOS.rds")
N_CAMELS <- readRDS("data/real/N_CAMELS.rds")

# Figure 1a: Model Fitting

#### beta binomial EXP
pa <- plot_fig1a(fit = fit4bbexpf,
                   data = data_sero)

ggsave(filename = "figs/1a_manuscript_4bbexp.png", 
       pa, height = 14, width = 18, units = "in")

# Figure 1b: posteriors - with mode
pb <- plot_fig1b(fit = fit4bbexpf, central_stat = "mode_ms")

ggsave(filename = "figs/1b_manuscript_mode_4bbexp_log.png", 
       pb, height = 17, width = 7, units = "in")

# Figure 1b: posteriors - with mean
pb <- plot_fig1b(fit = fit4bbexpf, central_stat = "mean")

ggsave(filename = "figs/1b_manuscript_mean_4bbexp.png", 
       pb, height = 14.5, width = 7, units = "in")
