# Title:      Power (Size) Plots
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-27
# Purpose:    This script loads the empirical size and power values
#             and creates power plots.
rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(Cairo)
library(abind)
library(here)

invisible(lapply(list.files(here("results/simulations/size_power"), pattern = "\\plot.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))
invisible(lapply(list.files(here("results/simulations/true_taus"), pattern = "iid\\.RData$|TS\\.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))

# Parameter Specification -------------------------------------------------
SampleSizes <- c("50", "200", "800")
rhos <- round(seq(-0.99, 0.99, 0.01), digits = 2)
size <- 10
linewidth <- 0.2

# IID: Prepare Datasets for Plotting -------------------------------------------
for (Ti in SampleSizes) {
  assign(paste("data", "Zipf", Ti, sep = "_"), data.frame(cbind(t(rbind(gammas_MvZipf_iid, colMeans(decision_array_zipf_iid_plot["tau",Ti,,]), colMeans(decision_array_zipf_iid_plot["r",Ti,,]), colMeans(decision_array_zipf_iid_plot["rho",Ti,,]))), "Zipf", Ti)))
  data.table::setnames(get(paste("data", "Zipf", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
  assign(paste("data", "Normal", Ti, sep = "_"), data.frame(cbind(t(rbind(2/pi*asin(rhos), colMeans(decision_array_norm_iid_plot["tau",Ti,,]), colMeans(decision_array_norm_iid_plot["r",Ti,,]), colMeans(decision_array_norm_iid_plot["rho",Ti,,]))), "Normal", Ti)))
  data.table::setnames(get(paste("data", "Normal", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
  assign(paste("data", "Cauchy", Ti, sep = "_"), data.frame(cbind(t(rbind(taus_t1_iid, colMeans(decision_array_cauchy_iid_plot["tau",Ti,,]), colMeans(decision_array_cauchy_iid_plot["r",Ti,,]), colMeans(decision_array_cauchy_iid_plot["rho",Ti,,]))), "Cauchy", Ti)))
  data.table::setnames(get(paste("data", "Cauchy", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
}

dataset <- rbind(data_Zipf_50, data_Zipf_200, data_Zipf_800, data_Cauchy_50, data_Cauchy_200, data_Cauchy_800, data_Normal_50, data_Normal_200, data_Normal_800)
dataset$Ti <- factor(dataset$Ti, levels = c(50, 200, 800))
dataset$DGP <- factor(dataset$DGP, levels = c("Cauchy", "Normal", "Zipf"))
class(dataset$kendall) <- "numeric"
class(dataset$pearson) <- "numeric"
class(dataset$spearman) <- "numeric"
class(dataset$gamma) <- "numeric"
dataset <- pivot_longer(dataset, cols = c("kendall", "pearson", "spearman"), names_to = "Power")

options(ggplot2.discrete.colour= c("blue", "red", "green4")) # Specify colors for the lines in the plots

n_names <- c(`50` = "T = 50",
             `200` = "T = 200",
             `800` = "T = 800")

fin_plot <- ggplot2::ggplot(data = dataset) +
  geom_line(mapping = aes(y = value, x = gamma, color = Power)) +
  annotate("point", x = 0, y = 0.1) +
  facet_grid(cols = vars(DGP), rows = vars(Ti), labeller = labeller(Ti=n_names)) +
  theme_bw(base_size = size) + ylim(c(0,1)) + xlim(c(-0.9, 0.9)) +
  theme(axis.ticks = element_line(color = "black"), legend.position = "bottom")

ggsave(filename = here("results/plots/PowerIID_CauchyNormZipf.pdf"), plot = fin_plot, device = cairo_pdf, width = 15, height = 15, units = "cm")

# TS: Prepare Datasets for Plotting -------------------------------------------
for (Ti in SampleSizes) {
  assign(paste("data", "Zipf", Ti, sep = "_"), data.frame(cbind(t(rbind(gammas_MvZipf_TS, colMeans(decision_array_zipf_TS_plot["tau",Ti,,]), colMeans(decision_array_zipf_TS_plot["r",Ti,,]), colMeans(decision_array_zipf_TS_plot["rho",Ti,,]))), "Zipf", Ti)))
  data.table::setnames(get(paste("data", "Zipf", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
  assign(paste("data", "Normal", Ti, sep = "_"), data.frame(cbind(t(rbind(taus_norm_TS, colMeans(decision_array_norm_TS_plot["tau",Ti,,]), colMeans(decision_array_norm_TS_plot["r",Ti,,]), colMeans(decision_array_norm_TS_plot["rho",Ti,,]))), "Normal", Ti)))
  data.table::setnames(get(paste("data", "Normal", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
  assign(paste("data", "Cauchy", Ti, sep = "_"), data.frame(cbind(t(rbind(taus_t1_TS, colMeans(decision_array_cauchy_TS_plot["tau",Ti,,]), colMeans(decision_array_cauchy_TS_plot["r",Ti,,]), colMeans(decision_array_cauchy_TS_plot["rho",Ti,,]))), "Cauchy", Ti)))
  data.table::setnames(get(paste("data", "Cauchy", Ti, sep = "_")), c("gamma", "kendall", "pearson", "spearman", "DGP", "Ti"))
}

dataset <- rbind(data_Zipf_50, data_Zipf_200, data_Zipf_800, data_Cauchy_50, data_Cauchy_200, data_Cauchy_800, data_Normal_50, data_Normal_200, data_Normal_800)
dataset$Ti <- factor(dataset$Ti, levels = c(50, 200, 800))
dataset$DGP <- factor(dataset$DGP, levels = c("Cauchy", "Normal", "Zipf"))
class(dataset$kendall) <- "numeric"
class(dataset$pearson) <- "numeric"
class(dataset$spearman) <- "numeric"
class(dataset$gamma) <- "numeric"
dataset <- pivot_longer(dataset, cols = c("kendall", "pearson", "spearman"), names_to = "Power")

options(ggplot2.discrete.colour= c("blue", "red", "green4")) # Specify colors for the lines in the plots

fin_plot <- ggplot2::ggplot(data = dataset) +
  geom_line(mapping = aes(y = value, x = gamma, color = Power)) +
  annotate("point", x = 0, y = 0.1) +
  facet_grid(cols = vars(DGP), rows = vars(Ti), labeller = labeller(Ti=n_names)) +
  theme_bw(base_size = size) + ylim(c(0,1)) + xlim(c(-0.9, 0.9)) +
  theme(axis.ticks = element_line(color = "black"), legend.position = "bottom")

ggsave(filename = here("results/plots/PowerTS_CauchyNormZipf.pdf"), plot = fin_plot, device = cairo_pdf, width = 15, height = 15, units = "cm")

