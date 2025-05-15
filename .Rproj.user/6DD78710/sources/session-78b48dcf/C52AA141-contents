# Title:      Coverage Plots
# Author:     Jan-Lukas Wermuth
# Date:       2025-05-09
# Purpose:    This script loads the empirical coverage rates
#             and creates coverage plots.
rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(Cairo)
library(abind)

invisible(lapply(list.files(here("results/simulations/coverages/kendall"), pattern = "\\plot.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))
invisible(lapply(list.files(here("results/simulations/coverages/spearman"), pattern = "\\plot.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))
invisible(lapply(list.files(here("results/simulations/true_taus"), pattern = "\\CIs.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))
invisible(lapply(list.files(here("results/simulations/true_rhos"), pattern = "\\CIs.RData$", full.names = TRUE), function(x) load(x, envir = globalenv())))

# Parameter Specification -------------------------------------------------
SampleSizes <- c("50", "200", "800")
rhos <- c(sort(-log10(seq(1,9.9999,0.0999))[-1]), log10(seq(1,9.9999,0.0999)))[2:180]
DGPs <- c("Normal", "Cauchy", "Zipf")
size <- 10
linewidth <- 0.2

# Prepare Datasets for Plotting -------------------------------------------
for (DGP in DGPs){
  for (Ti in SampleSizes) {
    if (DGP == "Normal") {
      assign(paste("data", DGP, Ti, sep = "_"), data.frame(cbind(t(rbind(2/pi*asin(rhos), 6/pi*asin(rhos/2), colMeans(decision_kendall_array_norm_iid_plot[Ti,,]), colMeans(decision_kendall_array_norm_iid_fis_plot[Ti,,]), colMeans(decision_spearman_array_norm_iid_plot[Ti,,]), colMeans(decision_spearman_array_norm_iid_fis_plot[Ti,,]))), DGP, Ti)))
      data.table::setnames(get(paste("data", DGP, Ti, sep = "_")), c("gamma", "rhob", "Kendall_Coverage", "Kendall_Coverage_Fis", "Spearman_Coverage", "Spearman_Coverage_Fis", "DGP", "Ti"))
    }
    if (DGP == "Zipf") {
      assign(paste("data", DGP, Ti, sep = "_"), data.frame(cbind(t(rbind(gammas_MvZipf_iid_CIs[1:90], rhobs_MvZipf_iid_CIs[1:90], colMeans(decision_kendall_array_MvZipf_iid_plot[Ti,,]), colMeans(decision_kendall_array_MvZipf_iid_fis_plot[Ti,,]), colMeans(decision_spearman_array_MvZipf_iid_plot[Ti,,]), colMeans(decision_spearman_array_MvZipf_iid_fis_plot[Ti,,]))), DGP, Ti)))
      data.table::setnames(get(paste("data", DGP, Ti, sep = "_")), c("gamma", "rhob", "Kendall_Coverage", "Kendall_Coverage_Fis", "Spearman_Coverage", "Spearman_Coverage_Fis", "DGP", "Ti"))
    }
    if (DGP == "Cauchy") {
      assign(paste("data", DGP, Ti, sep = "_"), data.frame(cbind(t(rbind(taus_t1_iid_CIs[2:180], rhos_t1_iid_CIs[2:180], colMeans(decision_kendall_array_t1_iid_plot[Ti,,]), colMeans(decision_kendall_array_t1_iid_fis_plot[Ti,,]), colMeans(decision_spearman_array_t1_iid_plot[Ti,,]), colMeans(decision_spearman_array_t1_iid_fis_plot[Ti,,]))), DGP, Ti)))
      data.table::setnames(get(paste("data", DGP, Ti, sep = "_")), c("gamma", "rhob", "Kendall_Coverage", "Kendall_Coverage_Fis", "Spearman_Coverage", "Spearman_Coverage_Fis", "DGP", "Ti"))
    }
  }
}

dataset <- rbind(data_Normal_50, data_Normal_200, data_Normal_800, data_Cauchy_50, data_Cauchy_200, data_Cauchy_800, data_Zipf_50, data_Zipf_200, data_Zipf_800)
dataset$Ti <- factor(dataset$Ti, levels = c(50, 200, 800))
dataset$DGP <- factor(dataset$DGP, levels = c("Cauchy", "Normal", "Zipf"))
class(dataset$Kendall_Coverage) <- "numeric"
class(dataset$Kendall_Coverage_Fis) <- "numeric"
class(dataset$Spearman_Coverage) <- "numeric"
class(dataset$Spearman_Coverage_Fis) <- "numeric"
class(dataset$gamma) <- "numeric"
class(dataset$rhob) <- "numeric"
dataset_kendall <- pivot_longer(dataset, cols = c("Kendall_Coverage", "Kendall_Coverage_Fis"), names_to = "CI")
dataset_kendall <- dataset_kendall %>% mutate(CI = if_else(CI == "Kendall_Coverage_Fis", "with Fisher transformation", "without Fisher transformation"))
dataset_spearman <- pivot_longer(dataset, cols = c("Spearman_Coverage", "Spearman_Coverage_Fis"), names_to = "CI")
dataset_spearman <- dataset_spearman %>% mutate(CI = if_else(CI == "Spearman_Coverage_Fis", "with Fisher transformation", "without Fisher transformation"))

options(ggplot2.discrete.colour= c("blue", "red")) # Specify colors for the lines in the plots

n_names <- c(`50` = "T = 50",
             `200` = "T = 200",
             `800` = "T = 800")

fin_plot_kendall <- ggplot2::ggplot(data = dataset_kendall) +
  geom_line(mapping = aes(y = value, x = gamma, color = CI)) +
  geom_hline(yintercept=c(0.9), col="black") +
  facet_grid(cols = vars(DGP), rows = vars(Ti), labeller = labeller(Ti=n_names)) +
  theme_bw(base_size = size) + ylim(c(0,1)) + xlim(c(-0.98, 0.98)) +
  theme(axis.ticks = element_line(color = "black"), legend.position = "bottom")

ggsave(filename = here("results/plots/Kendall_NormalCauchyZipf.pdf"), plot = fin_plot_kendall, device = cairo_pdf, width = 15, height = 15, units = "cm")

fin_plot_spearman <- ggplot2::ggplot(data = dataset_spearman) +
  geom_line(mapping = aes(y = value, x = rhob, color = CI)) +
  geom_hline(yintercept=c(0.9), col="black") +
  facet_grid(cols = vars(DGP), rows = vars(Ti), labeller = labeller(Ti=n_names)) +
  theme_bw(base_size = size) + ylim(c(0,1)) + xlim(c(-0.98, 0.98)) +
  theme(axis.ticks = element_line(color = "black"), legend.position = "bottom")

ggsave(filename = here("results/plots/Spearman_NormalCauchyZipf.pdf"), plot = fin_plot_spearman, device = cairo_pdf, width = 15, height = 15, units = "cm")

