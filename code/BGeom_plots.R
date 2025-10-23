library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)
library(Cairo)

tab.mu <- seq(0, 10, length.out = 200)
leg_title <- expression(Geom(1/(mu+1))~":")

# ---------- (a) Tie probabilities ----------
tie_geom1 <- function(mu) 1/(2*mu + 1)                  # P(X = X')
tie_geom2 <- function(mu) 1/(3*mu^2 + 3*mu + 1)         # P(X = X' = X'')

levels_a <- c("P(X=X')", "P(X=X'=X'')")

df_a <- tibble(
  mu = tab.mu,
  `P(X=X')`     = tie_geom1(tab.mu),
  `P(X=X'=X'')` = tie_geom2(tab.mu)
) |>
  pivot_longer(-mu, names_to = "measure", values_to = "value") |>
  mutate(measure = factor(measure, levels = levels_a))

p_a <- ggplot(df_a, aes(mu, value, linetype = measure, color = measure)) +
  geom_line(linewidth = 1.2) +
  scale_linetype_manual(
    values = c("P(X=X')" = "dotted", "P(X=X'=X'')" = "solid"),
    breaks = levels_a, name = leg_title
  ) +
  scale_color_manual(
    values = c("P(X=X')" = "black", "P(X=X'=X'')" = "grey40"),
    breaks = levels_a, name = leg_title
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression(mu), y = "Tie probabilities") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key.width = unit(1.3, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.8)),
         color    = "legend")   # keep color so colors appear in legend

# ---------- (b) Asymptotic variances ----------
rho_fun   <- function(mu) (1 - tie_geom2(mu))^2
tau_fun   <- function(mu) 4/9 * (1 - tie_geom2(mu))^2
gamma_fun <- function(mu) tau_fun(mu) / (1 - tie_geom1(mu))^4
taub_fun  <- function(mu) tau_fun(mu) / (1 - tie_geom1(mu))^2

levels_b <- c("rho", "tau", "gamma", "taub")

df_b <- tibble(
  mu = tab.mu,
  rho   = rho_fun(mu),
  tau   = tau_fun(mu),
  gamma = gamma_fun(mu),
  taub  = taub_fun(mu)
) |>
  pivot_longer(-mu, names_to = "measure", values_to = "value") |>
  mutate(measure = factor(measure, levels = levels_b))

p_b <- ggplot(df_b, aes(mu, value, linetype = measure, color = measure)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 1,   linewidth = 0.3) +
  geom_hline(yintercept = 4/9, linewidth = 0.3) +
  scale_linetype_manual(
    values = c(rho = "dotted", tau = "dotdash", gamma = "dashed", taub = "solid"),
    breaks = levels_b,
    labels = c(expression(rho), expression(tau), expression(gamma), expression(tau[b])),
    name = leg_title
  ) +
  scale_color_manual(
    values = c(rho = "black", tau = "grey20", gamma = "grey60", taub = "grey40"),
    breaks = levels_b, 
    labels = c(expression(rho), expression(tau), expression(gamma), expression(tau[b])),
    name = leg_title
  ) +
  coord_cartesian(ylim = c(0, 2.2)) +
  labs(x = expression(mu), y = "Asymptotic variance") +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key.width = unit(1.3, "cm"),
    legend.key.height = unit(0.45, "cm"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  ) +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.8)),
         color    = "legend")

# ---------- combine ----------
p <- (p_a | p_b)

ggsave(plot = p,
       filename = here("results/plots/BGeom_tieprob_asympVar.pdf"),
       device = cairo_pdf, width = 20, height = 10, units = "cm")
