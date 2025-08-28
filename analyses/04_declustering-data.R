# In this script, we will test the effects of declustering vs taking standardised 
# garch residuals and the resulting implications for our GPD fits

#===============================================================================
# load packages
#===============================================================================

library(here)
library(dplyr)
library(magrittr)
library(tidyr)
library(tseries)
library(ggplot2)
library(knitr)
library(kableExtra)
library(tidyverse)
library(rugarch)
library(bayesplot)
library(cmdstanr)
library(extRemes)
library(ismev)
library(data.table)
library(evir)
library(knitr)

#===============================================================================
# load data
#===============================================================================

log_returns <- read.csv(here("data", "derived", "log-returns.csv"))
log_returns$Date <- as.Date(log_returns$Date)

garch_returns <- read.csv(here("data", "derived", "garch-returns.csv"))
garch_returns$Date <- as.Date(garch_returns$Date)

# choose a country
garch_series <- garch_returns$CZE
returns <- log_returns$CZE


#===============================================================================
# decluster
#===============================================================================

#-------------------------------------------------------------------------------
# compute extremal index for multiple r and u values and check for stability 
#-------------------------------------------------------------------------------

plot_ext_index <- function(data, r_vals){
  u_values <- seq(quantile(data, probs = 0.75), 
                  quantile(data, probs = 0.99), 
                  length.out=10)
  
  for (r in r_vals){
    
    theta_estimates <- rep(0, length(u_values))
    ci_lower <- rep(0,length(u_values))
    ci_upper <- rep(0, length(u_values))
    
    for (i in 1:length(u_values)){
      res <- extremalindex(data, 
                           threshold = u_values[i], 
                           method = "runs", 
                           run = r, 
                           plot = FALSE)
      theta_estimates[i] <- res[1]
      ci_bounds <- ci(res)
      ci_lower[i] <- ci_bounds[1]
      ci_upper[i] <- ci_bounds[5]
    }
    
    df <- data.frame(u = u_values, 
                     theta = theta_estimates, 
                     ci_lower = ci_lower, 
                     ci_upper = ci_upper)
    
    
    plot <- ggplot(df, aes(x = u, y = theta)) +
      geom_line(color = "black") +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), 
                  fill = "grey", alpha = 0.5) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "#ED8221") +
      labs(
        title = paste("Run Length =", r),
        x = "Threshold", 
        y = "Extremal Index"
      ) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14)
      )
    
    print(plot)
    # save each plot individually
    ggsave(
      filename = here("outputs", paste0("extremal_index_", r, ".pdf")),
      plot = plot,
      width = 6, height = 4, dpi = 300)
  }
}


returns <- log_returns$CZE*100
r_vals <- c(5:8)
#plot_ext_index(returns, r_vals)


# pick r = 6 and u = 0.011 to have good bias-variance tradeoff
u <- 0.011*100
m <- 6
res <- extremalindex(returns, 
                     threshold = u, 
                     method = "runs", run = m, plot = FALSE)
print(res)

#-------------------------------------------------------------------------------
# decluster data and compute maxima of each cluster
#-------------------------------------------------------------------------------

declustered <- extRemes::decluster(returns, threshold = u, r = m)

print(declustered)

# calculate cluster maxima
cluster_data <- data.frame(
  exceedance = returns[returns > u],
  cluster_id = attributes(declustered)$clusters
)

# Compute cluster maxima
cluster_maxima <- cluster_data %>%
  group_by(cluster_id) %>%
  summarise(cluster_max = max(exceedance))

cluster_max <- cluster_maxima$cluster_max

#===============================================================================
# fit GPD
#===============================================================================

u <- 0.011*100
excesses <- cluster_max - u
avg_cluster_num <- max(attributes(declustered)$clusters)/ (length(returns)/255)

fit_declust <- gpd.fit(cluster_max, threshold = u, npy = avg_cluster_num)
gpd.diag(fit_declust)


#-------------------------------------------------------------------------------
# plot cluster maxima
#-------------------------------------------------------------------------------
u <- 1.1
df <- data.frame(
  Date = log_returns$Date,
  Returns = returns,
  Declustered = FALSE,
  Maxima = FALSE
)

# add exceedances and mark these TRUE
declust_idx <- which(returns > u)
df$Declustered[declust_idx] <- TRUE

# retrieve maxima of cluster indices
maxima_idx <- which(returns %in% cluster_max)
df$Maxima[maxima_idx] <- TRUE

df <- df %>%
  mutate(Type = case_when(
    Maxima ~ "Cluster maxima",
    Declustered ~ "Declustered exceedances",
    TRUE ~ "Returns below threshold")) 

cluster_p <- ggplot(df, aes(x = Date, y = Returns, colour = Type, alpha = Type, size = Type)) +
  geom_point() +
  geom_hline(yintercept = u, linetype = "dashed", colour = "#ED8221", linewidth = 0.8) +
  scale_colour_manual(
    name = "Return Type",
    values = c(
      "Returns below threshold" = "black",
      "Declustered exceedances" = "#22C77B", 
      "Cluster maxima" = "#196BB0"         
    )) +
  scale_size_manual(
    values = c(
      "Returns below threshold" = 1,        
      "Declustered exceedances" = 1.7,
      "Cluster maxima" = 1.7      
    )) +
  scale_alpha_manual(
    values = c(
      "Returns below threshold" = 0.8,    
      "Declustered exceedances" = 0.9,
      "Cluster maxima" = 0.8
    ))+
  guides(alpha = "none")+
  guides(size = "none")+
  labs(x = "Date", y = "Daily Log Returns (%)", title = "Declustered Exceedances and Cluster Maxima") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14))

cluster_p

ggsave(here("outputs", "clusters-plot.pdf"), plot = cluster_p, width = 8, height = 4)


#-------------------------------------------------------------------------------
# m-observation return level
#-------------------------------------------------------------------------------

xi <- fit_declust$mle[2]
sigma <- fit_declust$mle[1]

se_sigma <- fit_declust$se[1]
se_xi <- fit_declust$se[2]
cov <- fit_declust$cov

# total observations and exceedances
n_total <- length(returns)
n_u <- sum(returns > u)
zeta <- n_u / n_total  

# estimated extremal index
theta <- as.numeric(res["extremal.index"])


ny <- 255  # approx. number of trading days per year
N_years <- seq(0.1, 1000, by=1) # return periods 
m_obs <- N_years * ny        

# calculate return levels
return_level <- function(m) {
  if (abs(xi) > 1e-6) {
    u + (sigma/xi) * (((m * zeta * theta)^xi) - 1)
  } else {
    u + sigma * log(m * zeta * theta)
  }
}

rl_vals <- c()
for (i in 1:length(m_obs)){
  rl_vals[i] <- return_level(m_obs[i])
}


# bootstrap for CIs
set.seed(123)
B <- 1000 

n <- length(cluster_max)
boot_res_t <- matrix(NA, nrow = B, ncol = length(m_obs))

for (b in 1:B) {
  idx <- sample(1:n, size = n, replace = TRUE)
  boot_data <- cluster_max[idx]
  
  # fit GPD
  fit_b <- gpd.fit(boot_data, threshold = u, show = FALSE)
  sigma_b <- fit_b$mle[1]
  xi_b <- fit_b$mle[2]
  
  # calculate return levels for each m in m_obs
  for (j in 1:length(m_obs)) {
    m <- m_obs[j]
    if (abs(xi_b) > 1e-6) {
      boot_res_t[b, j] <- u + (sigma_b / xi_b) * (((m * zeta * theta)^xi_b) - 1)
    } else {
      boot_res_t[b, j] <- u + sigma_b * log(m * zeta * theta)
    }
  }
}

# compute CIs 
ci_low <- c()
ci_high <- c()

for (j in 1:length(m_obs)) {
  ci_low[j] <- quantile(boot_res_t[, j], 0.025)
  ci_high[j] <- quantile(boot_res_t[, j], 0.975)
}


plot_df <- data.frame(
  m_obs = m_obs,
  ReturnLevel = rl_vals,
  CI_low = ci_low,
  CI_high = ci_high
)

# calculate empirical return levels
empirical_rl_points <- function(data, ny) {
  n <- length(data)
  sorted_data <- sort(data, decreasing = TRUE)
  
  k <- 1:n
  T_years <- (n / k) / ny 
  
  data.frame(
    ReturnPeriod = T_years,
    EmpiricalRL  = sorted_data
  )
}

# declustered empirical return levels
emp_declust_df <- empirical_rl_points(cluster_max, ny = avg_cluster_num)
# plot
p <- ggplot(plot_df, aes(x = m_obs/ny, y = ReturnLevel)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", linewidth = 1) +
  geom_point(data = emp_declust_df, aes(x = ReturnPeriod, y = EmpiricalRL), colour = "#22C77B") +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() +
  theme(plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14)) +
  labs(x = "Return period (years)", y = "Return level", title = "Return Level Plot for Declustered Data ")
p

ggsave(here("outputs", "rl-declustered.pdf"), plot = p, width = 6, height = 4)


#===============================================================================
# fit GPD to garch residuals
#===============================================================================

par(mfrow = c(1,1))
garch_dat <- subset(garch_returns, select = -Date)

g_returns <- garch_dat$CZE

#for (i in 1:ncol(garch_dat)) {
#  series <- garch_dat[[i]]
#  
#  q_levels <- c(0.90, 0.95, 0.99)
#  q_vals <- quantile(series, probs = q_levels)
#  
#  mrl.plot(series, umin = 0, umax = quantile(series, probs = 0.995))
#  
#  cols <- c("blue", "darkgreen", "red")
#  for (j in seq_along(q_vals)) {
#    abline(v = q_vals[j], col = cols[j], lty = 2)
#  }
#  legend("topleft", legend = paste0(q_levels * 100, "%"),
#         col = cols, lty = 2, cex = 0.8, bty = "n")
#  threshrange.plot(series, 
#                   r = c(quantile(series, probs = 0.75), 
#                         quantile(series, probs = 0.99)), 
#                   nint=10)
#}
#
u_vals <- c(1.5, 1.5, 1.6, 1.7, 1.6, 1.5, 1.5)

gpd_fits <- list()
for (i in 1:ncol(garch_dat)){
 print(colnames(garch_dat)[i])
 gpd_fits[[i]] <-  gpd.fit(garch_dat[,i], u_vals[i] , npy=255)
}


ctry_codes <- colnames(garch_dat)

sample_size <- c()
sigma <- c()
sigma_se <- c()
xi <- c()
xi_se <-  c()

for (i in seq_along(ctry_codes)) {
  fit <- gpd_fits[[i]]
  
  sample_size[i] <- sum(garch_dat[,i] > u_vals[i])
  
  sigma[i]    <- round(fit$mle[1], 3)
  sigma_se[i] <- round(fit$se[1], 3)
  
  xi[i]    <- round(fit$mle[2], 3)
  xi_se[i] <- round(fit$se[2], 3)
}

table_df <- data.frame(
  Country = ctry_codes,
  Sample_size = sample_size,
  Sigma = sigma,
  Sigma_SE = sigma_se,
  Xi = xi,
  Xi_SE = xi_se
)

kbl(
  table_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "GPD parameter estimates for each country (SE = standard error).",
  col.names = c("Country", "Sample Size", "$\\sigma$", "SE($\\sigma$)", "$\\xi$", "SE($\\xi$)")
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "gpd-params-all-countries.tex"))



# GPD fit to residuals for one dataset
# use mrlplot here since it returns mrl values
g_returns <- garch_dat$CZE

u_max <- quantile(g_returns, 0.9997)

mrl_obj <- mrlplot(g_returns, nint = 10000, main = "", xlab = "Threshold",xlim = c(0,u_max))

mrl_df <- as.data.frame(mrl_obj)
colnames(mrl_df) <- c("Lower", "MRL", "Upper")
thresholds <- seq(min(g_returns), max(g_returns) - 1, length.out = 10000)
mrl_df$Threshold <- thresholds

mrl_df <- subset(mrl_df, Threshold >= 0 & Threshold <= u_max)

mrl_p <- ggplot(mrl_df, aes(x = Threshold, y = MRL)) +
  geom_line(colour = "black") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey", alpha = 0.4) +
  labs(
    title = "CZE Mean Residual Life Plot",
    x = "Threshold", y = "Mean Excess"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

mrl_p
ggsave(here("outputs", "mrl-plot.pdf"), plot = mrl_p, width = 6, height = 4)

thresh_data <- threshrange.plot(
  g_returns,
  r = c(0, quantile(g_returns, 0.99)),
  type = "GP",
  nint = 15,
  alpha = 0.05,
  set.panels = FALSE
)

thresh_df <- as.data.frame(thresh_data)
thresh_df$Threshold = seq(0, quantile(g_returns, 0.99), length.out = nrow(thresh_df))

colnames(thresh_df) <- c("Sigma_lower", "Xi_lower", "Sigma", "Xi", "Sigma_upper",
                         "Xi_upper", "Threshold")

# sigma plot
sigma_plot <- ggplot(thresh_df, aes(x = Threshold, y = Sigma)) +
  geom_line(colour = "black") +
  geom_point(shape = 1, size = 2, colour = "black") +
  geom_errorbar(aes(ymin = Sigma_lower, ymax = Sigma_upper), width = 0.05) +
  labs(
    title = expression("Reparameterised Scale Threshold Stability Plot"),
    x = "Threshold", y = paste0("Reparameterised scale")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

sigma_plot

# xi plot
xi_plot <- ggplot(thresh_df, aes(x = Threshold, y = Xi)) +
  geom_line(colour = "black") +
  geom_point(shape = 1, size = 2, colour = "black") +
  geom_errorbar(aes(ymin = Xi_lower, ymax = Xi_upper), width = 0.05) +
  labs(
    title = expression("Shape Threshold Stability Plot"),
    x = "Threshold", y = paste0("Shape")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

xi_plot

# save plots
ggsave(here("outputs", "cze-scale-plot.pdf"), sigma_plot, width = 6, height = 4, dpi = 300)
ggsave(here("outputs", "cze-xi-plot.pdf"), xi_plot, width = 6, height = 4, dpi = 300)



u <- 1.6
fit_resid <- gpd.fit(garch_dat$CZE, threshold = u, npy =255)
gpd.diag(fit_resid)
sigma_hat <- fit_resid$mle[1]
xi_hat    <- fit_resid$mle[2]

k <- sum(garch_dat$CZE > u)
zeta_u <- k / n
ny <- 255
N_years <- seq(0.1, 1000, by=1)
m_obs <- N_years * ny

# calculate return levels
return_level_resid <- function(m, sigma, xi) {
  if (abs(xi) > 1e-6) {
    u + (sigma / xi) * (((m * zeta_u)^xi) - 1)
  } else {
    u + sigma * log(m * zeta_u)
  }
}

rl_vals <- numeric(length(m_obs))
for (j in seq_along(m_obs)) {
  rl_vals[j] <- return_level_resid(m_obs[j], sigma_hat, xi_hat)
}

# bootstrap for CIs
set.seed(123)
B <- 1000
n <- length(garch_dat$CZE)
boot_res_t <- matrix(NA, nrow = B, ncol = length(m_obs))

for (b in 1:B) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  boot_data <- garch_dat$CZE[idx]
  
  # fit GPD
  fit_b <- gpd.fit(boot_data, threshold = u, show = FALSE)
  sigma_b <- fit_b$mle[1]
  xi_b <- fit_b$mle[2]
  
  # calculate return levels for this bootstrap
  for (j in 1:length(m_obs)) {
    boot_res_t[b, j] <- return_level_resid(m_obs[j], sigma_b, xi_b)
  }
}

# compute CIs 
ci_low <- c()
ci_high <- c()
for (j in 1:length(m_obs)) {
  ci_low[j] <- quantile(boot_res_t[, j], 0.025)
  ci_high[j] <- quantile(boot_res_t[, j], 0.975)
}


plot_df_resid <- data.frame(
  m_obs = m_obs,
  ReturnLevel = rl_vals,
  CI_low = ci_low,
  CI_high = ci_high
)

empirical_rl_exceed <- function(data, u, ny) {
  exceedances <- data[data > u] 
  Nu <- length(exceedances)
  n  <- length(data)
  
  sorted_exc <- sort(exceedances, decreasing = TRUE)
  k <- 1:Nu

  T_years <- (n / k) / ny
  
  data.frame(
    ReturnPeriod = T_years,
    EmpiricalRL  = sorted_exc
  )
}

emp_resid_df <- empirical_rl_exceed(garch_dat$CZE,u = u, ny = 255)

p1 <- ggplot(plot_df_resid, aes(x = m_obs/ny, y = ReturnLevel)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), fill = "grey80", alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  geom_point(data = emp_resid_df, aes(x = ReturnPeriod, y = EmpiricalRL), colour = "#22C77B") +
  scale_x_log10(labels = scales::label_number()) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  ) +
  labs(
    x = "Return period (years)",
    y = "Return level",
    title = "Return Level Plot for GARCH Residuals"
  )

p1
ggsave(here("outputs", "rl-garch-residuals.pdf"), plot = p1, width = 6, height = 4)

#===============================================================================
# table summarising GPD fits for garch vs declustered
#===============================================================================


u_declust <- 1.1
u_garch <- 1.6

garch_exc <- garch_dat$CZE[garch_dat$CZE>u_garch]

# extract MLEs + SEs and calculate RSEs
table_df <- data.frame(
  Method = c("Declustering", "GARCH"),
  Sample_size = c(length(cluster_max), length(garch_exc)),
  Sigma = round(c(fit_declust$mle[1], fit_resid$mle[1]), 3),
  Sigma_SE = round(c(fit_declust$se[1], fit_resid$se[1]), 3),
  Sigma_RSE = c(
    paste0(round(fit_declust$se[1] / fit_declust$mle[1] * 100, 1), "\\%"),
    paste0(round(fit_resid$se[1] / fit_resid$mle[1] * 100, 1), "\\%")
  ),
  
  Xi = round(c(fit_declust$mle[2], fit_resid$mle[2]), 3),
  Xi_SE = round(c(fit_declust$se[2], fit_resid$se[2]), 3),
  Xi_RSE = c(
    paste0(round(fit_declust$se[2] / fit_declust$mle[2] * 100, 1), "\\%"),
    paste0(round(fit_resid$se[2] / fit_resid$mle[2] * 100, 1), "\\%")
  )
)


# table
kbl(
  table_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c", 
  caption = "GPD parameter estimates for declustered data 
  and GARCH standardised residuals (SE = standard error, 
  RSE = relative standard error).",
  col.names = c("Method", "Sample Size", "$\\sigma$", "SE($\\sigma$)", 
                "RSE($\\sigma$)", "$\\xi$", "SE($\\xi$)", "RSE($\\xi$)")
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "gpd-params.tex"))



#===============================================================================
# plot declustered vs garch fit comparisons
#===============================================================================



unif_cond <- function(x, u, xi, sigma) {
  if (abs(xi) > 1e-6) {
    1 - (1 + xi * (x - u)/sigma)^(-1/xi)
  } else {
    1 - exp(-(x - u)/sigma)
  }
}

# declustered
u_declust <- 1.1
excess_declust <- cluster_max
fit_declust <- gpd.fit(cluster_max, threshold = u_declust, npy = avg_cluster_num)
sigma_declust <- fit_declust$mle[1]
xi_declust    <- fit_declust$mle[2]
u_common <- 1.6
excess_declust_common <- cluster_max[cluster_max > u_common]
sigma_declust_common <- sigma_declust + xi_declust * (u_common - u_declust)
unif_declust <- unif_cond(excess_declust_common, u_common,
                          xi_declust, sigma_declust_common)
# GARCH residuals
u_garch <- 1.6
excess_garch <- garch_dat$CZE[garch_dat$CZE > u_garch]
fit_resid <- gpd.fit(garch_dat$CZE, threshold = u_garch, npy =255)
sigma_garch <- fit_resid$mle[1]
xi_garch    <- fit_resid$mle[2]
unif_garch <- unif_cond(excess_garch, u_garch, xi_garch, sigma_garch)


get_unif_rl <- function(unif_values, B = 500, label = "Method") {
  m_obs   <- seq(1, 1000, by = 1)
  probs   <- 1 - 1/m_obs
  
  rl_hat <- quantile(unif_values, probs)
  
  # bootstrap CIs
  set.seed(123)
  n <- length(unif_values)
  boot_mat <- replicate(B, {
    U_boot <- sample(unif_values, n, replace = TRUE)
    quantile(U_boot, probs)
  })
  ci_low  <- apply(boot_mat, 1, quantile, 0.025)
  ci_high <- apply(boot_mat, 1, quantile, 0.975)
  
  data.frame(
    ReturnPeriod = m_obs,
    UniformRL = rl_hat,
    CI_low = ci_low,
    CI_high = ci_high,
    Method = label
  )
}

df_declust <- get_unif_rl(unif_declust, B = 500, label = "Declustered")
df_garch   <- get_unif_rl(unif_garch, B = 500, label = "GARCH residuals")

df_all <- rbind(df_declust, df_garch)


comp_plot <- ggplot(df_all, aes(x = ReturnPeriod, y = UniformRL, colour = Method, fill = Method)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2) +
  geom_line(aes(y = 1 - 1/ReturnPeriod), linetype = "dashed", colour = "black",
            linewidth = 1) +
  scale_x_log10() +
  scale_colour_manual(values = c("Declustered" = "#196BB0", 
                                 "GARCH residuals" = "#22C77B")) +
  scale_fill_manual(values   = c("Declustered" = "#196BB0", 
                                 "GARCH residuals" = "#22C77B"))+
  labs(x = "Return period (number of exceedances)",
       y = "Conditional uniform return level",
       title = "Conditional Uniform Return Level Comparison") +
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14))

comp_plot
ggsave(here("outputs", "comparison-plot.pdf"), comp_plot, width = 8, height = 4)


qq_with_ci <- function(u_vals, B = 500) {
  n <- length(u_vals)
  probs <- ppoints(n)        
  theor_q <- probs          
  # empirical quantiles
  emp_q <- quantile(u_vals, probs)
  
  # bootstrap quantiles
  set.seed(123)
  boot_mat <- replicate(B, {
    U_b <- sample(u_vals, n, replace = TRUE)
    quantile(U_b, probs)
  })
  
  ci_low  <- apply(boot_mat, 1, quantile, 0.025)
  ci_high <- apply(boot_mat, 1, quantile, 0.975)
  
  data.frame(
    theor_q = theor_q,
    emp_q = as.numeric(emp_q),
    ci_low = ci_low,
    ci_high = ci_high
  )
}


df_garch_q <- qq_with_ci(unif_garch)

ggplot(df_garch_q, aes(x = theor_q, y = emp_q)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "#196BB0") +
  geom_point(size = 1, colour = "#22C77B") +
  labs(x = "Uniform Quantiles",
       y = "Empirical Quantiles",
       title = "QQ-plot with 95% Confidence Intervals (GARCH)") +
  theme_minimal()


df_declust_q <- qq_with_ci(unif_declust)

ggplot(df_declust_q, aes(x = theor_q, y = emp_q)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linetype = "dashed") +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), alpha = 0.2, fill = "#196BB0") +
  geom_point(size = 1, colour = "#22C77B") +
  labs(x = "Uniform Quantiles",
       y = "Empirical Quantiles",
       title = "QQ-plot with 95% Confidence Intervals (Declustered)") +
  theme_minimal()



