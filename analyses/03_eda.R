# In this script, we complete an exploratory data analysis for our exchange rate
# data sets for each country

#===============================================================================
# load packages
#===============================================================================

library(here)
library(dplyr)
library(magrittr)
library(tidyr)
library(tseries)
library(ggplot2)
library(moments)
library(knitr)
library(kableExtra)
library(extRemes)

#===============================================================================
# load data
#===============================================================================

currencies <- read.csv(here("data", "derived", "currency-data.csv"))

currencies$Date <- as.Date(currencies$Date)

log_returns <- read.csv(here("data", "derived", "log-returns.csv"))
log_returns$Date <- as.Date(log_returns$Date)

garch_returns <- read.csv(here("data", "derived", "garch-returns.csv"))
garch_returns$Date <- as.Date(garch_returns$Date)

#===============================================================================
# exploratory data analysis
#===============================================================================

#-------------------------------------------------------------------------------
# visualise exchange rate for each country
#-------------------------------------------------------------------------------

data_long <- pivot_longer(currencies, 
                          cols = -Date, 
                          names_to = "Country", 
                          values_to = "Rate")

data_long$Date <- as.Date(data_long$Date)

ggplot(data_long, aes(x = Date, y = Rate)) +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  labs(title = "Exchange Rate to GBP over Time",
       x = "Date", y = "Exchange Rate") +
  geom_line() +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y")+
  theme_bw()


#-------------------------------------------------------------------------------
# calculate and visualise log returns
#-------------------------------------------------------------------------------

head(log_returns)

logret_long <- pivot_longer(log_returns, cols = -Date,
                            names_to = "Country", values_to = "LogReturn") %>%
  mutate(Country = gsub("_logret", "", Country))

logret_long <- as.data.frame(logret_long)

ggplot(logret_long, aes(x = Date, y = LogReturn, group = 1)) +
  geom_line() +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  labs(title = "Daily Log Returns of Exchange Rates",
       x = "Date", y = "Log Return") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y")+
  theme_bw()

#-------------------------------------------------------------------------------
# summary of statistics
#-------------------------------------------------------------------------------

summary_stats <- subset(log_returns, select = -Date) %>%
  pivot_longer(everything())%>%
  group_by(name) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            skewness = skewness(value),
            kurtosis = kurtosis(value))

summary_stats %>%
  kbl(caption = "Summary Statistics for Log Returns") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


garch_dat <- subset(garch_returns, select = -Date)

summary_stats_g <- garch_dat %>%
  pivot_longer(everything())%>%
  group_by(name) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            skewness = skewness(value),
            kurtosis = kurtosis(value))

summary_stats_g %>%
  kbl(caption = "Summary Statistics for Garch Standardised Residuals") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


#-------------------------------------------------------------------------------
# compare to a normal distribution
#-------------------------------------------------------------------------------

# compute mean, sd, min and max for each currency
stats <- logret_long %>%
  group_by(Country) %>%
  summarise(mean = mean(LogReturn),
            sd = sd(LogReturn),
            min = min(LogReturn),
            max = max(LogReturn))

normal_curves <- matrix(nrow = nrow(stats), ncol = 1000)

for (i in 1:nrow(stats)){
  x <- seq(stats$min[i], stats$max[i], length.out = 1000)
  normal_curves[i,] <- dnorm(x, mean = stats$mean[i], sd = stats$sd[i])
}

normal_df <- data.frame()

for (i in 1:nrow(stats)) {
  df_temp <- data.frame(
    x = seq(stats$min[i], stats$max[i], length.out = 1000),
    y = normal_curves[i, ],
    Country = stats$Country[i]
  )
  normal_df <- rbind(normal_df, df_temp)
}

logret_long <- logret_long %>%
  filter(!is.na(LogReturn)) %>%
  mutate(Country = factor(Country))


logret_only_long <- subset(logret_long, select =-Date)
str(logret_only_long$LogReturn)
summary(logret_only_long)

# plot normal density vs empirical distribution
ggplot(logret_only_long, aes(x = LogReturn)) +
  geom_histogram(aes(y= after_stat(density)), 
                 bins = 75, 
                 fill = "darkgrey", 
                 alpha = 0.6) +
  geom_density(color = "black", linewidth = 1) +
  geom_line(data = normal_df, aes(x = x, y = y), color = "red") +
  facet_wrap(~ Country, scales = "free") +
  theme_bw() +
  labs(title = "Distribution of Daily Log Returns vs Standard Normal Distribution",
       x = "Log Return",
       y = "Density")

# QQ plots
ggplot(logret_long, aes(sample = scale(LogReturn))) +
  stat_qq() +
  stat_qq_line(color = "red") +
  facet_wrap(~ Country, scales = "free") +
  labs(title = "QQ-Plots of Log Returns vs Standard Normal Distribution") +
  theme_bw()

#-------------------------------------------------------------------------------
# ADF test
#-------------------------------------------------------------------------------

log_returns_only <- colnames(subset(log_returns, select = -Date))

adf_results <- list()

for (i in log_returns_only) {
  # run ADF test on log returns for each currency
  adf_test <- adf.test(log_returns[[i]], alternative = "stationary")
  adf_results[[i]] <- adf_test
}

test_stats <- c()
p_vals <- c()

for (i in 1:length(names(adf_results))) {
  test_stats[i] <- adf_results[[i]]$statistic
  p_vals[i] <- adf_results[[i]]$p.value
}


adf_summary <- data.frame(
  Log_Returns = log_returns_only,
  Test_Statistic = test_stats,
  P_Value = p_vals,
  row.names = NULL)

# show results in a table
adf_summary %>%
  kbl(caption = "ADF Test Results for Log Returns") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


#-------------------------------------------------------------------------------
# investigate volatility clustering
#-------------------------------------------------------------------------------

abs_log_returns <- log_returns %>%
  mutate(across(-Date, ~ abs(.), .names = "{.col}_abs")) 

abs_log_returns <- abs_log_returns %>%
  dplyr::select(Date, ends_with("_abs"))

absret_long <- pivot_longer(abs_log_returns, cols = -Date,
                            names_to = "Country", values_to = "AbsLogRet") %>%
  mutate(Country = gsub("_abs", "", Country))

# plot absolute log returns
ggplot(absret_long, aes(x = Date, y = AbsLogRet, group = 1)) +
  geom_line(color = "black") +
  labs(title = "Volatility Clustering: Absolute Log Returns",
       y = "Absolute Return") +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  theme_bw()

# plot acf of absolute log returns 
par(mfrow = c(3, 2))

for (country in names(subset(abs_log_returns, select = -Date))) { 
  acf(abs_log_returns[[country]],
      main = paste(country))
}


#-------------------------------------------------------------------------------
# investigate extreme currency depreciation
#-------------------------------------------------------------------------------

right_tail_exceedances <- log_returns %>%
  summarise(across(-Date, ~sum(. > quantile(., probs = 0.95))))

print(right_tail_exceedances)

logret_long <- logret_long %>%
  group_by(Country) %>%
  mutate(threshold_95 = quantile(LogReturn, 0.95),
         is_extreme = LogReturn > threshold_95)

ggplot(logret_long, aes(x = Date, y = LogReturn)) +
  geom_line(alpha = 0.6) +
  geom_point(data = subset(logret_long, is_extreme == TRUE),
             aes(x = Date, y = LogReturn), color = "red", size = 1) +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Extreme Log Returns (Currency Depreciation Events)",
       y = "Log Return", x = "Date")


#===============================================================================
# do all plots
#===============================================================================

cze_dat <- currencies %>% select(Date, CZE)
price_p <- ggplot(cze_dat, aes(x = Date, y = CZE)) +
  geom_line(color = "black") +
  labs(title = "CZE Exchange Rate to GBP", x = "Date", y = "Exchange Rate") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

price_p

ggsave(here("outputs", "cze-prices.pdf"), price_p, width = 6, height = 4, dpi = 300)


cze_logret <- log_returns %>%
  mutate(CZE_logret = CZE * 100) %>%
  select(Date, CZE_logret)  

lr_p <- ggplot(cze_logret, aes(x = Date, y = CZE_logret)) +
  geom_line(color = "black") +
  labs(title = "CZE Daily Log Returns", x = "Date", y = "Log Return (%)") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

lr_p

ggsave(here("outputs", "cze-logret.pdf"), lr_p, width = 6, height = 4, dpi = 300)

acf_data <- acf(abs(cze_logret$CZE_logret), plot = FALSE, lag.max =30)
acf_ci <- qnorm((1 + 0.95)/2)/sqrt(length(cze_logret$CZE_logret))
acf_df <- data.frame(Lag = acf_data$lag, ACF = acf_data$acf)
acf_p <- ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "darkgrey", colour = "black") +
  labs(title = "ACF of Absolute Log Returns", x = "Lag", y = "ACF") +
  theme_bw() +
  theme(plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14))+
  geom_hline(yintercept = acf_ci, colour = "#ED8221", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -acf_ci, colour = "#ED8221", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -0, colour = "black", linetype = "dashed", linewidth = 0.4)

acf_p  

ggsave(here("outputs", "abs-cze-acf.pdf"), acf_p, width = 6, height = 4, dpi = 300)


summary_stats <- log_returns %>%
  select(-Date) %>%
  pivot_longer(everything(), names_to = "Country", values_to = "LogReturn") %>%
  group_by(Country) %>%
  summarise(
    Mean = round(mean(LogReturn), 3),
    "Standard Deviation" = round(sd(LogReturn), 3),
    Skewness = round(skewness(LogReturn), 3),
    Kurtosis = round(kurtosis(LogReturn),2)
  )

kbl(
  summary_stats,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c", 
  caption = "Summary Statistics for Daily Log Returns") %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "stats-logret.tex"))


summary_stats_g <- garch_returns %>%
  select(-Date) %>%
  pivot_longer(everything(), names_to = "Country", values_to = "Garch_resid") %>%
  group_by(Country) %>%
  summarise(
    Mean = round(mean(Garch_resid), 3),
    "Standard Deviation" = round(sd(Garch_resid), 3),
    Skewness = round(skewness(Garch_resid), 3),
    Kurtosis = round(kurtosis(Garch_resid),2)
  )

kbl(
  summary_stats_g,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c", 
  caption = "Summary Statistics GARCH Standardised Residuals") %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "stats-garchret.tex"))



g_acf_data <- acf(abs(garch_returns$CZE), plot = FALSE, lag.max =30)
g_acf_ci <- qnorm((1 + 0.95)/2)/sqrt(length(garch_returns$CZE))
g_acf_df <- data.frame(Lag = g_acf_data$lag, ACF = g_acf_data$acf)
g_acf_p <- ggplot(g_acf_df, aes(x = Lag, y = ACF)) +
  geom_bar(stat = "identity", fill = "darkgrey", colour = "black") +
  labs(title = "ACF of Absolute GARCH Returns", x = "Lag", y = "ACF") +
  theme_bw() +
  theme(plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 14))+
  geom_hline(yintercept = acf_ci, colour = "#ED8221", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -acf_ci, colour = "#ED8221", linetype = "dashed", linewidth = 0.7) +
  geom_hline(yintercept = -0, colour = "black", linetype = "dashed", linewidth = 0.4)

g_acf_p 

ggsave(here("outputs", "abs-cze-garch-acf.pdf"), g_acf_p, width = 6, height = 4, dpi = 300)

ljung_df <- data.frame(Lag = c(),
                      "$p$-value" = c(),
                      "Reject $H_0$" = c())

for (l in 1:10) {
  test <- Box.test(garch_returns$CZE, lag = l, type = "Ljung-Box")
  pval <- test$p.value
  reject <- ifelse(pval < 0.05, "Rejected", "Not Rejected")
  
  ljung_df <- rbind(ljung_df, data.frame(Lag = l,
                                         "$p$-value" = round(pval, 3),
                                         "Reject $H_0$" = reject))
}

ljung_df %>%
  kbl(format = "latex",
      booktabs = TRUE,
      escape = FALSE,
      align = "c", 
      caption = "Ljung–Box Test on CZE GARCH Standardised Residuals ($H_0$ = Independent, $H_1$ = Dependent)") %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "box-test-garch.tex"))


u_cze <- 1.6
ei_cze <- extremalindex(garch_returns$CZE, threshold = u_cze, method = "runs", run = 1)



mrl_all <- data.frame()

for (country in names(garch_dat)) {
  g_returns <- garch_dat[[country]]
  u_max <- quantile(g_returns, 0.9997)
  
  mrl_obj <- mrlplot(g_returns, nint = 10000, main = "", 
                     xlab = "Threshold", xlim = c(0, u_max))
  
  mrl_df <- as.data.frame(mrl_obj)
  colnames(mrl_df) <- c("Lower", "MRL", "Upper")
  thresholds <- seq(min(g_returns), max(g_returns) - 1, length.out = 10000)
  mrl_df$Threshold <- thresholds
  
  mrl_df <- subset(mrl_df, Threshold >= 0 & Threshold <= u_max)
  mrl_df$Country <- country   
  
  mrl_all <- rbind(mrl_all, mrl_df)
}

# faceted ggplot
mrl_p_all <- ggplot(mrl_all, aes(x = Threshold, y = MRL)) +
  geom_line(colour = "black") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey", alpha = 0.4) +
  facet_wrap(~ Country, scales = "free", ncol = 2) +
  labs(x = "Threshold", y = "Mean Excess",
       title = "Mean Residual Life Plots by Country") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 13.5),
    axis.title = element_text(size = 14)
  )

mrl_p_all
ggsave(here("outputs", "mrl-plots-all.pdf"), plot = mrl_p_all, width = 8, height = 6)


extreme_results <- data.frame()
u_values <- c(1.7, 1.6, 1.6, 1.8, 1.9, 1.8, 2)
for (i in 1:ncol(garch_dat)) {
  g_returns <- garch_dat[[colnames(garch_dat)[i]]]
  u <- u_values[i]
  
  ei <- extremalindex(g_returns, threshold = u, method = "runs", run = 1)
  extreme_results <- rbind(extreme_results, data.frame(
    Country = colnames(garch_dat)[i],
    Threshold = u,
    "Extremal Index" = round(as.numeric(ei[1]), 3)
  ))
}

kbl(
  extreme_results,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c", 
  caption = "Extremal Index Estimates for GARCH Standardised Residuals"
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "extremal-indices.tex"))


