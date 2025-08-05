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

#===============================================================================
# load data
#===============================================================================

currencies <- read.csv(here("data", "derived", "currency-data.csv"))

currencies$Date <- as.Date(currencies$Date)

log_returns <- read.csv(here("data", "derived", "log-returns.csv"))
log_returns$Date <- as.Date(log_returns$Date)

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

ggplot(data_long, aes(x = Date, y = Rate, group = 1)) +
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
  kbl(caption = "Descriptive Statistics for Log Returns") %>%
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
  filter(!is.na(LogReturn), is.finite(LogReturn)) %>%  # remove any bad values
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
  # Run ADF test on log returns for each currency
  adf_test <- adf.test(log_returns[[i]], alternative = "stationary")
  adf_results[[i]] <- adf_test
}

adf_summary <- data.frame(
  Log_Returns = log_returns_only,
  Test_Statistic = sapply(adf_results, function(x) x$statistic),
  P_Value = sapply(adf_results, function(x) x$p.value),
  row.names = NULL)

# show results in a table
adf_summary %>%
  kbl(caption = "ADF Test Results for Log Returns") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))


#-------------------------------------------------------------------------------
# investigate volatility clustering
#-------------------------------------------------------------------------------

log_returns_squared <- log_returns %>%
  mutate(across(-Date, ~ .^2, .names = "{.col}_squared")) 

log_returns_squared <- log_returns_squared %>%
  select(Date, ends_with("_squared"))

retsquared_long <- pivot_longer(log_returns_squared, cols = -Date,
                            names_to = "Country", values_to = "LogReturnSquared") %>%
  mutate(Country = gsub("_squared", "", Country))

# plot log returns squared
ggplot(retsquared_long, aes(x = Date, y = LogReturnSquared, group = 1)) +
  geom_line(color = "black") +
  labs(title = "Volatility Clustering: Squared Log Returns",
       y = "Squared Return") +
  facet_wrap(~Country, scales = "free_y", ncol = 2) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  theme_bw()

# plot acf of log returns squared
par(mfrow = c(3, 2))

# plot ACF of standardised residuals
for (country in names(subset(log_returns_squared, select = -Date))) { 
  acf(log_returns_squared[[country]],
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


