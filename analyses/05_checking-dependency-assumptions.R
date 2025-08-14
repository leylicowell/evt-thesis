# In this script, we verify asymptotic dependence assumptions in our 
#financial data

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
library(graphicalExtremes)

#===============================================================================
# load data
#===============================================================================

source(here("analyses","01_load_packages.R"))

data <- read.csv(here("data", "derived", "garch-returns.csv"))
data <- subset(data, select = -Date)
ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(data) <- ctry_codes

dat <- as.matrix(data)
head(dat)

#===============================================================================
# verify unique root node assumption
#===============================================================================

q_vals<- seq(0.75, 0.99, length.out = 150)
B <- 20  # number of bootstrap resamples for CI estimation

# focus on root node vs other nodes to check that we don't have asymptotic 
# independence

pairs <- list(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(1, 6),
  c(1, 7)
)

asym_dep_test <- function(data, pairs, quantiles, bootstrap_num){
  n <- nrow(data)
  # initialise vectors and dataframe to store results
  q_all <- c()
  chi_all <- c()
  curr_i_all <- c()
  curr_j_all <- c()
  chi_df <- data.frame()
  
  for (pair in pairs) {
    i <- pair[1]
    j <- pair[2]
    
    for (q in quantiles) {
      chi_hat <- emp_chi(dat[, c(i, j)], p = q)[1, 2]
      chi_boot <- c()
      for (b in 1:bootstrap_num) {
        idx <- sample(seq_len(n), size = n, replace = TRUE)
        chi_boot[b] <- emp_chi(dat[idx, c(i, j)], p = q)[1, 2]
      }
      
      ci <- quantile(chi_boot, probs = c(0.025, 0.975), na.rm = TRUE)
      
      chi_df <- rbind(
        chi_df,
        data.frame(
          q = q,
          chi_hat = chi_hat,
          chi_lower = ci[1],
          chi_upper = ci[2],
          curr_i = colnames(dat)[i],
          curr_j = colnames(dat)[j]
        )
      )
      
      q_all <- c(q_all, q)
      chi_all <- c(chi_all, chi_hat)
      curr_i_all <- c(curr_i_all, colnames(dat)[i])
      curr_j_all <- c(curr_j_all, colnames(dat)[j])
    }
  }
  return(chi_df)
}


# calculate extremal correlation for different quantiles 
chi_df <- asym_dep_test(dat, pairs, q_vals, B)

# plot three different quantiles
q_df <- data.frame(
  q = c(0.85, 0.90, 0.95),
  label = c("q = 0.85", "q = 0.90", "q = 0.95")
)

# plot
ggplot(chi_df, aes(x = q, y = chi_hat)) +
  geom_ribbon(aes(ymin = chi_lower, ymax = chi_upper), 
              alpha = 0.6, 
              fill = "grey") +
  geom_line() +
  geom_vline(data = q_df,
             aes(xintercept = q, colour = label),
             linetype = "dashed") +
  facet_wrap(~ paste0(curr_i, "/", curr_j), scales = "free_y") +
  labs(x = "Quantile, q",
       y = expression(hat(chi)),
       colour = "Quantiles") +
  theme_bw()



#-------------------------------------------------------------------------------
# repeat for 4 random pairings
#-------------------------------------------------------------------------------

# for reproducibility
set.seed(42)

# sample 4 random pairs of exchange rates to plot, excluding EUR currency
d <- ncol(dat)
random_pairs <- replicate(4, sample(2:7, 2, replace = FALSE), simplify = FALSE)


random_chi_df <-  asym_dep_test(dat, random_pairs, q_vals, B)

ggplot(random_chi_df, aes(x = q, y = chi_hat)) +
  geom_ribbon(aes(ymin = chi_lower, ymax = chi_upper), 
              alpha = 0.6, 
              fill = "grey") +
  geom_line() +
  geom_vline(data = q_df,
             aes(xintercept = q, colour = label),
             linetype = "dashed") +
  facet_wrap(~ paste0(curr_i, "/", curr_j), scales = "free_y") +
  labs(x = "Quantile, q",
       y = expression(hat(chi))) +
  theme_bw()
