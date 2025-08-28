# In this script we run our structure learning algorithm in the bulk 
# on gaussian data (data the model expects) and student-t data (model
# mispecified, heavy tails like our financial data)

#===============================================================================
# load packages and functions
#===============================================================================

library(here)
library(dplyr)
library(igraph)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(causalXtreme)
library(bnlearn)
library(graphicalExtremes)
library(tidyverse)
library(latex2exp)
library(rugarch)
library(dagitty)
library(MASS)  
library(Rgraphviz)
library(mvtnorm)
library(causalDisco)

source(here("analyses","01_load_packages.R"))

#===============================================================================
# simulate data function
#===============================================================================

sample_bulk_SCM <- function(n, B_adj, W, Sigma = diag(nrow(B_adj)),
                            noise = c("gaussian", "student")) {
  
  d <- nrow(B_adj)
  
  # determine causal ordering
  caus_order <- causalXtreme:::compute_caus_order(B_adj) 
  
  # simulate noise variables
  if (noise == "gaussian") {
    eps <- mvtnorm::rmvnorm(n, sigma = Sigma)
  } else if (noise == "student") {
    # student t distribution with similar kurtosis to our financial data
    eps <- mvtnorm::rmvt(n, sigma = Sigma, df = 7)
  }
  
  # store samples
  X <- matrix(0, nrow = n, ncol = d)
  
  # define root node
  X[, caus_order[1]] <- eps[, caus_order[1]]
  
  # simulate other nodes
  for (j in caus_order[-1]) {
      X[, j] <- (X %*% W[, j, drop = FALSE]) + eps[, j]
  }
  
  colnames(X) <- paste0("X", 1:d)
  return(X)
}

#===============================================================================
# simulation parameters
#===============================================================================


set.seed(123)
d_vals <- c(4, 6, 8)
n_vals <- c(3000, 5000, 10000)
n_sim <- 50 

all_results <- data.frame()


#===============================================================================
# run simulations
#===============================================================================

for (r in 1:n_sim) {
  for (d in d_vals) {
    # generate random DAG
    true_dag <- random.graph(nodes = paste0("X", 1:d))
    B_adj <- amat(true_dag)
    
    # assign weights for edges
    W <- matrix(0, d, d)
    W[B_adj == 1] <- runif(sum(B_adj), 0.5, 2)
    
    for (n in n_vals) {
      
      # simulate Gaussian data (correctly specified)
      X_gauss <- sample_bulk_SCM(n = n, B_adj = B_adj, W = W,
                                 noise = "gaussian")
      bn_gauss <- tabu(as.data.frame(X_gauss))
      
      # simulate Student-t data (misspecified)
      X_student <- sample_bulk_SCM(n = n, B_adj = B_adj, W = W,
                                   noise = "student")
      bn_student <- tabu(as.data.frame(X_student))
      
      # compute SHD
      shd_gauss   <- pcalg::shd(as.graphNEL(amat(bn_gauss)), as.graphNEL(B_adj))
      shd_student <- pcalg::shd(as.graphNEL(amat(bn_student)), as.graphNEL(B_adj))
      
      # store results
      all_results <- rbind(all_results,
                           data.frame(d = d, n = n, rep = r,
                                      model = "Gaussian",
                                      shd = shd_gauss))
      all_results <- rbind(all_results,
                           data.frame(d = d, n = n, rep = r, model = "Student-t",
                                      shd = shd_student))
      
      if (r %% 10 == 0) {
        cat("d =", d, "n =", n, "rep =", r, "SHD Gaussian =", shd_gauss,
            "SHD Student =", shd_student, "\n")
      }
    }
  }
}

#===============================================================================
# bar plot of SHD proportions
#===============================================================================

prop_df <- all_results %>%
  group_by(d, n, model, shd) %>%
  summarise(freq = n(), .groups = "drop") %>%
  group_by(d, n, model) %>%
  mutate(prop = freq / sum(freq))



shd_p <- ggplot(prop_df, aes(x = factor(shd), y = prop, fill = factor(n))) +
  geom_col(position = "dodge") +
  facet_grid(model ~ d, labeller = labeller(
    d     = function(x) paste0("d = ", x),
    model = function(x) ifelse(x == "Gaussian", "Gaussian", "Student-t")
  )) +
  scale_fill_manual(
    name = "Sample size",
    values = c("3000" = "#ED8221", 
               "5000" = "#196BB0", 
               "10000" = "#22C77B")
  ) +
  labs(x = "Structural Hamming Distance (SHD)", 
       y = "Proportion of Runs",
       title = "Distribution of SHD across Simulations") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14),
    strip.text = element_text(face = "bold", size = 13.5)
  )

shd_p

ggsave(here("outputs", "tabu-robustness.pdf"), plot = shd_p, width = 12, height = 5)
