# In this script we run our pruning structure learning algorithms in the 
# the extremes and very the learned structure fits our data 
# for our initial analysis we choose to generate an HR-SCM dataset where
# our HR-based pruning algorithm should give us a close to perfect fit

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

source(here("analyses","01_load_packages.R"))

#===============================================================================
# simulate HR-SCM
#===============================================================================

set.seed(123)
tau <- 0.9

b2 <- 0.4
b3 <- 0.6
n <- 1000
B <- matrix(0, 4, 4)
B[1, 2] <- 1
B[1, 3] <- 1
B[2, 4] <- b2
B[3, 4] <- b3

d <- nrow(B)

L <- matrix(0, 3, 4)
L[1, 1] <- -1
L[1,2] <- 1
L[2,1] <- -1
L[2,3] <- 1
L[3, 2] <- -b2
L[3, 3] <- -b3
L[3,4] <- 1

nu2 <- c(0.5, 0.7, 0.4)
D_eps <- diag(nu2)

# precision matrix 
Theta <- t(L) %*% solve(D_eps) %*% L
Gamma <- Theta2Gamma(Theta)
mu_eps <- -0.5 * L %*% Gamma[,1]
# simulate Gaussian noise
eps <- mvrnorm(n, mu = as.numeric(mu_eps), Sigma = D_eps)

R <- rexp(n)
root <- 1
d <- ncol(B)
Y <- matrix(0, n, d)
Y[, root] <- R
Y[,2] <- B[1,2] * Y[,1] + eps[,1]
Y[,3] <- B[1,3] * Y[,1] + eps[,2]
Y[,4] <- B[2,4] * Y[,2] + B[3,4] * Y[,3] + eps[,3]

Y_mp <- data2mpareto(Y, p = tau)

Gamma_hat <- emp_vario(Y_mp)
Gamma_hat
Gamma


#-------------------------------------------------------------------------------
# assume a fully-connected DAG based on causal order (1,2,3,4)
#-------------------------------------------------------------------------------

full_dag <- rbind(
  c(0, 1, 1, 1),  
  c(0, 0, 1, 1),  
  c(0, 0, 0, 1),  
  c(0, 0, 0, 0)
)


k <- nrow(Y_mp)
alpha <- 0.05

ci_test <- function(x, y, S, ...) {
  HR_CItest(x, y, S, suffStat = list(Gamma = Gamma_hat, n = k), ...)
}

prunbl_prs <- prunable_pairs(full_dag) %>% transpose()

B_pruned <- prune_edge_algo(
  B_full = full_dag,
  prunbl_prs = prunbl_prs,
  ci_test = ci_test,
  alpha = alpha,
  base_pruning_algo = prune_edge_fast,
  separating_set_selector = mb_selector
)

B_pruned

#-------------------------------------------------------------------------------
#  train DAG algorithm and verify inferred CIs hold on test data
#-------------------------------------------------------------------------------

set.seed(123)

# split data into train and test set
n_train <- floor(0.7 * nrow(Y_mp))
train_idx <- sample(seq_len(nrow(Y_mp)), n_train)

Y_train <- Y_mp[train_idx, ]
Y_test  <- Y_mp[-train_idx, ]

# estimate gamma on train set
Gamma_train <- emp_vario(Y_train)

# prune DAG on train set
k_train <- nrow(Y_train)
ci_test_train <- function(x, y, S, ...) {
  HR_CItest(x, y, S, suffStat = list(Gamma = Gamma_train, n = k_train), ...)
}

prunbl_prs_train <- prunable_pairs(full_dag) %>% transpose()

B_pruned_train <- prune_edge_algo(
  B_full = full_dag,
  prunbl_prs = prunbl_prs_train,
  ci_test = ci_test_train,
  alpha = alpha,
  base_pruning_algo = prune_edge_fast,
  separating_set_selector = mb_selector
)

# test fit on test data
Gamma_test <- emp_vario(Y_test)
k_test <- nrow(Y_test)

ci_test_test <- purrr::partial(
  HR_CItest,
  suffStat = list(Gamma = Gamma_test, n = k_test)
)

res_test <- perform_citests(
  B_pruned_train,
  ci_test_fun <- function(x, y, S, ...) {
    HR_CItest(x, y, S, suffStat = list(Gamma = Gamma_test, n = k_test), ...)
  }
)
res_test

# summarise results
alpha <- 0.05
res_summary <- res_test %>%
  mutate(reject = pvalue < alpha) %>%
  group_by(is_dsep) %>%
  summarise(
    n_tests = n(),
    prop_reject = mean(reject),
    .groups = "drop"
  )
res_summary


#===============================================================================
# compute conditional probabilities from extremal DAG
#===============================================================================

retrieve_parents <- function(dag){
  # convert adjacency matrix to igraph
  g_dir <- graph_from_adjacency_matrix(B_pruned, mode = "directed")
  
  parents_list <- list()
  
  # find parent nodes
  for (j in seq_along(nodes)) {
    parent_indices <- which(dag[, j] != 0)   
    parents_list[[j]] <- nodes[parent_indices]
  }
  
  names(parents_list) <- nodes
  
  # remove root node since this has no conditional probabilities
  parents_list[[nodes[1]]] = NULL
  
  return(parents_list)
}

nodes <- c("A", "B", "C", "D")

cond_density_HR <- function(dag, data, Gamma, nodes) {
  
  parents_list <- retrieve_parents(dag)
  colnames(Gamma) <- rownames(Gamma) <- nodes
  
  results_list <- list()
  colnames(data) <- nodes
  
  for (v in names(parents_list)) {
    pa <- parents_list[[v]]
    yv <- data[, v]
    ypa <- data[,pa, drop = FALSE]
    p <- length(pa)
    
    G_CC <- Gamma[pa, pa, drop = FALSE]
    M <- rbind(
      cbind(-0.5 * G_CC, matrix(1, nrow = p, ncol = 1)),
      c(rep(1, p), 0)
    )
    Minv <- solve(M)
    G_AC <- matrix(Gamma[v, pa, drop = FALSE], nrow = 1)
    
    # mu_star for all rows
    mu_star <- as.numeric(
      cbind(-0.5 * G_AC, 1) %*% Minv %*% rbind(t(ypa), 1)
    )
    
    # sd_star constant
    U <- c(v, pa)
    Theta_sub <- Gamma2Theta(Gamma[U, U])
    sd_star <- sqrt(1 / Theta_sub[1, 1])
    
    
    z_scores <- (yv - mu_star) / sd_star
    
    pit_vals <- pnorm(z_scores)
    ks_p <- ks.test(pit_vals, "punif")$p.value
    
    results_list[[v]] <- list(
      z_scores = z_scores,
      pit_vals = pit_vals,
      ks_pvalue = ks_p
    )
  }
  
  return(results_list)
}


results_sim <- cond_density_HR(
  dag   = B_pruned,   
  data = Y,       
  Gamma = Gamma_hat,   
  nodes = nodes
)

hist(results_sim[[3]]$pit_vals, breaks = 10,
     main = paste("PIT histogram for", names(results_sim)[1]),
     xlab = "PIT value", freq = FALSE)

results_sim[[1]]$ks_pvalue

pit_vals1 <- results_sim[[1]]$pit_vals
qq_df <- data.frame(
  emp_q = sort(pit_vals1),
  uni_q = qunif(ppoints(length(pit_vals1))))


p_qq <- ggplot(qq_df, aes(x = uni_q, y = emp_q)) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = paste("QQ Plot vs Uniform(0,1)"), 
       x = "Uniform Quantiles", 
       y = "Empirical Quantiles")

p_qq




