# In this script we run our pruning structure learning algorithms in the 
# the extremes and very the learned structure fits our data 
# for this analysis, we choose to generate an SCM data set 

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
# simulate SCM
#===============================================================================

set.seed(42)
d <- 4
tau <- 0.975
alpha = 0.05
nodes <- c("A", "B", "C", "D")
n <- 10000
Sigma <- diag(c(1, 3, 1, 2))

# define true extreme DAG
B <- rbind(
  c(0, 1, 1, 0),
  c(0, 0, 0, 1),
  c(0, 0, 0, 0),
  c(0, 0, 0, 0)
)


B_0_w <- B
B_0 <- B
B_resit <- B
B_resit[3, 4] <- 1

Y <- sample_SCM(n = n,
                B_full = B_resit, 
                B_0 = B_0, 
                B_0_w = B_0_w, 
                Sigma = Sigma)

prunbl_prs <- prunable_pairs(B_resit) %>%
  arrange(desc(y)) %>% transpose()

  
Y_mp <- graphicalExtremes::data2mpareto(Y, p = tau)
  
# estimate variogram and HR CI-test
Gamma_hat <- emp_vario(Y_mp)
ci_test_ext <- function(x, y, S, ...) {
  HR_CItest(x, y, S, suffStat = list(Gamma = Gamma_hat, n = nrow(Y_mp)), ...)
}
  
# estimate extremal DAG via pruning
B_pruned <- prune_edge_algo(
  B_full = B_resit,
  prunbl_prs = prunbl_prs,
  ci_test = ci_test_ext,
  alpha = alpha,
  base_pruning_algo = purrr::partial(prune_edge_fast, verbose = FALSE),
  separating_set_selector = mb_selector,
  shuffle = TRUE
)
  
colnames(B_pruned) <- nodes
rownames(B_pruned) <- nodes

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

prunbl_prs_train <- prunable_pairs(B_resit) %>% transpose()

B_pruned_train <- prune_edge_algo(
  B_full = B_resit,
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
  ci_test_fun = ci_test_test
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
# compute conditional probabilities from inferred DAG
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
      ks_pvalue = ks_p,
      mu = mu_star,
      sigma = sd_star,
      yv = yv,
      ypa = ypa
    )
  }
  
  return(results_list)
}


results_sim <- cond_density_HR(
  dag   = B_pruned,   
  data = Y_mp,       
  Gamma = Gamma_hat,   
  nodes = nodes
)

hist(results_sim[[2]]$pit_vals, breaks = 10,
     main = paste("PIT histogram for", names(results_sim)[1]),
     xlab = "PIT value", freq = FALSE)

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


