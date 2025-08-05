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
library(timeSeries)
library(fGarch)
library(tidyverse)
library(latex2exp)

source(here("analyses","load_packages.R"))

data <- read.csv(here("data", "derived", "garch_returns.csv"))

ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(data) <- ctry_codes

z = FALSE
# parameters
set.seed(123)
n_boot <- 200
thresholds <- c(0.75, 0.85, 0.90, 0.95, 0.99)
alpha = 0.05

dat <- as.matrix(data)
head(dat)

# estimate DAG on all observations
nodes <- ctry_codes
d <- length(nodes)

all_pruned_dags <- list()
all_bulk_dags <- list()

for (tau in thresholds) {
  cat("Threshold =", tau, "\n")
  pruned_dags_tau <- list()
  bulk_dags_tau <- list()
  
  # Bootstrap loop
  for (b in 1:n_boot) {
    
    # Bootstrap sample (with replacement)
    boot_idx <- sample(1:nrow(data), replace = TRUE)
    data_boot <- data[boot_idx, ]
    dat <- as.matrix(data_boot)

    binned_returns <- as.data.frame(dat)
    for (colname in names(binned_returns)) {
      binned_returns[[colname]] <- cut(
        binned_returns[[colname]],
        breaks = quantile(binned_returns[[colname]],
                          probs = c(0, 1/3, 2/3, 1), 
                          na.rm = TRUE),
        labels = c("Low", "Medium", "High"),
        include.lowest = TRUE
      )
    }
    
    blacklist <- data.frame(
      from = nodes[nodes != "EUR"],
      to = rep("EUR", length(nodes) - 1)
    )
    
    dag_binned <- hc(binned_returns, score = "bde", blacklist = blacklist)
    
    full_adjmat <- amat(dag_binned)
    
    bulk_dags_tau[[b]] <- full_adjmat
    
    # set up prunable pairs
    full_adjmat <- full_adjmat
    prunbl_prs <- prunable_pairs(full_adjmat) %>%
      arrange(desc(y)) %>%
      mutate(pair = map2(x, y, ~list(x = .x, y = .y))) %>%
      pull(pair)
    
    Y_mp <- data2mpareto(dat, p = tau)
    # estimate variogram and HR CI-test
    Gamma_hat <- emp_vario(Y_mp)
    ci_test_ext <- purrr::partial(
      HR_CItest,
      suffStat = list(Gamma = Gamma_hat, n = nrow(Y_mp)))
    
    # estimate extremal DAG via pruning
    B_pruned <- prune_edge_algo(
      B_full = full_adjmat,
      prunbl_prs = prunbl_prs,
      ci_test = ci_test_ext,
      alpha = alpha,
      base_pruning_algo = purrr::partial(prune_edge_fast, verbose = FALSE),
      separating_set_selector = mb_selector,
      shuffle = TRUE
    )
    
    colnames(B_pruned) <- ctry_codes
    rownames(B_pruned) <- ctry_codes
    pruned_dags_tau[[b]] <- B_pruned
    
    if (z == TRUE){ 
      cat("  Bootstrap", b, "of", n_boot, "tau = ", tau, "\n")
      
      g <- graph_from_adjacency_matrix(full_adjmat, mode = "directed")
      
      layout <- layout_with_sugiyama(g)$layout
      plot(
        g,
        layout = layout,
        vertex.label = V(g)$name,
        vertex.size = 30,
        vertex.label.cex = 1.2,
        vertex.label.color = "black",
        edge.arrow.size = 0.5,
        main = "Bulk DAG")
      
      g <- graph_from_adjacency_matrix(B_pruned, mode = "directed")
      
      layout <- layout_with_sugiyama(g)$layout
      plot(
        g,
        layout = layout,
        vertex.label = V(g)$name,
        vertex.size = 30,
        vertex.label.cex = 1.2,
        vertex.label.color = "black",
        edge.arrow.size = 0.5,
        main = "Extremal Pruned DAG"
      )}
    }
    all_pruned_dags[[as.character(tau)]] <- pruned_dags_tau
    all_bulk_dags[[as.character(tau)]] <- bulk_dags_tau
}


d <- length(ctry_codes)
edge_freq <- matrix(0, nrow = d, ncol = d)
colnames(edge_freq) <- rownames(edge_freq) <- ctry_codes

for (dag in all_pruned_dags[["0.95"]]) {
  edge_freq <- edge_freq + dag
}

edge_freq <- edge_freq / n_boot

edge_freq



final_dags <- list()

for (tau in thresholds) {
  
  d <- length(ctry_codes)
  edge_freq <- matrix(0, nrow = d, ncol = d)
  colnames(edge_freq) <- rownames(edge_freq) <- ctry_codes
  

  for (dag in all_pruned_dags[[as.character(tau)]]) {
    edge_freq <- edge_freq + dag
  }
  
  edge_freq <- edge_freq / n_boot  # Convert to frequency
  

  final_mat <- ifelse(edge_freq > 0.5, 1, 0)
  colnames(final_mat) <- rownames(final_mat) <- ctry_codes
  final_dags[[as.character(tau)]] <- final_mat
  

  g <- graph_from_adjacency_matrix(final_mat, mode = "directed", diag = FALSE)
  plot(
    g,
    main = paste("Majority DAG (tau =", as.character(tau), ")"),
    edge.arrow.size = 0.4,
    vertex.size = 30,
    vertex.label.cex = 0.8,
    layout = layout_with_sugiyama(g)$layout
  )
}

edge_freq <- matrix(0, nrow = d, ncol = d)
for (tau in thresholds) {
  bulk_list_tau <- all_bulk_dags[[as.character(tau)]]
  for (dag in bulk_list_tau) {
    edge_freq <- edge_freq + dag
    
  }
}

print(edge_freq)
edge_freq <- edge_freq / (n_boot * length(thresholds))

# Threshold to get final bulk DAG
final_bulk <- ifelse(edge_freq > 0.5, 1, 0)
colnames(final_bulk) <- rownames(final_bulk) <- ctry_codes

g <- graph_from_adjacency_matrix(final_bulk, mode = "directed", diag = FALSE)
plot(
  g,
  main = paste("Majority Bulk DAG"),
  edge.arrow.size = 0.4,
  vertex.size = 30,
  vertex.label.cex = 0.8,
  layout = layout_with_sugiyama(g)$layout
)
