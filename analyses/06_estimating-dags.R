# In this script we run our structure learning algorithms both in the bulk and 
# and in the extremes, using the pruning algorithm for the latter
# we also analyse the learned structure fit to our data

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
# load data
#===============================================================================

data <- read.csv(here("data", "derived", "garch-returns.csv"))
data <- subset(data, select = -Date)
ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(data) <- ctry_codes

#===============================================================================
# estimate bulk DAG
#===============================================================================

nodes <- colnames(data)

blacklist <- data.frame(
  from = nodes[nodes != "EUR"],
  to = rep("EUR", length(nodes) - 1))


whitelist <- data.frame(from = c("NOR","NOR", "CZE"), to = c("CZE","HUN", "SWE"))

set.seed(123)
full_dag <- tabu(data, blacklist = blacklist, whitelist = whitelist)

# plot full dag
graphviz.plot(full_dag)


#-------------------------------------------------------------------------------
# check that conditional independencies inferred from DAG match our data
#-------------------------------------------------------------------------------

# check conditional independence statements match data
arcs_df <- arcs(full_dag)
dagitty_str <- paste0(
  "dag {", 
  paste(apply(arcs_df, 1, function(x) paste0(x[1], " -> ", x[2])), 
        collapse = "\n"),
  "}"
)

# get implied independencies
implied_CIs <- impliedConditionalIndependencies(dagitty_str)
implied_CIs 
# test them against our data
ci_results <- localTests(
  x = dagitty_str,
  data = data,
  type = "cis"       
)
ci_results

#-------------------------------------------------------------------------------
# bootstrap for robustness and check fit of bootstrapped DAG
#-------------------------------------------------------------------------------

set.seed(123)
# bootstrap dags
boot_dag <- boot.strength(
  data = data,
  R = 200,  # number of bootstrap replications
  algorithm = "tabu",
  algorithm.args = list(blacklist = blacklist, whitelist = whitelist)
)

# keep edges with strength > 0.5
avg_dag <- averaged.network(boot_dag, threshold = 0.5)
avg_mat <- amat(avg_dag)

# plot average dag
avg_g <- graph_from_adjacency_matrix(avg_mat, mode = "directed")

layout <- layout_with_sugiyama(avg_g)$layout

plot(avg_g,
     layout = layout,
     vertex.label = V(avg_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     main = "Average Bulk DAG",
     vertex.color = "skyblue",
     edge.color = "black"
)

edge_direction_strength <- boot_dag[order(-boot_dag$strength), ]


arcs_df <- arcs(avg_dag)
dagitty_str <- paste0(
  "dag {", 
  paste(apply(arcs_df, 1, function(x) paste0(x[1], " -> ", x[2])), collapse = "\n"),
  "}"
)
dagitty_obj <- dagitty(dagitty_str)
# get implied independencies
implied_CIs <- impliedConditionalIndependencies(dagitty_str)
implied_CIs 
# test them against our data
ci_results <- localTests(
  x = dagitty_str,
  data = data,
  type = "cis"       
)
ci_results

#===============================================================================
# estimate pruned DAG
#===============================================================================

set.seed(123)
alpha = 0.05
full_adjmat <- amat(avg_dag)

# set up prunable pairs
prunbl_prs <- prunable_pairs(full_adjmat) %>%
  arrange(desc(y)) %>%
  mutate(pair = map2(x, y, ~list(x = .x, y = .y))) %>%
  pull(pair)

dat <- as.matrix(data)
head(dat)
tau <- 0.95
Y_mp <- data2mpareto(dat, p = tau)
cat("tau =", tau, ", n =", nrow(Y_mp), "\n")
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

extreme_g <- graph_from_adjacency_matrix(B_pruned, mode = "directed")

layout <- layout_with_sugiyama(extreme_g)$layout

plot(extreme_g,
     layout = layout,
     vertex.label = V(extreme_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     main = "Pruned DAG (Extreme Structure)",
     vertex.color = "skyblue",
     edge.color = "black"
)


# to assess model fit, we can't use same CI tests as in bulk DAG as 
# MPD not a product space so "standard" CI tests don't work

#===============================================================================
# test vs training set
#===============================================================================

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

prunbl_prs_train <- prunable_pairs(full_adjmat) %>% transpose()

B_pruned_train <- prune_edge_algo(
  B_full = full_adjmat,
  prunbl_prs = prunbl_prs_train,
  ci_test = ci_test_train,
  alpha = alpha,
  base_pruning_algo = prune_edge_fast,
  separating_set_selector = mb_selector
)

# test fit on test data
Gamma_test <- emp_vario(Y_test)
k_test <- nrow(Y_test)

colnames(B_pruned_train) <- NULL
rownames(B_pruned_train) <- NULL

ci_test_test <- purrr::partial(
  HR_CItest,
  suffStat = list(Gamma = Gamma_test, n = k_test)
)


res_test <- perform_citests(B_pruned_train, ci_test_fun = ci_test_test)

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




