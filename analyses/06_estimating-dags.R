# In this script we run our structure learning algorithms both in the bulk and 
# and in the extremes, using the pruning algorithm for the latter
# we also analyse the learned structure fit to our data
# for our pruning algorithm, we refer to https://github.com/nicolagnecco/extremeSCM

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
library(kableExtra)
library(knitr)

source(here("analyses","01_load_packages.R"))

#===============================================================================
# load data
#===============================================================================

data <- read.csv(here("data", "derived", "garch-returns.csv"))
data <- subset(data, select = -Date)
ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(data) <- ctry_codes

#===============================================================================
# plot dependencies in my dataset
#===============================================================================

set.seed(231)

# pick 4 random pairs from our dataset
n_pairs <- 4
rand_pairs <- replicate(n_pairs, {
  sample(colnames(data), 2, replace = FALSE)
}, simplify = FALSE)

plots_list <- list()
# plot
for (i in 1:n_pairs) {
  xvar <- rand_pairs[[i]][1]
  yvar <- rand_pairs[[i]][2]
  
  p <- ggplot(data, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.5, colour = "black", size = 1) +
    labs(
      title = paste0(xvar, " versus ", yvar),
      x = paste0(xvar, " Residuals"),
      y = paste0(yvar, " Residuals")
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 14)
    )
  
  plots_list[[i]] <- p
}


plot1 <- plots_list[[1]]
plot2 <- plots_list[[2]]
plot3 <- plots_list[[3]]
plot4 <- plots_list[[4]]

plot1
plot2
plot3
plot4

ggsave(here("outputs", "scatter1.pdf"), plot1, width = 6, height = 4)
ggsave(here("outputs", "scatter2.pdf"), plot2, width = 6, height = 4)
ggsave(here("outputs", "scatter3.pdf"), plot3, width = 6, height = 4)
ggsave(here("outputs", "scatter4.pdf"), plot4, width = 6, height = 4)

#===============================================================================
# estimate bulk DAG
#===============================================================================

nodes <- colnames(data)

blacklist <- data.frame(
  from = nodes[nodes != "EUR"],
  to = rep("EUR", length(nodes) - 1))

set.seed(123)
full_dag <- tabu(data, blacklist = blacklist)

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
  algorithm.args = list(blacklist = blacklist)
)

# keep edges with strength > 0.5
avg_dag <- averaged.network(boot_dag, threshold = 0.5)
avg_mat <- amat(avg_dag)

# plot average dag
avg_g <- graph_from_adjacency_matrix(avg_mat, mode = "directed")

layout <- layout_with_kk(avg_g)

pdf(here("outputs", "avg_bulk_dag.pdf"), width = 6, height = 6)
plot(avg_g,
     layout = layout,
     vertex.label = V(avg_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     main = "Bulk DAG (Bootstrapped)",
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
     edge.color = "black"
)
dev.off()

plot(avg_g,
     layout = layout,
     vertex.label = V(avg_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     main = "Bulk Network (Averaged across Bootstraps)",
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
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


ci_table <- ci_results %>%
  mutate(
    "p-value" = signif(p.value, 3),
    Decision = ifelse(p.value < 0.05, "Reject", "Do not reject")
  ) %>%
  dplyr::select("p-value", Decision)

kbl(
  ci_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "Tests of implied conditional independencies from the averaged DAG. 
  The null hypothesis is that the stated conditional independence holds in the data."
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "dag-ci-tests.tex"))

#-------------------------------------------------------------------------------
# assess model fit through conditional distributions
#-------------------------------------------------------------------------------

avg_dag_fit <- bn.fit(avg_dag,data)
bn.fit.histogram(avg_dag_fit)
bn.fit.qqplot(avg_dag_fit)
avg_dag_fit

#-------------------------------------------------------------------------------
# plot these (using my ggplot2 formatting)
#-------------------------------------------------------------------------------

resid_df <- data.frame()

for (var in names(avg_dag_fit)) {
  node <- avg_dag_fit[[var]]
  
  tmp <- data.frame(variable = var,
                    residuals = node$residuals,
                    fitted = node$fitted.values)
  resid_df <- rbind(resid_df, tmp)
}

norm_df <- data.frame()
for (var in unique(resid_df$variable)) {
  resid_data <- subset(resid_df, variable == var)
  mu <- mean(resid_data$residuals)
  sigma <- avg_dag_fit[[var]]$sd
  
  x_vals <- seq(min(resid_data$residuals), 
                max(resid_data$residuals), 
                length.out = 200)

  y_vals <- dnorm(x_vals, mean = mu, sd = sigma)
  
  tmp <- data.frame(
    variable = var,
    x = x_vals,
    y = y_vals
  )
  
  norm_df <- rbind(norm_df, tmp)
}

hist_plot <- ggplot(resid_df, aes(x = residuals)) +
  geom_histogram(bins = 50, aes(y = after_stat(density)), fill = "grey", colour = "black") +
  geom_line(data = norm_df, aes(x = x, y = y), colour = "#ED8221", linewidth = 1) +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Residuals", y = "Density") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 13.5),
    axis.title = element_text(size = 14)
  )

hist_plot

ggsave(here("outputs", "cond-dens-hist.pdf"), hist_plot, width = 10, height = 6)

qq_df <- data.frame()

for (var in unique(resid_df$variable)) {
  resid_data <- subset(resid_df, variable == var)
  n <- nrow(resid_data)
  
  sigma_hat <- avg_dag_fit[[var]]$sd
  
  emp_q <- sort(resid_data$residuals)
  theor_q <- qnorm(ppoints(n), 0, sigma_hat)
  
  tmp <- data.frame(
    variable = var,
    theor_q = theor_q,
    sample_q = emp_q
  )
  
  qq_df <- rbind(qq_df, tmp)
}

qq_plot <- ggplot(qq_df, aes(x = theor_q, y = sample_q)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#ED8221",
              linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    strip.text = element_text(face = "bold", size = 13.5),
    axis.title = element_text(size = 14)
  )

qq_plot

ggsave(here("outputs", "cond-dens-qq-plots.pdf"), qq_plot, width = 10, height = 6)

# plot table of edge frequencies

edge_df <- boot_dag[, c("from", "to", "strength", "direction")] %>%
  arrange(desc(strength), desc(direction)) %>%  
  mutate(strength  = round(strength, 3),
         direction = round(direction, 3))

colnames(edge_df) <- c("From", "To", "Frequency", "Direction Probability")

kbl(
  edge_df,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "Complete bootstrap edge frequencies and direction probabilities."
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "dag-edge-freqs-all.tex"))

# table of kept edges according to our 50% threshold
avg_arcs <- as.data.frame(arcs(avg_dag))

edge_df_main <- edge_df %>%
  inner_join(avg_arcs, by = c("From" = "from", "To" = "to"))

kbl(
  edge_df_main,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "Bootstrap edge frequencies and direction probabilities for edges retained in the averaged network (threshold = 0.5)."
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "dag-edge-freqs-main.tex"))

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
tau <- 0.925
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


pdf(here("outputs", "extreme_dag_925.pdf"), width = 6, height = 6)
layout <- layout_with_kk(extreme_g)
layout[which(V(extreme_g)$name == "CZE"), ] <- 
  layout[which(V(extreme_g)$name == "CZE"), ] + c(0, 0.4)
layout[which(V(extreme_g)$name == "HUN"), ] <- 
  layout[which(V(extreme_g)$name == "HUN"), ] + c(0, 0.2)
plot(extreme_g,
     layout = layout,
     vertex.label = V(extreme_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
     edge.color = "black"
)


dev.off()




# BOOTSTRAPPING PRUNING ALGORITHM
alpha <- 0.05
tau <- 0.925

set.seed(42)
n_boot <- 200  
boot_results <- list()

for (b in 1:n_boot) {
  idx <- sample(1:nrow(Y_mp), replace = TRUE)
  Y_boot <- Y_mp[idx, ]
  
  Gamma_boot <- emp_vario(Y_boot)
  ci_test_boot <- purrr::partial(HR_CItest,
                                 suffStat = list(Gamma = Gamma_boot, n = nrow(Y_boot)))
  
  # prune DAG
  B_pruned_boot <- prune_edge_algo(
    B_full = full_adjmat,
    prunbl_prs = prunbl_prs,
    ci_test = ci_test_boot,
    alpha = alpha,
    base_pruning_algo = purrr::partial(prune_edge_fast, verbose = FALSE),
    separating_set_selector = mb_selector,
    shuffle = TRUE
  )
  boot_results[[b]] <- B_pruned_boot
}

sum_mat <- matrix(
  0,
  nrow = nrow(boot_results[[1]]),
  ncol = ncol(boot_results[[1]]),
  dimnames = dimnames(boot_results[[1]])
)


for (b in seq_len(n_boot)) {
  sum_mat <- sum_mat + boot_results[[b]]
}


edge_freq <- sum_mat / n_boot
edge_freq

# keep edges that show up > 50% of time
B_majority <- ifelse(edge_freq >= 0.5, 1, 0)

B_majority

avg_extreme_g <- graph_from_adjacency_matrix(B_majority, mode = "directed")


pdf(here("outputs", "avg_extreme_dag_925.pdf"), width = 6, height = 6)
layout <- layout_with_kk(avg_extreme_g)
layout[which(V(avg_extreme_g)$name == "CZE"), ] <- 
  layout[which(V(avg_extreme_g)$name == "CZE"), ] + c(0, 0.4)
layout[which(V(avg_extreme_g)$name == "HUN"), ] <- 
  layout[which(V(avg_extreme_g)$name == "HUN"), ] + c(0, 0.2)
plot(avg_extreme_g,
     layout = layout,
     vertex.label = V(avg_extreme_g)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     main = expression(bold(paste("Extremal DAG (Bootstrapped), ", tau , " = 0.925"))),
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
     edge.color = "black"
)

dev.off()
# to assess model fit, we can't use same CI tests as in bulk DAG as 
# MPD not a product space so "standard" CI tests don't work

#===============================================================================
# test vs training set
#===============================================================================

set.seed(321)

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


#===============================================================================
# plots and tables for extremal DAG
#===============================================================================


edge_list <- list()

for (i in seq_len(nrow(edge_freq))) {
  for (j in seq_len(ncol(edge_freq))) {
    if (edge_freq[i, j] > 0) {
      edge_list[[length(edge_list) + 1]] <- list(
        From = rownames(edge_freq)[i],
        To   = colnames(edge_freq)[j],
        Frequency = edge_freq[i, j]
      )
    }
  }
}

edge_df_ext <- data.frame(
  From = sapply(edge_list, function(x) x$From),
  To = sapply(edge_list, function(x) x$To),
  Frequency = sapply(edge_list, function(x) x$Frequency)
) %>%
  arrange(desc(Frequency))

kbl(
  edge_df_ext,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "Bootstrap edge frequencies for the extremal DAG (threshold = 0.5)."
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "extremal-dag-edge-freqs.tex"))


cv_table <- res_summary %>%
  mutate(is_dsep = ifelse(is_dsep, "Implied independence", "Implied dependence"),
         prop_reject = round(prop_reject, 3)) %>%
  rename("DAG Implication" = is_dsep,
         "Number of Tests" = n_tests,
         "Proportion Rejected" = prop_reject)

kbl(
  cv_table,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = "c",
  caption = "Train/test validation of extremal DAG: proportion of conditional independence tests rejected at $\\alpha=0.05$."
) %>%
  kable_styling(latex_options = "hold_position", position = "center") %>%
  save_kable(here("outputs", "extremal-dag-cv.tex"))


set.seed(123)
B <- 500
boot_chi <- replicate(B, {
  idx <- sample(1:nrow(Y_mp), replace = TRUE)
  emp_chi(Y_mp[idx, ])
}, simplify = "array")

chi_low  <- apply(boot_chi, c(1,2), quantile, probs = 0.025)
chi_high <- apply(boot_chi, c(1,2), quantile, probs = 0.975)

get_upper <- function(mat) {
  idx <- which(upper.tri(mat), arr.ind = TRUE)
  data.frame(
    i = idx[,1],
    j = idx[,2],
    value = mat[idx]
  )
}

g <- graph_from_adjacency_matrix(B_majority, mode = "directed", diag = FALSE)

theor_chi <- Gamma2chi(complete_Gamma(Gamma=Gamma_hat, graph=g))
empirical_chi = emp_chi(Y_mp)

chi_df <- get_upper(theor_chi) %>%
  rename(chi_theor = value) %>%
  left_join(get_upper(empirical_chi) %>% rename(chi_emp = value),
            by = c("i","j")) %>%
  left_join(get_upper(chi_low) %>% rename(chi_low = value),
            by = c("i","j")) %>%
  left_join(get_upper(chi_high) %>% rename(chi_high = value),
            by = c("i","j"))


chi_plot <- ggplot(chi_df, aes(x = chi_theor, y = chi_emp)) +
  geom_errorbar(aes(ymin = chi_low, ymax = chi_high),
                width = 0, colour = "grey", alpha = 0.7) +
  geom_point(size = 2, alpha = 0.8, colour = "black") +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "#ED8221", linewidth = 0.9) +
  labs(x = TeX("DAG-implied $\\chi$"),
       y = TeX("Empirical $\\hat{\\chi}$"),
       title = "Pairwise Extremal Correlations") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14.5, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14)
  )

chi_plot

ggsave(here("outputs", "pairwise-chi-plots.pdf"), chi_plot, width = 6, height = 4)


#===============================================================================
# test algorithm for different taus
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
tau <- 0.9
Y_mp <- data2mpareto(dat, p = tau)
cat("tau =", tau, ", n =", nrow(Y_mp), "\n")
# estimate variogram and HR CI-test
Gamma_hat <- emp_vario(Y_mp)
ci_test_ext <- purrr::partial(
  HR_CItest,
  suffStat = list(Gamma = Gamma_hat, n = nrow(Y_mp)))

# estimate extremal DAG via pruning
B_pruned_90 <- prune_edge_algo(
  B_full = full_adjmat,
  prunbl_prs = prunbl_prs,
  ci_test = ci_test_ext,
  alpha = alpha,
  base_pruning_algo = purrr::partial(prune_edge_fast, verbose = FALSE),
  separating_set_selector = mb_selector,
  shuffle = TRUE
)

colnames(B_pruned_90) <- ctry_codes
rownames(B_pruned_90) <- ctry_codes

extreme_g_90 <- graph_from_adjacency_matrix(B_pruned_90, mode = "directed")

pdf(here("outputs", "extreme_dag_90.pdf"), width = 6, height = 6)
layout <- layout_with_kk(extreme_g_90)
plot(extreme_g_90,
     layout = layout,
     vertex.label = V(extreme_g_90)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
     edge.color = "black"
)


dev.off()




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
B_pruned_95 <- prune_edge_algo(
  B_full = full_adjmat,
  prunbl_prs = prunbl_prs,
  ci_test = ci_test_ext,
  alpha = alpha,
  base_pruning_algo = purrr::partial(prune_edge_fast, verbose = FALSE),
  separating_set_selector = mb_selector,
  shuffle = TRUE
)

colnames(B_pruned_95) <- ctry_codes
rownames(B_pruned_95) <- ctry_codes

extreme_g_95 <- graph_from_adjacency_matrix(B_pruned_95, mode = "directed")

pdf(here("outputs", "extreme_dag_95.pdf"), width = 6, height = 6)
layout <- layout_with_kk(extreme_g_95)
plot(extreme_g_95,
     layout = layout,
     vertex.label = V(extreme_g_95)$name,
     vertex.size = 30,
     vertex.label.cex = 1.2,
     vertex.label.color = "black",
     edge.arrow.size = 0.5,
     vertex.color = "#70B4EB",
     vertex.label.font = 2,
     edge.color = "black"
)


dev.off()

#===============================================================================
# check implied CIs of extremal average dag
#===============================================================================

B_majority_mat <- as.matrix(B_majority)
colnames(B_majority_mat) <- ctry_codes
rownames(B_majority_mat) <- ctry_codes

edges <- which(B_majority==1, arr.ind=TRUE)
arcs <- apply(edges, 1, function(x) {
  paste0(rownames(B_majority)[x[1]], " -> ", colnames(B_majority)[x[2]])
})

dag_str <- paste0("dag {", paste(arcs, collapse="\n"), "}")

dag_obj <- dagitty(dag_str)

# list implied conditional independencies of extremal bootstrapped DAG
impliedConditionalIndependencies(dag_obj)
