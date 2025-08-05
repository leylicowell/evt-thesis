library(here)
library(dplyr)
library(igraph)
library(tidyr)
library(ggplot2)
library(igraph)
library(ggraph)
library(causalXtreme)
library(bnlearn)
library(pdfetch)
library(graphicalExtremes)
library(timeSeries)
library(fGarch)
library(tidyverse)
library(latex2exp)
library(rugarch)

source(here("analyses","load_packages.R"))

BoE_exchange_rates <- as.matrix(pdfetch_BOE(identifiers = c("XUDLERS", "XUDLDKS", "XUDLBK25","XUDLBK33", "XUDLBK47","XUDLSKS", "XUDLNKS"), 
                                            from =  "2005-10-01", to = "2020-09-30"))

save(BoE_exchange_rates, file = "~/Desktop/BoE_exchange_rates.RData")

load(file = "~/Desktop/BoE_exchange_rates.RData")
data <- BoE_exchange_rates
ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(BoE_exchange_rates) <- ctry_codes

# log returns
log_returns <-  apply(data, 2, FUN = function(x) diff(log(x))) #ndifflog

d <-  dim(log_returns)[2]
n <-  dim(log_returns)[1]

spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                   mean.model = list(armaOrder = c(0,2), include.mean = TRUE),
                   distribution.model = "norm")


fit <- lapply(as.data.frame(log_returns), function(x) {
  ugarchfit(spec = spec, data = x)
})

for (i in 1:ncol(log_returns)){
  print(fit[i])
}

garch_returns <- list()
for (name in names(fit)) {
  std_residuals <- residuals(fit[[name]], standardize = TRUE)
  garch_returns[[name]] <- as.numeric(std_residuals)
}
garch_returns <- as.data.frame(garch_returns)

write.csv(garch_returns, 
          file = file.path(here("data", "derived"), "garch_returns.csv"), 
          row.names = FALSE)

data <- garch_returns
colnames(data) <- ctry_codes

# parameters
set.seed(123)
thresholds<- c(0.95)
alpha = 0.05

dat <- as.matrix(data)
head(dat)

# estimate DAG on all observations
nodes <- ctry_codes

binned_returns <- as.data.frame(data)
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



# set up prunable pairs
full_adjmat <- full_adjmat
prunbl_prs <- prunable_pairs(full_adjmat) %>%
  arrange(desc(y)) %>%
  mutate(pair = map2(x, y, ~list(x = .x, y = .y))) %>%
  pull(pair)

    
tau <- 0.9
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
    main = "Pruned DAG (Extreme Structure)"
  )
  

