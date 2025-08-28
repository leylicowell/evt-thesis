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
library(mvtnorm)

source(here("analyses","01_load_packages.R"))

#===============================================================================
# simulate data from a bulk DAG SCM
#===============================================================================



#-------------------------------------------------------------------------------
# Simulate data from a bulk DAG with weighted edges
#-------------------------------------------------------------------------------
sample_bulk_SCM <- function(n, B_adj, W, Sigma = diag(nrow(B_adj)),
                            nonlinear = FALSE,
                            noise = c("gaussian", "student")) {
  
  d <- nrow(B_adj)
  noise <- match.arg(noise)
  
  # Causal ordering from binary adjacency
  caus_order <- causalXtreme:::compute_caus_order(B_adj) 
  
  # Noise generation
  if (noise == "gaussian") {
    eps <- mvtnorm::rmvnorm(n, sigma = Sigma)
  } else if (noise == "student") {
    eps <- mvtnorm::rmvt(n, sigma = Sigma, df = 7)
  }
  
  # Data matrix
  X <- matrix(0, nrow = n, ncol = d)
  
  # Root node
  X[, caus_order[1]] <- eps[, caus_order[1]]
  
  # Other nodes
  for (j in caus_order[-1]) {
    if (nonlinear) {
      X_transformed <- 1 / (1 + X^2)
      X[, j] <- (X %*% W[, j, drop = FALSE]) +
        (X_transformed %*% W[, j, drop = FALSE]) +
        eps[, j]
    } else {
      X[, j] <- (X %*% W[, j, drop = FALSE]) + eps[, j]
    }
  }
  
  colnames(X) <- paste0("X", 1:d)
  return(X)
}

#------------------------------------------------------
# Example usage
#------------------------------------------------------

# Binary adjacency
B_adj <- matrix(0, 4, 4)
B_adj[1,2] <- 1  # A->B
B_adj[1,3] <- 1  # A->C
B_adj[2,4] <- 1  # B->D
B_adj[3,4] <- 1  # C->D

# Weights (can be >1 to strengthen edge
W <- matrix(0, 4, 4)
W[1,2] <- 1   # A->B
W[1,3] <- 3   # A->C
W[2,4] <- 1   # B->D
W[3,4] <- 3   # C->D

# Gaussian linear case (correctly specified for Tabu+BIC)


# Heavy-tailed nonlinear case (misspecified, like GARCH residuals)
X_student <- sample_bulk_SCM(n = 2000, B = B, nonlinear = TRUE, noise = "student")



#------------------------------------------------------
# learn DAG with Tabu algorithm
#------------------------------------------------------

# Learn from Gaussian data
X_gauss <- sample_bulk_SCM(n = 3000, B_adj = B_adj, W= W, nonlinear = FALSE, noise = "gaussian")
bn_gauss <- tabu(as.data.frame(X_gauss))
plot(bn_gauss)   # quick DAG plot

# Learn from Student-t data
X_student <- sample_bulk_SCM(n = 10000, B = B_adj, W= W, nonlinear = FALSE, noise = "student")
bn_student <- tabu(as.data.frame(X_student))
plot(bn_student)

#------------------------------------------------------
# Compare learned vs true DAG
#------------------------------------------------------

# True DAG as bn object
true_bn <- empty.graph(nodes = colnames(X_gauss))
amat(true_bn) <- B   # insert adjacency matrix

# Structural Hamming Distance (SHD)
shd_gauss   <- shd(amat(bn_gauss), true_bn)
shd_student <- shd(bn_student, true_bn)

cat("SHD (Gaussian, linear):", shd_gauss, "\n")
cat("SHD (Student-t, nonlinear):", shd_student, "\n")


