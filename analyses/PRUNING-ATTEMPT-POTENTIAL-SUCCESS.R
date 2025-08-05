library(pdfetch)
BoE_exchange_rates <- as.matrix(pdfetch_BOE(identifiers = c("XUDLERS", "XUDLDKS", "XUDLBK25","XUDLBK33", "XUDLBK47","XUDLSKS", "XUDLNKS"), 
                                            from =  "2005-10-01", to = "2020-09-30"))


save(BoE_exchange_rates, file = "~/Desktop/BoE_exchange_rates.RData")

load(file = "~/Desktop/BoE_exchange_rates.RData")
data <- BoE_exchange_rates
ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(BoE_exchange_rates) <- ctry_codes
library("igraph")
library("graphicalExtremes")
library("timeSeries")
library("fGarch")
library("tidyverse")
library("latex2exp")

# create negative log ratios
nlr <-  apply(data, 2, FUN = function(x) diff(log(x))) #ndifflog

d <-  dim(nlr)[2]
n <-  dim(nlr)[1]

x <- nlr
fit.fGarch <- list()
residus <-  matrix(NA,n,d)
colnames(residus) <-  colnames(nlr)

for(i in 1:d){
  ga <- NULL
  form <- paste("~arma(0,2)+garch(1,1)")
  ga <- garchFit(formula=as.formula(form), data=x[,i], trace=FALSE, cond.dist="norm")
  residus[,i] <- ga@residuals/ga@sigma.t
  fit.fGarch[[i]] <- ga
  cat("Done for ",i,"/",d,"\n")
}


colnames(residus) <-  colnames(nlr)
abs.residus <- abs(residus)
acf(abs.residus[,7])

data <- residus
graph.full <- make_full_graph(d)






library(bnlearn)

colnames(data) <- ctry_codes
nodes <- c("EUR", "DNK", "CZE", "HUN", "POL", "SWE", "NOR")

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

graphviz.plot(dag_binned, main = "binned data")

head(data)

dag_init <- hc(as.data.frame(data),blacklist = blacklist)

graphviz.plot(dag_init, main = "continuous data")
