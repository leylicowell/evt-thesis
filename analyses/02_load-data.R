# In this script, we load our data sets and create data 
# frames containing each time series, their respective log returns and
# their garch standardised residuals

#===============================================================================
# load packages
#===============================================================================

library(here)
library(dplyr)
library(magrittr)
library(tidyr)
library(pdfetch)
library(rugarch)

#===============================================================================
# load data
#===============================================================================

data <- as.matrix(pdfetch_BOE(identifiers = c("XUDLERS", 
                                              "XUDLDKS", 
                                              "XUDLBK25",
                                              "XUDLBK33", 
                                              "XUDLBK47",
                                              "XUDLSKS", 
                                              "XUDLNKS"), 
                              from =  "2005-10-01", 
                              to = "2020-09-30"))

ctry_codes <- c("EUR","DNK", "CZE","HUN","POL","SWE", "NOR")
colnames(data) <- ctry_codes

#===============================================================================
# calculate log returns and extract garch residuals
#===============================================================================

# log returns
log_returns <-  apply(data, 2, FUN = function(x) diff(log(x))) 

d <-  dim(log_returns)[2]
n <-  dim(log_returns)[1]

# garch residuals
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

head(garch_returns)
#===============================================================================
# we choose to save our new data sets as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#===============================================================================

data <- data.frame(Date = as.Date(rownames(data)), data)
log_returns <- data.frame(Date = as.Date(rownames(data)[2:nrow(data)]), log_returns)
garch_returns <- data.frame(Date = as.Date(rownames(data)[2:nrow(data)]), garch_returns)

write.csv(data, 
          file = file.path(here("data", "derived"), "currency-data.csv"), 
          row.names = FALSE)

write.csv(log_returns, 
          file = file.path(here("data", "derived"), "log-returns.csv"), 
          row.names = FALSE)

write.csv(garch_returns, 
          file = file.path(here("data", "derived"), "garch-returns.csv"), 
          row.names = FALSE)


