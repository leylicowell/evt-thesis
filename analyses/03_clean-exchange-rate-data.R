# In this script, we deal with missing values in our datasets and create a data 
# frame containing each time series

#===============================================================================
# load packages
#===============================================================================

library(here)
library(dplyr)
library(magrittr)
library(tidyr)

#===============================================================================
# load data
#===============================================================================

sar <- read.csv(here("data", "raw", "saudi-arabia.csv"))
egp <- read.csv(here("data", "raw", "egypt.csv"))
ils <- read.csv(here("data", "raw", "israel.csv"))
try <- read.csv(here("data", "raw", "turkey.csv"))
tnd <- read.csv(here("data", "raw", "tunisia.csv"))

#===============================================================================
# define function that filters each dataset
#===============================================================================

clean_data <- function(df, country_name) {
  # rename columns
  colnames(df) <- c("date", country_name)
  
  # convert date column to date type
  df$date <- as.Date(df$date, format = "%d/%m/%Y")
  
  
  # find all weekday dates for our data from 1995 to 2024
  all_dates <- seq(from = as.Date("2005-01-03"), 
                   to = as.Date("2024-12-31"), 
                   by = "day")
  weekday_dates <- all_dates[weekdays(all_dates) %in% c("Monday", 
                                                        "Tuesday", 
                                                        "Wednesday", 
                                                        "Thursday", 
                                                        "Friday")]
  
  # filter to the date range of interest
  df <- df %>%
    filter(date %in% weekday_dates) %>%
    arrange(date)
  
  # fill missing dates with NA
  all_dates <- data.frame(date = weekday_dates)
  
  df <- full_join(all_dates, df, by = "date") %>%
    arrange(date)
  
  
  return(df)
}

#===============================================================================
# clean datasets
#===============================================================================

sar_clean <- clean_data(sar, "SAR")
egp_clean <- clean_data(egp, "EGP")
ils_clean <- clean_data(ils, "ILS")
try_clean <- clean_data(try, "TRY")
tnd_clean <- clean_data(tnd, "TND")

#===============================================================================
# combine all datasets into one dataframe 
#===============================================================================

exchange_rates <- sar_clean %>%
  left_join(ils_clean, by = "date") %>%
  left_join(try_clean, by = "date") %>%
  left_join(tnd_clean, by = "date")


#===============================================================================
# replace NA values with last traded exchange rate
#===============================================================================

exchange_rates_filled <- exchange_rates %>%
  arrange(date) %>%
  fill(SAR, ILS, TRY, TND, .direction = "down")


#===============================================================================
# we choose to save our new data set as a .csv to facilitate readability and to 
# facilitate using this data set on different platforms/with different software
#===============================================================================

write.csv(exchange_rates_filled, 
          file = file.path(here("data", "derived"), "currency-data.csv"), 
          row.names = FALSE)


