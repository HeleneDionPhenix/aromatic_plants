# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Formatting data
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-05

# Libraries --------------------------------------------------------------

library(dplyr)
library(lubridate)

# Import data -------------------------------------------------------------

rm(list = ls())

# Information on collected nests
nest <- read.csv2("data/exp/0-data_nest.csv",
                  sep = ";",
                  dec = ",")

# Information on each samples
sample <- read.csv2("data/exp/0-data_sample.csv",
                    sep = ";",
                    dec = ",")

# Format data ---------------------------------------------------------------

# merge nest and sample info

data <- left_join(sample, nest, by = "nichoir")
rm(nest, sample)

# format variables
data$type <- as.factor(data$type)
data$nichoir <- as.factor(data$nichoir)
data$traitement <- as.factor(data$traitement)
data$category <- as.factor(data$category)
data$population <- as.factor(data$population)
data$station <- as.factor(data$station)
data$stage_pheno <- as.factor(data$stage_pheno)
data$stage_chick <- as.factor(data$stage_chick)
data$failure <- as.factor(data$failure)

#reorder factor treatments
data$traitement <- relevel(data$traitement, ref = "ctr")

#date
data$date_j1 <- ymd(data$date_j1)
data$date_pret <- ymd(data$date_pret)
data$date_recup <- ymd(data$date_recup)

# Number of days between nest ready and collection ------------------------

data$nday <- yday(data$date_recup) - yday(data$date_pret)
hist(data$nday)

rownames(data) <- data$sample.names

# Export data -------------------------------------------------------------

save.image("data/exp/0-data.RData")

