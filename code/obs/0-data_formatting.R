# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Format nest composition data
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-20

# Librairies --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)

# Import data -------------------------------------------------------------

rm(list = ls())

# nest composition observations in 2023
compo <- read.csv2("data/obs/0-nest_composition.csv",
                   row.names = 1,
                   sep = ";",
                   dec = ",")


# Summarize composition  --------------------------------------------------

# Change codes for medians of the categories
# Classes (code : description : median)
#     0.5 : 1-2 fragments, < 5%  : 1%
#       1 : 3-5 fragments, < 5%  : 3%
#       2 : 5%-25%               : 15%
#       3 : 25%-50%              : 37.5%
#       4 : 50%-75%              : 62.5%
#       5 : 75%-100%             : 87.5%

for (i in 9:13){
  compo[which(compo[,i] == 5), i] <- 87.5
  compo[which(compo[,i] == 4), i] <- 62.5
  compo[which(compo[,i] == 3), i] <- 37.5
  compo[which(compo[,i] == 2), i] <- 15
  compo[which(compo[,i] == 1), i] <- 3
  compo[which(compo[,i] == 0.5), i] <- 1
}

# Sum the five focal species of aromatic plants
compo$ARO <- NA 
compo$ARO <- with(compo, ACH + MEN + IMM + LAV + PUL)

# remove species
compo <- select(compo,
                -ACH, -MEN, -IMM, -LAV, -PUL)

# Format variables 
compo$station <- as.factor(compo$station)
compo$pop <- as.factor(compo$pop)
compo$date <- ymd(compo$date)
compo$stade <- as.factor(compo$stade)

# Create composition data per clutch per stage ----------------------------

compo.stage <- select(compo,
                      -date,
                      -obs)

compo.stage.ARO <- compo.stage %>%
  group_by(id) %>%
  summarize(ARO = mean(ARO))

compo.stage.ARO <- ungroup(compo.stage.ARO)

compo.stage <- select(compo.stage, -ARO)
compo.stage <- left_join(compo.stage.ARO, compo.stage, by = "id")
compo.stage <- distinct(compo.stage)

compo.stage <- select(compo.stage,
                      id, clutch, nichoir, station, pop, stade, ARO)
compo.stage <- as.data.frame(compo.stage)

rownames(compo.stage) <- compo.stage$id

rm(compo.stage.ARO, i)


# clean and save ----------------------------------------------------------

save.image("data/obs/0-data.RData")
