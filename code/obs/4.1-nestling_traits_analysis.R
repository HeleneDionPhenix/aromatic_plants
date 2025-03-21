# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nestling traits
# Nestling traits analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(lubridate)
library(vegan)
library(lmerTest)
library(effects)

# Import data -------------------------------------------------------------

rm(list = ls())

# nest composition
load("data/obs/0-data.RData")

# data on nestlings
chick <- read.csv2("data/obs/0-data_nestlings.csv",
                   row.names = 1)

set.seed(4)

# Format data -------------------------------------------------------------

# Keep only aromatic plant quantity post-hatching
compo.ele <- filter(compo.stage, stade == "elevage")

# merge tables
compo.ele <- left_join(compo.ele, chick, by = "clutch")
#date as day of the year
compo.ele$date_eclo <- yday(compo.ele$date_eclo)
compo.ele$date_mesure <- yday(compo.ele$date_mesure)

#calculate real age at measure
compo.ele$age_mesure <- compo.ele$date_mesure - compo.ele$date_eclo
hist(compo.ele$age_mesure)

#calculate difference between feather age and real age
#as indication of feather development
compo.ele$devel <- compo.ele$age_plume - compo.ele$age_mesure

#select variables
data.cond <- select(compo.ele,
                    chick_id,
                    poids, tarsed, docilite, devel,
                    ARO, pop, clutch)
summary(data.cond)
data.cond <- drop_na(data.cond)

#examine distribution
hist(data.cond$poids)
hist(data.cond$tarsed)
hist(data.cond$docilite)
hist(data.cond$devel)

#Scale ARO
hist(data.cond$ARO)
hist(log(data.cond$ARO + 1))

data.cond$ARO.log <- log(data.cond$ARO + 1)
data.cond$ARO.s <- scale(data.cond$ARO.log)[,1]

#Create nestling trait matrix (scaled)
cond <- select(data.cond,
               chick_id,
               poids, tarsed, docilite, devel)
rownames(cond) <- cond$chick_id
cond <- select(cond,
               -chick_id)
cond <- scale(cond)

#Create explanatory variable table
data.cond <- select(data.cond,
                    chick_id, pop, ARO.log, ARO.s, clutch)
rownames(data.cond) <- data.cond$chick_id
data.cond <- data.cond[rownames(cond),]
summary(data.cond)


# PCoA ------------------------------------------------------------------

pcoa <- capscale(cond~1, distance = "eucl")
pcoa
plot(pcoa$CA$eig/sum(pcoa$CA$eig))
plot(pcoa, type="text", scaling = 2)

# Model axis 1 ------------------------------------------------------------

#set reference level
data.cond$pop <- relevel(data.cond$pop, ref = "DMuro")

#extract nestling scores on first axis
data.cond$pcoa1 <- -scores(pcoa, 
                           choices = c(1), 
                           display = "sites",
                           scaling = 1)

mod.1 <- lmer(pcoa1 ~ pop*ARO.s + 
                (1|clutch),
              data = data.cond)
summary(mod.1)
plot(allEffects(mod.1))

#confidence intervals
confint.1 <- confint.merMod(mod.1)

#R2
rpt.1 <- rptGaussian(pcoa1 ~ pop*ARO.s + 
                       (1|clutch),
                     grname = c("clutch", "Fixed", "Residual"),
                     data = data.cond,
                     nboot = 100, npermut = 0, ratio = TRUE, adjusted = F)
rpt.1.adj <- rptGaussian(pcoa1 ~ pop*ARO.s + 
                       (1|clutch),
                     grname = c("clutch", "Residual"),
                     data = data.cond,
                     nboot = 100, npermut = 0, ratio = TRUE, adjusted = T)

# Model axis 2 ------------------------------------------------------------

data.cond$pcoa2 <- scores(pcoa, 
                          choices = c(2), 
                          display = "sites",
                          scaling = 1)

mod.2 <- lmer(pcoa2 ~ ARO.s*pop +
                (1|clutch),
              data = data.cond)
summary(mod.2)
plot(allEffects(mod.2))

#confidence intervals
confint.2 <- confint.merMod(mod.2)

#R2
rpt.2 <- rptGaussian(pcoa2 ~ pop*ARO.s + 
                       (1|clutch),
                     grname = c("clutch", "Fixed", "Residual"),
                     data = data.cond,
                     nboot = 100, npermut = 0, ratio = TRUE, adjusted = F)
rpt.2.adj <- rptGaussian(pcoa2 ~ pop*ARO.s + 
                       (1|clutch),
                     grname = c("clutch", "Residual"),
                     data = data.cond,
                     nboot = 100, npermut = 0, ratio = TRUE, adjusted = T)

# Clean and save environment --------------------------------------------------------

rm(compo.ele)

save.image("data/obs/4-output_models_nestling_traits.RData")
