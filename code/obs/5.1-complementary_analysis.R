# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nestling traits
# Nestling traits analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-19

# Librairies --------------------------------------------------------------

library(dplyr)
library(vegan)
library(lmerTest)
library(effects)
library(rptR)
library(lubridate)

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

data.cond <- compo.ele
rm(compo.ele)

#Scale ARO
hist(data.cond$ARO)
hist(log(data.cond$ARO + 1))

data.cond$ARO.log <- log(data.cond$ARO + 1)
data.cond$ARO.s <- scale(data.cond$ARO.log)[,1]

# scale date by pop
scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

data.cond <- 
  data.cond %>%
  group_by(pop) %>%
  mutate(date_eclo.s = scale_this(date_eclo))
rm(scale_this)

# Model S5.1 - Developmental index ~ hatching date ----------------------------

hist(data.cond$devel)

mod.devel <- lmer(devel ~ date_eclo.s+pop +
                   (1|clutch),
                 data = data.cond)
summary(mod.devel)
plot(allEffects(mod.devel))


# Confidence intervals

confint.devel <- confint.merMod(mod.devel)

# R2
rpt.devel <- rptGaussian(devel ~ date_eclo.s+pop +
                           (1|clutch),
                           grname = c("clutch", "Fixed", "Residual"),
                           data = data.cond,
                           nboot = 100, npermut = 0, ratio = T, adjusted = F)
rpt.devel.adj <- rptGaussian(devel ~ date_eclo.s+pop +
                           (1|clutch),
                         grname = c("clutch", "Fixed", "Residual"),
                         data = data.cond,
                         nboot = 100, npermut = 0, ratio = T, adjusted = T)

# Model S5.2 - Aromatic plant quantity ~ hatching date ------------------------

data.nest <- select(data.cond,
                    ARO.log, ARO.s, date_eclo.s, date_eclo, pop)
data.nest <- distinct(data.nest)

mod.aro <- lm(ARO.log ~ date_eclo.s + pop,
              data = data.nest)

summary(mod.aro)
plot(allEffects(mod.aro))


# Confidence intervals

confint.aro <- confint.lm(mod.aro)


# Model S5.3 - Fledging success over the years ------------------------------

## Import data -------------------------------------------------------------

data.fledg <- read.csv2("data/obs/0-data_fledging.csv")


## Format data -------------------------------------------------------------

data.fledg$pop <- as.factor(data.fledg$pop)
data.fledg$an <- as.factor(data.fledg$an)

mod.fledg <- glm(cbind(pulenv,echec_fledg) ~ pop*an,
                 data = data.fledg,
                 family = "binomial")
summary(mod.fledg)
plot(allEffects(mod.fledg))

# Model S5.4 - Mean fledging success in the last 20 years ---------------------

mod.fledg.m <- glm(cbind(pulenv,echec_fledg) ~ pop,
                   data = data.fledg,
                   family = "binomial")
summary(mod.fledg.m)
plot(allEffects(mod.fledg.m))

# mean over 20 years
newdat.m <- expand.grid(pop = levels(data.fledg$pop))

pred <- predict(mod.fledg.m,
                newdat.m,
                se.fit = TRUE,
                type = "link")

inv_link <- family(mod.fledg.m)$linkinv

newdat.m$fledg <- inv_link(pred$fit)
newdat.m$lowerCI <- inv_link(pred$fit - 1.96*pred$se.fit)
newdat.m$upperCI <- inv_link(pred$fit + 1.96*pred$se.fit)

mean.fledg.tab <- newdat.m
rm(newdat.m)

# Save environment --------------------------------------------------------

save.image("data/obs/5-output_models_complementary.RData")
