# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Diversity analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03


# Librairies --------------------------------------------------------------

library(dplyr)
library(vegan)
library(lmerTest)
library(effects)
library(rptR)

# Import data -------------------------------------------------------------

rm(list = ls())

# nest composition
load("data/obs/0-data.RData")

# nest bacterial microbiota
load("data/obs/1-phyloseq_objects_nest.RData")
rm(comm.nest, mdata.nest, mdata.nest.r, 
   ps.nest, ps.nest.r, taxo.nest, taxo.nest.r)

# eggshell bacterial microbiota
load("data/obs/1-phyloseq_objects_eggshell.RData")
rm(comm.eggshell, mdata.eggshell, mdata.eggshell.r, 
   ps.eggshell, ps.eggshell.r, taxo.eggshell, taxo.eggshell.r)

set.seed(8)

# Nest bacterial diversity ------------------------------------------------

## Calcul of Shannon diversity ---------------------------------------------

all(rownames(comm.nest.r) %in% rownames(compo.stage))
compo.stage <- compo.stage[rownames(comm.nest.r),]

compo.stage$shannon <- vegan::diversity(comm.nest.r)
hist(compo.stage$shannon, breaks = 20)

sum(is.na(compo.stage))
hist(compo.stage$ARO)
hist(log(compo.stage$ARO + 1))
compo.stage$ARO.log <- log(compo.stage$ARO + 1)
compo.stage$ARO.s <- scale(compo.stage$ARO.log)[,1]

## Models ------------------------------------------------------------------

# Set reference levels
compo.stage$stade <- relevel(compo.stage$stade, ref = "elevage")
compo.stage$pop <- relevel(compo.stage$pop, ref = "EPirio")

# model
mod.shannon.nest <- lmer(shannon ~ ARO.s*pop*stade +
                  (1|clutch),
                data = compo.stage)

summary(mod.shannon.nest)
plot(allEffects(mod.shannon.nest))


## Confidence intervals and repeatabilities --------------------------------

confint.nest <- confint.merMod(mod.shannon.nest)

rpt.shannon.nest.adj <- rptGaussian(shannon ~ ARO.s*pop*stade +
              (1|clutch), 
            grname = c("clutch", "Residual"),
            data = compo.stage,
            nboot = 100, npermut = 0, ratio = TRUE, adjusted = T)
rpt.shannon.nest <- rptGaussian(shannon ~ ARO.s*pop*stade +
                                  (1|clutch), 
                                grname = c("clutch", "Fixed", "Residual"),
                                data = compo.stage,
                                nboot = 100, npermut = 0, ratio = TRUE, adjusted = F)


# Eggshell bacterial diversity --------------------------------------------

## Calcul of Shannon diversity ---------------------------------------------

all(rownames(comm.eggshell.r) %in% compo.stage$id)
compo.incu <- compo.stage[rownames(comm.eggshell.r),]

compo.incu$shannon <- vegan::diversity(comm.eggshell.r)
hist(compo.incu$shannon, breaks = 20)

sum(is.na(compo.incu))
hist(compo.incu$ARO)
hist(log(compo.incu$ARO + 1))
compo.incu$ARO.log <- log(compo.incu$ARO + 1)
compo.incu$ARO.s <- scale(compo.incu$ARO.log)[,1]

## Models ------------------------------------------------------------------

# set reference levels
compo.incu$pop <- relevel(compo.incu$pop, ref = "EPirio")

mod.shannon.eggshell <- lm(shannon ~ ARO.s*pop,
                           data = compo.incu)

summary(mod.shannon.eggshell)
plot(allEffects(mod.shannon.eggshell))


## Confidence intervals ----------------------------------------------------

confint.eggshell <- confint(mod.shannon.eggshell)


# Aromatic plant quantities between populations ------------------------------

hist(compo$ARO)
hist(log(compo$ARO + 1))
compo$ARO.log <- log(compo$ARO + 1)


## Model -------------------------------------------------------------------

# set reference levels
compo$pop <- relevel(compo$pop, ref = "DMuro")
compo$stade <- relevel(compo$stade, ref = "elevage")

mod.aro.pop <- lmer(ARO.log ~ pop*stade +
                  (1|clutch),
                data = compo)

summary(mod.aro.pop)
plot(allEffects(mod.aro.pop))


# Confidence intervals and repeatabilities --------------------------------

confint.aro.pop <- confint.merMod(mod.aro.pop)

rpt.aro.pop.adj <- rptGaussian(ARO.log ~ pop*stade +
              (1|clutch),
            grname = c("clutch", "Residual"),
            data = compo,
            nboot = 100, npermut = 0, ratio = T, adjusted = T)
rpt.aro.pop <- rptGaussian(ARO.log ~ pop*stade +
                             (1|clutch),
                           grname = c("clutch", "Fixed", "Residual"),
                           data = compo,
                           nboot = 100, npermut = 0, ratio = T, adjusted = F)


# clean and save environment --------------------------------------------------------

rm(comm.nest.r, comm.eggshell.r)

save.image("data/obs/2-output_models_diversity.RData")
