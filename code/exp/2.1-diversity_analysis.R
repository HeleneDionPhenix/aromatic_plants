# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Diversity analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-18

# Libraries --------------------------------------------------------------

library(dplyr)
library(vegan)
library(lmerTest)
library(effects)
library(rptR)
library(car)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/exp/0-data.RData")
load("data/exp/1-phyloseq_objects.RData")

rm(comm,comm.r, data, data.r, ps, ps.r, taxo, taxo.nest, taxo.r, taxo.plant,
   comm.plant, data.plant, ps.plant, ps.nest)

set.seed(1)

# Calculate Shannon index ----------------------------------------------------

# Verify concordance between rows
all(rownames(comm.nest)==rownames(data.nest))
data.nest <- data.nest[rownames(comm.nest),]

# Add shannon diversity to data.nest
data.nest$shannon <- vegan::diversity(comm.nest)

# Distribution
hist(data.nest$shannon, breaks = 20)

# Model 1 --------------------------------------------------------------------

# Set reference factors
data.nest$traitement <- factor(data.nest$traitement,
                               levels = c("ctr", 
                                          "ach", "pul", "men", "lav", "imm", 
                                          "lavmen", "lavimm", "menimm"))
data.nest$stage_pheno <- factor(data.nest$stage_pheno, 
                                levels = c("breed", "brood", "prebrood"))

mod.trt <- lmer(shannon ~  traitement*stage_pheno +
               (1|nichoir),
             data = data.nest)
summary(mod.trt)
plot(allEffects(mod.trt))
(confint.trt <- confint.merMod(mod.trt))

rpt.trt.adj <- rptGaussian(shannon ~  traitement *stage_pheno + 
              (1|nichoir), 
            grname = c("nichoir", "Residual"),
            data = data.nest,
            nboot = 100, npermut = 1000, ratio = TRUE, adjusted = T)

rpt.trt <- rptGaussian(shannon ~  traitement *stage_pheno + 
                         (1|nichoir), 
                       grname = c("Fixed", "nichoir", "Residual"),
                       data = data.nest,
                       nboot = 100, npermut = 1000, ratio = TRUE, adjusted = F)


# Levene test ------------------------------------------------------------

lev.stage <- leveneTest(shannon ~ stage_pheno,
                        data = data.nest)
lev.stage

# Save models output ------------------------------------------------------

save.image("data/exp/2-output_models_diversity.RData")
