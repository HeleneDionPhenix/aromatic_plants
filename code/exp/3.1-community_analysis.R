# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Community analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-18

# Libraries --------------------------------------------------------------

library(dplyr)
library(vegan)
library(MicEco)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/1-phyloseq_objects.RData")

rm(comm,comm.r, data, data.r, ps, ps.r, taxo, taxo.nest, taxo.r, taxo.plant)

set.seed(2)

# Plant bacterial microbiota ----------------------------------------------

# Hellinger distance matrix
dist.plant <- vegdist(comm.plant, method = "hellinger")

## Tests of multivariate dispersion ----------------------------------------

disper.plant <- betadisper(dist.plant, data.plant$traitement)
anova(disper.plant)
plot(disper.plant)

## Model S3.1 - dbRDA -------------------------------------------------------

mod.pl <- capscale(comm.plant ~ traitement,
                   data = data.plant,
                   distance = "hellinger")

perm.pl <- adonis2(dist.plant ~ traitement,
                   data = data.plant,
                   by = "term")
perm.pl.adj <- adonis_OmegaSq(perm.pl)

mod.pl
p.value.pl <- anova(mod.pl)
p.value.pl <- p.value.pl$`Pr(>F)`[1]
perm.pl.adj
RsquareAdj(mod.pl)

perm.pl.adj

# Nest bacterial microbiota -----------------------------------------------

# Hellinger distance matrix
dist.nest <- vegdist(comm.nest, method = "hellinger")

## Tests of multivariate dispersion ----------------------------------------

# Nest box
table(data.nest$nichoir)
disper.nest <- betadisper(dist.nest, data.nest$nichoir)
anova(disper.nest)
plot(disper.nest, label = F)

# Treatment
table(data.nest$traitement)
disper.trt <- betadisper(dist.nest, data.nest$traitement)
anova(disper.trt)
plot(disper.trt, label = T)

#Treatment vs controls
data.nest$trtctrl <- as.character(data.nest$traitement)
data.nest[which(data.nest$traitement != "ctr"), "trtctrl"] <- "trt"
data.nest$trtctrl <- as.factor(data.nest$trtctrl)
table(data.nest$trtctrl)
disper.trt <- betadisper(dist.nest, data.nest$trtctrl)
anova(disper.trt)
plot(disper.trt, label = T)

# Stage pheno
table(data.nest$stage_pheno)
disper.stage <- betadisper(dist.nest, data.nest$stage_pheno)
anova(disper.stage)
plot(disper.stage, label = T)

# Model 2 - dbRDA -------------------------------------------------------------

mod.trt <- capscale(comm.nest ~ trtctrl*stage_pheno + 
                      traitement*stage_pheno + 
                      nichoir,
                    data = data.nest,
                    distance = "hellinger")

perm.trt <- adonis2(dist.nest ~ trtctrl*stage_pheno + 
                      traitement*stage_pheno +
                      nichoir,
                    data = data.nest,
                    by = "term")
perm.trt.adj <- adonis_OmegaSq(perm.trt, partial = TRUE)

mod.trt
p.value.trt <- anova(mod.trt)
p.value.trt <- p.value.trt$`Pr(>F)`[1]
RsquareAdj(mod.trt)
perm.trt.adj

# Save environment --------------------------------------------------------

save.image("data/exp/3-output_models_community.RData")
