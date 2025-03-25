# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Community analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(vegan)
library(MicEco)

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

set.seed(10)

# Format plant cover variables --------------------------------------------

hist(compo.stage$ARO)
hist(log(compo.stage$ARO + 1))
compo.stage$ARO.log <- log(compo.stage$ARO + 1)
compo.stage$ARO.s <- scale(compo.stage$ARO.log)[,1]

# Nest --------------------------------------------------------------------

# Hellinger distance matrix
compo.stage <- compo.stage[rownames(comm.nest.r),]
dist.comm.nest <- vegdist(comm.nest.r, method = "hellinger")

## Tests of multivariate dispersion ----------------------------------------

# between phenological stages
disper.nest.stage <- betadisper(dist.comm.nest, compo.stage$stade)
anova(disper.nest.stage)
plot(disper.nest.stage)

# between populations
disper.nest.pop <- betadisper(dist.comm.nest, compo.stage$pop)
anova(disper.nest.pop)
plot(disper.nest.pop)

## Model 6 - dbRDA ------------------------------------------------------------

mod.nest <- capscale(comm.nest.r ~ ARO.s*pop*stade +
                      clutch,
                     data = compo.stage,
                     distance = "hellinger")

perm.nest <- adonis2(dist.comm.nest ~ ARO.s*pop*stade +
                      clutch,
                      data = compo.stage,
                    by = "term")
perm.nest.adj <- adonis_OmegaSq(perm.nest, partial = TRUE)

mod.nest
p.nest <- anova(mod.nest)
p.nest <- p.nest$`Pr(>F)`[1]
anova(mod.nest, by = "term")
RsquareAdj(mod.nest)

perm.nest.adj

# Eggshell ----------------------------------------------------------------

# Hellinger distance matrix
compo.couv <- compo.stage[rownames(comm.eggshell.r),]
dist.comm.eggshell <- vegdist(comm.eggshell.r, method = "hellinger")

## Tests of multivariate dispersion ----------------------------------------

# between populations
disper.eggshell.pop <- betadisper(dist.comm.eggshell, compo.couv$pop)
anova(disper.eggshell.pop)
plot(disper.eggshell.pop)

## Model 7 - dbRDA ------------------------------------------------------------

mod.eggshell <- capscale(comm.eggshell.r ~ ARO.s*pop,
                    data = compo.couv,
                    distance = "hellinger")

perm.eggshell <- adonis2(dist.comm.eggshell ~ ARO.s*pop,
                    data = compo.couv,
                    by = "term")
perm.eggshell.adj <- adonis_OmegaSq(perm.eggshell, partial = TRUE)

mod.eggshell
p.eggshell <- anova(mod.eggshell)
p.eggshell <- p.eggshell$`Pr(>F)`[1]
anova(mod.eggshell, by = "term")
RsquareAdj(mod.eggshell)

perm.eggshell
perm.eggshell.adj
plot(mod.eggshell)

# Clean and save environment ----------------------------------------------

rm(comm.eggshell.r, comm.nest.r, 
   disper.eggshell.pop, disper.nest.pop, disper.nest.stage,
   perm.eggshell, perm.nest, dist.comm.eggshell, dist.comm.nest)

save.image("data/obs/3-output_models_community.RData")      
