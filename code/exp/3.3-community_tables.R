# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Community analysis - Table creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-19

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/exp/3-output_models_community.RData")


# Nest samples ------------------------------------------------------------


tab <- as.data.frame(perm.trt.adj)

tab$Effect <- c("Treatment vs control", 
                "Phenological stage",
                "Plant species",
                "Nest",
                "Treatment vs control:Phenological stage",
                "Plant species:Phenological stage",
                "Residual",
                "Total")


tab <- select(tab,
              Effect, Df, SumOfSqs, F, parOmegaSq, 'Pr(>F)')

#reorder rows
tab <- tab[c(1,3,2,4:8),]

#round
tab$SumOfSqs <- format(round(as.numeric(tab$SumOfSqs), digits=2), nsmall = 2)
tab$F <- format(round(as.numeric(tab$F), digits=2), nsmall = 2)
tab$parOmegaSq <- format.pval(as.numeric(tab$parOmegaSq), digits = 1,eps = 0.001)
tab$`Pr(>F)` <- format.pval(as.numeric(tab$`Pr(>F)`), digits = 1,eps = 0.001)

tab.f <- flextable(tab) %>%
  set_header_labels(Effect = "Effect",
                    Df = "Df",
                    SumOfSqs = "Sum of squares",
                    "F" = "Pseudo-F",
                    parOmegaSq = "Adjusted-R2",
                    "Pr(>F)" = "P-value") %>%
  font(fontname='Times New Roman') %>%
  autofit()

tab.f

save_as_docx(
  "Table 2" = tab.f,
  path = "table/Table_2.docx")


# Plant sample ------------------------------------------------------------


tab <- as.data.frame(perm.pl.adj)

tab$Effect <- c("Plant species",
                "Residual",
                "Total")


tab <- select(tab,
              Effect, Df, SumOfSqs, F, parOmegaSq, 'Pr(>F)')

#round
tab$SumOfSqs <- format(round(as.numeric(tab$SumOfSqs), digits=2), nsmall = 2)
tab$F <- format(round(as.numeric(tab$F), digits=2), nsmall = 2)
tab$parOmegaSq <- format.pval(as.numeric(tab$parOmegaSq), digits = 2,eps = 0.001)
tab$`Pr(>F)` <- format.pval(as.numeric(tab$`Pr(>F)`), digits = 2,eps = 0.001)

tab.f <- flextable(tab) %>%
  set_header_labels(Effect = "Effect",
                    Df = "Df",
                    SumOfSqs = "Sum of squares",
                    "F" = "Pseudo-F",
                    parOmegaSq = "Adjusted-R2",
                    "Pr(>F)" = "P-value") %>%
  font(fontname='Times New Roman') %>%
  autofit()

tab.f

save_as_docx(
  "Table S2.1" = tab.f,
  path = "table/Table_S2.1.docx")
