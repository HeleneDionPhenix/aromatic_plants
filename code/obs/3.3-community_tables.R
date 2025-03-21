# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Community tables
# Author: Hélène Dion-Phénix
# Last edition: 2024-03-03

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/obs/3-output_models_community.RData")

# Nest --------------------------------------------------------------------

tab <- as.data.frame(perm.nest.adj)

tab$Effect <- c("Aromatic plant quantity", 
                "Population",
                "Phenological stage",
                "Clutch",
                "Aromatic plant quantity:Population",
                "Aromatic plant quantity:Phenological stage",
                "Population:Phenological stage",
                "Aromatic plant quantity:Population:Phenological stage",
                "Residual",
                "Total")


tab <- select(tab,
              Effect, Df, SumOfSqs, F, parOmegaSq, 'Pr(>F)')


#round
tab$SumOfSqs <- format(round(as.numeric(tab$SumOfSqs), digits=2), nsmall = 2)
tab$F <- format(round(as.numeric(tab$F), digits=2), nsmall = 2)
tab$parOmegaSq <- format.pval(as.numeric(tab$parOmegaSq), digits = 1, eps = 0.001)
tab$`Pr(>F)` <- format.pval(as.numeric(tab$`Pr(>F)`), digits = 1, eps = 0.001)

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
  "Table 3" = tab.f,
  path = "table/Table_3.docx")

# Eggshell ----------------------------------------------------------------


tab <- as.data.frame(perm.eggshell.adj)

tab$Effect <- c("Aromatic plant quantity", 
                "Population",
                "Aromatic plant quantity:Population",
                "Residual",
                "Total")


tab <- select(tab,
              Effect, Df, SumOfSqs, F, parOmegaSq, 'Pr(>F)')


#round
tab$SumOfSqs <- format(round(as.numeric(tab$SumOfSqs), digits=2), nsmall = 2)
tab$F <- format(round(as.numeric(tab$F), digits=2), nsmall = 2)
tab$parOmegaSq <- format.pval(as.numeric(tab$parOmegaSq), digits = 1, eps = 0.001)
tab$`Pr(>F)` <- format.pval(as.numeric(tab$`Pr(>F)`), digits = 1, eps = 0.001)

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
  "Table 4" = tab.f,
  path = "table/Table_4.docx")
