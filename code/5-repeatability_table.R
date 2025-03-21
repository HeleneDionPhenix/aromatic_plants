# Description -------------------------------------------------------------
# Experimental & observational approach
# Create repeatability table for nest effect
# Sequencing: Illumina mi-seq June 2024
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-19

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/2-output_models_diversity.RData")
load("data/exp/3-output_models_community.RData")
load("data/obs/2-output_models_diversity.RData")
load("data/obs/3-output_models_community.RData")

rpt.extract <- function(rpt.object){
  R <- rpt.object$R[1]
  CI <- rpt.object$CI_emp[1,]
  CI <- paste("[",
              format(round(CI[1], digits=2), nsmall = 2),
              ", ",
              format(round(CI[2], digits=2), nsmall = 2),
              "]",
              sep = "")
  p <- rpt.object$P[1,1]
  return(c(R, CI, p))
}

tab <- data.frame(c(rpt.extract(rpt.trt.adj),rpt.extract(rpt.trt)))
names(tab) <- c("R2adj", "CIadj", "Padj", "R2", "CI", "P")


tab[2,] <- c("", "", "", perm.trt.adj$parOmegaSq[4], "", perm.trt.adj$`Pr(>F)`[4])
tab[3,] <- c(rpt.extract(rpt.shannon.nest.adj), rpt.extract(rpt.shannon.nest))
tab[4,] <- c("", "", "", perm.nest.adj$parOmegaSq[4], "", perm.nest.adj$`Pr(>F)`[4])
tab[5,] <- c(rpt.extract(rpt.aro.pop.adj), rpt.extract(rpt.aro.pop))

#round
tab$R2adj <- format(round(as.numeric(tab$R2adj), digits=2), nsmall = 2)
tab$R2 <- format(round(as.numeric(tab$R2), digits=2), nsmall = 2)
tab$P <- format.pval(as.numeric(tab$P), digits = 1,eps = 0.001)

tab$Model <- c("Experimental - Nest diversity",
               "Experimental - Nest composition",
               "Observational - Nest diversity",
               "Observational - Nest composition",
               "Observational - Aromatic plant cover")

tab <- select(tab,
              Model, R2adj, CIadj, R2, CI, P)

tabf <- flextable(tab)           


tabf <- tabf %>% 
  set_header_labels(Model = "Model",
                    R2adj = "R2",
                    CIadj = "95% CI",
                    R2 = "R2",
                    CI = "95% CI",
                    P = "p-value") %>%
  add_header_row(colwidths = c(1, 2, 2,1),
                 values = c("", "Adjusted","Unadjusted", "")) %>% 
  font(fontname='Times New Roman') %>%
  autofit()

tabf


save_as_docx(
  "Table 5" = tabf,
  path = "table/Table_5.docx")

