# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Diversity analysis - Figure creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-18


# Libraries --------------------------------------------------------------

library(flextable)
library(dplyr)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/2-output_models_diversity.RData")


# Table -------------------------------------------------------------------

summary.trt <- summary(mod.trt)

tab.trt <- data.frame(summary.trt$coefficients)

tab.trt <- select(tab.trt,
                  -df)

#add CI
confint.trt <- data.frame(confint.trt)
names(confint.trt) <- c("L", "H")
confint.trt$CI <- paste("[",
                        format(round(confint.trt$L, digits=2), nsmall = 2),
                        ", ",
                        format(round(confint.trt$H, digits=2), nsmall = 2),
                        "]",
                        sep = "")
confint.trt <- confint.trt[rownames(tab.trt),]
tab.trt$CI <- confint.trt$CI

#rename
names(tab.trt) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.trt$Estimate <- format(round(tab.trt$Estimate, digits=2), nsmall = 2)
tab.trt$SE <- format(round(tab.trt$SE, digits=2), nsmall = 2)
tab.trt$t <- format(round(tab.trt$t, digits=2), nsmall = 2)
tab.trt$p <- format.pval(tab.trt$p, digits = 1,eps = 0.001)


#reorder rows
tab.trt <- tab.trt[c(1:10,12:19,11,20:27),]

#add effect names
tab.trt$Effect <- c("Intercept (Ref. Post-hatching - Control)",
                    "Achillea",
                    "Pulicaria",
                    "Mentha",
                    "Lavandula",
                    "Helichrysum",
                    "Lavandula - Mentha",
                    "Lavandula - Helichrysum",
                    "Mentha - Helichrysum",
                    "Control",
                    "Achillea",
                    "Pulicaria",
                    "Mentha",
                    "Lavandula",
                    "Helichrysum",
                    "Lavandula - Mentha",
                    "Lavandula - Helichrysum",
                    "Mentha - Helichrysum",
                    "Control",
                    "Achillea",
                    "Pulicaria",
                    "Mentha",
                    "Lavandula",
                    "Helichrysum",
                    "Lavandula - Mentha",
                    "Lavandula - Helichrysum",
                    "Mentha - Helichrysum")

#reorder
tab.trt <- select(tab.trt,
                  Effect, Estimate, SE, t, p, CI)

#add random effect
tab.trt[28,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.trt[29,] <- c("Random", "", "", "", "", "")
tab.trt[30,] <- c("Nest", 
                  format(round(rpt.trt.adj$R$nichoir, digits=2), nsmall = 2),
                  format(round(rpt.trt.adj$se[1,], digits=2), nsmall = 2),
                  "", 
                  format.pval(rpt.trt.adj$P[1,1], digits = 1,eps = 0.001),
                  paste("[",
                        format(round(rpt.trt.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                        ", ",
                        format(round(rpt.trt.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                        "]",
                        sep = ""))
tab.trt[31,] <- c("Residual", 
                  format(round(rpt.trt.adj$R$Residual, digits=2), nsmall = 2),
                  "", "", "", "")

flex.trt <- flextable(tab.trt)

flex.trt <- flex.trt %>% 
  set_header_labels(Effect = "Effect",
                    Estimate = "Estimate",
                    SE = "SE",
                    t = "t-value",
                    p = "p-value",
                    CI = "95% CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()
flex.trt


save_as_docx(
  "Table S3.1" = flex.trt,
  path = "table/Table_S3.1.docx")
