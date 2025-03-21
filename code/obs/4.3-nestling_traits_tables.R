# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nestling traits
# Nestling traits tables
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/obs/4-output_models_nestling_traits.RData")


# PCoA 1 --------------------------------------------------------------

summary.1 <- summary(mod.1)

tab.1 <- data.frame(summary.1$coefficients)

tab.1 <- select(tab.1,
                  -df)

#add CI
confint.1 <- data.frame(confint.1)
names(confint.1) <- c("L", "H")
confint.1$CI <- paste("[",
                        format(round(confint.1$L, digits=2), nsmall = 2),
                        ", ",
                        format(round(confint.1$H, digits=2), nsmall = 2),
                        "]",
                        sep = "")
confint.1 <- confint.1[rownames(tab.1),]
tab.1$CI <- confint.1$CI

#rename
names(tab.1) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.1$Estimate <- format(round(tab.1$Estimate, digits=2), nsmall = 2)
tab.1$SE <- format(round(tab.1$SE, digits=2), nsmall = 2)
tab.1$t <- format(round(tab.1$t, digits=2), nsmall = 2)
tab.1$p <- format.pval(tab.1$p, digits = 1,eps = 0.001)

tab.1$Effect <- c("Intercept (Ref. DMuro",
                  "EMuro",
                  "EPirio",
                  "Aromatic plants",
                  "Aromatic plant:EMuro",
                  "Aromatic plant:EPirio")

tab.1 <- select(tab.1,
                Effect, Estimate, SE, t, p, CI)


#add random effect
tab.1[7,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.1[8,] <- c("Random", "", "", "", "", "")
tab.1[9,] <- c("Nest", 
                     format(round(rpt.1.adj$R$clutch, digits=2), nsmall = 2),
                     format(round(rpt.1.adj$se[1,], digits=2), nsmall = 2),
                     "", 
                     format.pval(rpt.1.adj$P[1,1], digits = 1,eps = 0.001),
                     paste("[",
                           format(round(rpt.1.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                           ", ",
                           format(round(rpt.1.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                           "]",
                           sep = ""))
tab.1[10,] <- c("Residual", 
                      format(round(rpt.1.adj$R$Residual, digits=2), nsmall = 2),
                      "", "", "", "")


# PCoA 2 -------------------------------------------------------------------

summary.2 <- summary(mod.2)

tab.2 <- data.frame(summary.2$coefficients)

tab.2 <- select(tab.2,
                -df)

#add CI
confint.2 <- data.frame(confint.2)
names(confint.2) <- c("L", "H")
confint.2$CI <- paste("[",
                      format(round(confint.2$L, digits=2), nsmall = 2),
                      ", ",
                      format(round(confint.2$H, digits=2), nsmall = 2),
                      "]",
                      sep = "")
confint.2 <- confint.2[rownames(tab.2),]
tab.2$CI <- confint.2$CI

#rename
names(tab.2) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.2$Estimate <- format(round(tab.2$Estimate, digits=2), nsmall = 2)
tab.2$SE <- format(round(tab.2$SE, digits=2), nsmall = 2)
tab.2$t <- format(round(tab.2$t, digits=2), nsmall = 2)
tab.2$p <- format.pval(tab.2$p, digits = 1,eps = 0.001)

tab.2$Effect <- c("Intercept (Ref. DMuro",
                  "Aromatic plants",
                  "EMuro",
                  "EPirio",
                  "Aromatic plant:EMuro",
                  "Aromatic plant:EPirio")

tab.2 <- select(tab.2,
                Effect, Estimate, SE, t, p, CI)


#add random effect
tab.2[7,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.2[8,] <- c("Random", "", "", "", "", "")
tab.2[9,] <- c("Nest", 
               format(round(rpt.2.adj$R$clutch, digits=2), nsmall = 2),
               format(round(rpt.2.adj$se[1,], digits=2), nsmall = 2),
               "", 
               format.pval(rpt.2.adj$P[1,1], digits = 1,eps = 0.001),
               paste("[",
                     format(round(rpt.2.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                     ", ",
                     format(round(rpt.2.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                     "]",
                     sep = ""))
tab.2[10,] <- c("Residual", 
                format(round(rpt.2.adj$R$Residual, digits=2), nsmall = 2),
                "", "", "", "")


tab.pcoa <- rbind(tab.1, tab.2)

tab.f <- flextable(tab.pcoa) %>%
  autofit()

tab.f

save_as_docx(
  "Table S3.4" = tab.f,
  path = "table/Table_S3.4.docx")
