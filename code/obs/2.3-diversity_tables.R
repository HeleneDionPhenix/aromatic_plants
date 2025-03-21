# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Diversity analysis - Table creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03


# Librairies --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(flextable)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/obs/2-output_models_diversity.RData")

# Nest bacterial diversity -----------------------------------------------

summary.shannon.nest <- summary(mod.shannon.nest)

tab.shannon.nest <- data.frame(summary.shannon.nest$coefficients)

tab.shannon.nest <- select(tab.shannon.nest,
                  -df)

#add CI
confint.shannon.nest <- data.frame(confint.nest)
names(confint.shannon.nest) <- c("L", "H")
confint.shannon.nest$CI <- paste("[",
                        format(round(confint.shannon.nest$L, digits=2), nsmall = 2),
                        ", ",
                        format(round(confint.shannon.nest$H, digits=2), nsmall = 2),
                        "]",
                        sep = "")
confint.shannon.nest <- confint.shannon.nest[rownames(tab.shannon.nest),]
tab.shannon.nest$CI <- confint.shannon.nest$CI

#rename
names(tab.shannon.nest) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.shannon.nest$Estimate <- format(round(tab.shannon.nest$Estimate, digits=2), nsmall = 2)
tab.shannon.nest$SE <- format(round(tab.shannon.nest$SE, digits=2), nsmall = 2)
tab.shannon.nest$t <- format(round(tab.shannon.nest$t, digits=2), nsmall = 2)
tab.shannon.nest$p <- format.pval(tab.shannon.nest$p, digits = 1,eps = 0.001)

tab.shannon.nest$Effect <- c("Intercept (Ref. EPirio - Post-hatching",
                             "Aromatic plant quantity",
                             "DMuro",
                             "EMuro",
                             "Incubation",
                             "Aromatic plant quantity:DMuro",
                             "Aromatic plant quantity:EMuro",
                             "Aromatic plant quantity:Incubation",
                             "DMuro:Incubation",
                             "EMuro:Incubation",
                             "Aromatic plant quantity:DMuro:Incubation",
                             "Aromatic plant quantity:EMuro:Incubation")

tab.shannon.nest <- select(tab.shannon.nest,
                           Effect, Estimate, SE, t, p, CI)

#add random effect
tab.shannon.nest[13,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.shannon.nest[14,] <- c("Random", "", "", "", "", "")
tab.shannon.nest[15,] <- c("Nest", 
                  format(round(rpt.shannon.nest.adj$R$clutch, digits=2), nsmall = 2),
                  format(round(rpt.shannon.nest.adj$se[1,], digits=2), nsmall = 2),
                  "", 
                  format.pval(rpt.shannon.nest.adj$P[1,1], digits = 1,eps = 0.001),
                  paste("[",
                        format(round(rpt.shannon.nest.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                        ", ",
                        format(round(rpt.shannon.nest.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                        "]",
                        sep = ""))
tab.shannon.nest[16,] <- c("Residual", 
                  format(round(rpt.shannon.nest.adj$R$Residual, digits=2), nsmall = 2),
                  "", "", "", "")

# Eggshell bacterial diversity --------------------------------------------

summary.shannon.eggshell <- summary(mod.shannon.eggshell)

tab.shannon.eggshell <- data.frame(summary.shannon.eggshell$coefficients)

#add CI
confint.shannon.eggshell <- data.frame(confint.eggshell)
names(confint.shannon.eggshell) <- c("L", "H")
confint.shannon.eggshell$CI <- paste("[",
                                     format(round(confint.shannon.eggshell$L, digits=2), nsmall = 2),
                                     ", ",
                                     format(round(confint.shannon.eggshell$H, digits=2), nsmall = 2),
                                     "]",
                                     sep = "")
confint.shannon.eggshell <- confint.shannon.eggshell[rownames(tab.shannon.eggshell),]
tab.shannon.eggshell$CI <- confint.shannon.eggshell$CI

#rename
names(tab.shannon.eggshell) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.shannon.eggshell$Estimate <- format(round(tab.shannon.eggshell$Estimate, digits=2), nsmall = 2)
tab.shannon.eggshell$SE <- format(round(tab.shannon.eggshell$SE, digits=2), nsmall = 2)
tab.shannon.eggshell$t <- format(round(tab.shannon.eggshell$t, digits=2), nsmall = 2)
tab.shannon.eggshell$p <- format.pval(tab.shannon.eggshell$p, digits = 1,eps = 0.001)

tab.shannon.eggshell$Effect <- c("Intercept (Ref. EPirio)",
                                 "Aromatic plant quantity",
                                 "DMuro",
                                 "EMuro",
                                 "Aromatic plant quantity:DMuro",
                                 "Aromatic plant quantity:EMuro")

tab.shannon.eggshell <- select(tab.shannon.eggshell,
                               Effect, Estimate, SE, t, p, CI)


#combine table
tab.shannon <- rbind(tab.shannon.nest, tab.shannon.eggshell)

tab.f <- flextable(tab.shannon) %>%
  set_header_labels(Effect = "Effect",
                    Estimate = "Estimate",
                    SE = "SE",
                    t = "t-value",
                    p = "p-value",
                    CI = "95% CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()

tab.f

save_as_docx(
  "Table S3.2" = tab.f,
  path = "table/Table_S3.2.docx")

# Amount of aromatic plant between pop ------------------------------------

summary.aro.pop <- summary(mod.aro.pop)

tab.aro.pop <- data.frame(summary.aro.pop$coefficients)

tab.aro.pop <- select(tab.aro.pop,
                               -df)

#add CI
confint.aro.pop <- data.frame(confint.aro.pop)
names(confint.aro.pop) <- c("L", "H")
confint.aro.pop$CI <- paste("[",
                                     format(round(confint.aro.pop$L, digits=2), nsmall = 2),
                                     ", ",
                                     format(round(confint.aro.pop$H, digits=2), nsmall = 2),
                                     "]",
                                     sep = "")
confint.aro.pop <- confint.aro.pop[rownames(tab.aro.pop),]
tab.aro.pop$CI <- confint.aro.pop$CI

#rename
names(tab.aro.pop) <- c("Estimate", "SE", "t", "p", "CI")

#round
tab.aro.pop$Estimate <- format(round(tab.aro.pop$Estimate, digits=2), nsmall = 2)
tab.aro.pop$SE <- format(round(tab.aro.pop$SE, digits=2), nsmall = 2)
tab.aro.pop$t <- format(round(tab.aro.pop$t, digits=2), nsmall = 2)
tab.aro.pop$p <- format.pval(tab.aro.pop$p, digits = 1,eps = 0.001)

tab.aro.pop$Effect <- c("Intercept (Ref. DMuro - Post-hatching)",
                        "EMuro",
                        "EPirio",
                        "Incubation",
                        "EMuro:Incubation",
                        "EPirio:Incubation")

tab.aro.pop <- select(tab.aro.pop,
                               Effect, Estimate, SE, t, p, CI)

#add random effect
tab.aro.pop[7,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.aro.pop[8,] <- c("Random", "", "", "", "", "")
tab.aro.pop[9,] <- c("Nest", 
                           format(round(rpt.aro.pop.adj$R$clutch, digits=2), nsmall = 2),
                           format(round(rpt.aro.pop.adj$se[1,], digits=2), nsmall = 2),
                           "", 
                           format.pval(rpt.aro.pop.adj$P[1,1], digits = 1,eps = 0.001),
                           paste("[",
                                 format(round(rpt.aro.pop.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                                 ", ",
                                 format(round(rpt.aro.pop.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                                 "]",
                                 sep = ""))
tab.aro.pop[10,] <- c("Residual", 
                           format(round(rpt.aro.pop.adj$R$Residual, digits=2), nsmall = 2),
                           "", "", "", "")


tab.f <- flextable(tab.aro.pop) %>%
  set_header_labels(Effect = "Effect",
                    Estimate = "Estimate",
                    SE = "SE",
                    t = "t-value",
                    p = "p-value",
                    CI = "95% CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()

tab.f

save_as_docx(
  "Table S3.3" = tab.f,
  path = "table/Table_S3.3.docx")

