# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nestling traits
# Nestling traits analysis
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/obs/5-output_models_complementary.RData")

# Table S4.1. Developmental index ~ hatching date --------------------------

summary.devel <- summary(mod.devel)

tab.devel <- data.frame(summary.devel$coefficients)

tab.devel <- select(tab.devel,
                  -df)

#add CI
confint.devel <- data.frame(confint.devel)
names(confint.devel) <- c("L", "H")
confint.devel$CI <- paste("[",
                        format(round(confint.devel$L, digits=2), nsmall = 2),
                        ", ",
                        format(round(confint.devel$H, digits=2), nsmall = 2),
                        "]",
                        sep = "")
confint.devel <- confint.devel[rownames(tab.devel),]
tab.devel$CI <- confint.devel$CI

#rename
names(tab.devel) <- c("Estimate", "SE", "t", "p-value", "95% CI")

#round
tab.devel$Estimate <- format(round(tab.devel$Estimate, digits=2), nsmall = 2)
tab.devel$SE <- format(round(tab.devel$SE, digits=2), nsmall = 2)
tab.devel$t <- format(round(tab.devel$t, digits=2), nsmall = 2)
tab.devel$p <- format.pval(tab.devel$p, digits = 1,eps = 0.001)


#add effect names
tab.devel$Effect <- c("Intercept (Ref. DMuro)",
                    "Hatching date",
                    "EMuro",
                    "EPirio")

#reorder
tab.devel <- select(tab.devel,
                  Effect, Estimate, SE, t, p, CI)


#add random effect
tab.devel[5,] <- c("Effect", "R2 ajusted", "SE", "", "p-value", "CI")
tab.devel[6,] <- c("Random", "", "", "", "", "")
tab.devel[7,] <- c("Nest", 
                  format(round(rpt.devel.adj$R$clutch, digits=2), nsmall = 2),
                  format(round(rpt.devel.adj$se[1,], digits=2), nsmall = 2),
                  "", 
                  format.pval(rpt.devel.adj$P[1,1], digits = 1,eps = 0.001),
                  paste("[",
                        format(round(rpt.devel.adj$CI_emp$`2.5%`[1], digits=2), nsmall = 2),
                        ", ",
                        format(round(rpt.devel.adj$CI_emp$`97.5%`[1], digits=2), nsmall = 2),
                        "]",
                        sep = ""))
tab.devel[8,] <- c("Residual", 
                  format(round(rpt.devel.adj$R$Residual, digits=2), nsmall = 2),
                  "", "", "", "")

flex.devel <- flextable(tab.devel)

flex.devel <- flex.devel %>% 
  set_header_labels(Effect = "Effect",
                    Estimate = "Estimate",
                    SE = "SE",
                    t = "t-value",
                    p = "p-value",
                    CI = "95% CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()
flex.devel

save_as_docx(
  "Table S4.1" = flex.devel,
  path = "table/Table_S4.1.docx")

# Table S4.2. Aromatic plant quantity ~ hatching date ----------------------

summary.aro <- summary(mod.aro)

tab.aro <- data.frame(summary.aro$coefficients)


#add CI
confint.aro <- data.frame(confint.aro)
names(confint.aro) <- c("L", "H")
confint.aro$CI <- paste("[",
                          format(round(confint.aro$L, digits=2), nsmall = 2),
                          ", ",
                          format(round(confint.aro$H, digits=2), nsmall = 2),
                          "]",
                          sep = "")
confint.aro <- confint.aro[rownames(tab.aro),]
tab.aro$CI <- confint.aro$CI

#rename
names(tab.aro) <- c("Estimate", "SE", "t", "p-value", "CI")

#round
tab.aro$Estimate <- format(round(tab.aro$Estimate, digits=2), nsmall = 2)
tab.aro$SE <- format(round(tab.aro$SE, digits=2), nsmall = 2)
tab.aro$t <- format(round(tab.aro$t, digits=2), nsmall = 2)
tab.aro$p <- format.pval(tab.aro$p, digits = 1,eps = 0.001)


#add effect names
tab.aro$Effect <- c("Intercept (Ref. DMuro)",
                      "Hatching date",
                      "EMuro",
                      "EPirio")

#reorder
tab.aro <- select(tab.aro,
                    Effect, Estimate, SE, t, p, CI)




flex.aro <- flextable(tab.aro)

flex.aro <- flex.aro %>% 
  set_header_labels(Effect = "Effect",
                    Estimate = "Estimate",
                    SE = "SE",
                    t = "t-value",
                    p = "p-value",
                    CI = "95% CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()
flex.aro

save_as_docx(
  "Table S4.2" = flex.aro,
  path = "table/Table_S4.2.docx")

# Table S4.3. Fledging success over the years ------------------------------

f.dm <- filter(fledg.tab, pop == "DMuro")
f.em <- filter(fledg.tab, pop == "EMuro")
f.ep <- filter(fledg.tab, pop == "EPirio")

f.dm <- select(f.dm,
               -pop)
f.em <- select(f.em,
               -pop, -an)
f.ep <- select(f.ep,
               -pop, -an)

tab.fledg <- data.frame(an = f.dm$an,
                        f.dm = format(round(f.dm$fledg, digits=2), nsmall = 2),
                        CI.dm = paste("[",
                                      format(round(f.dm$lowerCI, digits=2), nsmall = 2),
                                      ", ",
                                      format(round(f.dm$upperCI, digits=2), nsmall = 2),
                                      "]",
                                      sep = ""),
                        f.em = format(round(f.em$fledg, digits=2), nsmall = 2),
                        CI.em = paste("[",
                                      format(round(f.em$lowerCI, digits=2), nsmall = 2),
                                      ", ",
                                      format(round(f.em$upperCI, digits=2), nsmall = 2),
                                      "]",
                                      sep = ""),
                        f.ep = format(round(f.ep$fledg, digits=2), nsmall = 2),
                        CI.ep = paste("[",
                                      format(round(f.ep$lowerCI, digits=2), nsmall = 2),
                                      ", ",
                                      format(round(f.ep$upperCI, digits=2), nsmall = 2),
                                      "]",
                                      sep = ""))

tab.fledg$an <- as.character(tab.fledg$an)

tab.fledg[21,] <- c("Mean", 
                    format(round(mean.fledg.tab$fledg[1], digits=2), nsmall = 2),
                    paste("[",
                          format(round(mean.fledg.tab$lowerCI[1], digits=2), nsmall = 2),
                          ", ",
                          format(round(mean.fledg.tab$upperCI[1], digits=2), nsmall = 2),
                          "]",
                          sep = ""),
                    format(round(mean.fledg.tab$fledg[2], digits=2), nsmall = 2),
                    paste("[",
                          format(round(mean.fledg.tab$lowerCI[2], digits=2), nsmall = 2),
                          ", ",
                          format(round(mean.fledg.tab$upperCI[2], digits=2), nsmall = 2),
                          "]",
                          sep = ""),
                    format(round(mean.fledg.tab$fledg[3], digits=2), nsmall = 2),
                    paste("[",
                          format(round(mean.fledg.tab$lowerCI[3], digits=2), nsmall = 2),
                          ", ",
                          format(round(mean.fledg.tab$upperCI[3], digits=2), nsmall = 2),
                          "]",
                          sep = ""))

tab.fledg <- tab.fledg[c(21, 1:20),]                    

tab.f <- flextable(tab.fledg) %>%
  add_header_row(values = c("", "DMuro", "EMuro", "EPirio"),
                 colwidths = c(1,2,2,2), top = TRUE) %>%
  set_header_labels(an = "Year",
                    f.dm = "Fledging success",
                    CI.dm = "CI",
                    f.em = "Fledging success",
                    CI.em = "CI",
                    f.ep = "Fledging success",
                    CI.ep = "CI") %>%
  font(fontname='Times New Roman') %>%
  autofit()

tab.f

save_as_docx(
  "Table S4.3" = tab.f,
  path = "table/Table_S4.3.docx")
