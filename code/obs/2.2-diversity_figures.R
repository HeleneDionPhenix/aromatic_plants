# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Diversity analysis - Figure creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/obs/2-output_models_diversity.RData")

# Nest --------------------------------------------------------------------

newdat.nest <- expand.grid(ARO.s = seq(min(compo.stage$ARO.s),
                                  max(compo.stage$ARO.s),
                                  length = 100),
                           pop = levels(compo.stage$pop),
                           stade = levels(compo.stage$stade),
                           clutch = 1)

pred.shannon.nest <- predict(mod.shannon.nest,
                             newdat.nest,
                             type = "link",
                             se.fit = T,
                             allow.new.levels =T,
                             re.form=NA)

pred.shannon.nest <- data.frame(newdat.nest,
                                shannon = pred.shannon.nest$fit,
                                lowerCI = pred.shannon.nest$fit+(qnorm(0.025)*pred.shannon.nest$se.fit),
                                upperCI = pred.shannon.nest$fit+(qnorm(0.975)*pred.shannon.nest$se.fit))

#range per pop and stage
range(compo.stage$ARO.s)
DMC <- range(compo.stage[compo.stage$pop == "DMuro" & compo.stage$stade == "couvaison", "ARO.s"])
DME <- range(compo.stage[compo.stage$pop == "DMuro" & compo.stage$stade == "elevage", "ARO.s"])
EMC <- range(compo.stage[compo.stage$pop == "EMuro" & compo.stage$stade == "couvaison", "ARO.s"])
EME <- range(compo.stage[compo.stage$pop == "EMuro" & compo.stage$stade == "elevage", "ARO.s"])
EPC <- range(compo.stage[compo.stage$pop == "EPirio" & compo.stage$stade == "couvaison", "ARO.s"])
EPE <- range(compo.stage[compo.stage$pop == "EPirio" & compo.stage$stade == "elevage", "ARO.s"])

#Truncate predicted values to avoid extrapolation
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "DMuro" & stade == "couvaison" & ARO.s > max(DMC)))
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "DMuro" & stade == "elevage" & ARO.s > max(DME)))
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "EMuro" & stade == "couvaison" & ARO.s > max(EMC)))
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "EMuro" & stade == "elevage" & ARO.s > max(EME)))
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "EPirio" & stade == "couvaison" & ARO.s > max(EPC)))
pred.shannon.nest <- filter(pred.shannon.nest, !(pop == "EPirio" & stade == "elevage" & ARO.s > max(EPE)))

#unscale ARO.log
pred.shannon.nest$ARO.log <- pred.shannon.nest$ARO.s*sd(compo.stage$ARO.log) + mean(compo.stage$ARO.log)

#order pop
pred.shannon.nest$pop <- factor(pred.shannon.nest$pop, levels=c("DMuro", "EMuro", "EPirio"))
compo.stage$pop <- factor(compo.stage$pop, levels=c("DMuro", "EMuro", "EPirio"))

#plot per stage
pred.shannon.nest.c <- filter(pred.shannon.nest, stade == "couvaison")
pred.shannon.nest.e <- filter(pred.shannon.nest, stade == "elevage")
compo.couv <- filter(compo.stage, stade == "couvaison")
compo.ele <- filter(compo.stage, stade == "elevage")

## Nest in incubation ------------------------------------------------------

p.shannon.nest.c <- ggplot(data = pred.shannon.nest.c, 
                      aes(x = ARO.log,
                          y = shannon,
                          group = pop)) +
  geom_point(data = compo.couv,
             aes(x = ARO.log,
                 y = shannon,
                 color = pop,
                 shape = pop),
             size = 1) +
  
  geom_ribbon(aes(x = ARO.log,
                  y = shannon,
                  ymin = lowerCI,
                  ymax = upperCI,
                  fill = pop),
              alpha = 0.3) +
  geom_line(aes(x = ARO.log,
                y = shannon,
                color = pop),
            linetype = 5,
            size=0.6) +
  
  scale_color_manual(name = "",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  scale_fill_manual(name = "",
                    labels = c(EMuro = "EMuro",
                               DMuro = "DMuro",
                               EPirio = "EPirio"),
                    values = c("DMuro" = "steelblue1",
                               "EMuro" = "green3",
                               "EPirio" = "forestgreen"))+
  scale_shape_manual(name = "",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = 15,
                                "EMuro" = 16,
                                "EPirio" = 17))+
  annotate("text",
           x = 2.5,
           y = 2,
           label = paste("R2=",
                         round(rpt.shannon.nest$R$Fixed, 2),
                         sep = ""),
           family = "Times New Roman")+

  ggtitle("Nest microbiota - Incubation") +
  ylab("Bacterial Shannon diversity") + 
  xlab("Quantity of aromatic plants (log(%))") + 
  
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.15,0.2),
        legend.text = element_text(size=12), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))

p.shannon.nest.c

## Nest in post-hatching ---------------------------------------------------

p.shannon.nest.e <- ggplot(data = pred.shannon.nest.e, 
                      aes(x = ARO.log,
                          y = shannon,
                          group = pop)) +
  geom_point(data = compo.ele,
             aes(x = ARO.log,
                 y = shannon,
                 color = pop,
                 shape = pop),
             size = 1) +
  
  geom_ribbon(aes(x = ARO.log,
                  y = shannon,
                  ymin = lowerCI,
                  ymax = upperCI,
                  fill = pop),
              alpha = 0.3) +
  geom_line(aes(x = ARO.log,
                y = shannon,
                color = pop,
                linetype = pop),
            size=0.6) +
  
  scale_color_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  scale_shape_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = 15,
                                "EMuro" = 16,
                                "EPirio" = 17))+
  scale_fill_manual(name = "Population",
                    labels = c(EMuro = "EMuro",
                               DMuro = "DMuro",
                               EPirio = "EPirio"),
                    values = c("DMuro" = "steelblue1",
                               "EMuro" = "green3",
                               "EPirio" = "forestgreen"))+
  scale_linetype_manual(name = "Population",
                        labels = c(EMuro = "EMuro",
                                   DMuro = "DMuro",
                                   EPirio = "EPirio"),
                        values = c("DMuro" = 5,
                                   "EMuro" = 5,
                                   "EPirio" = 1))+
  ggtitle("Nest microbiota - Post-hatching") +
  ylab("Bacterial Shannon diversity") + 
  xlab("Quantity of aromatic plants (log(%))") + 
  
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position= "none")
p.shannon.nest.e


# Eggshell ----------------------------------------------------------------

newdat.egg <- expand.grid(ARO.s = seq(min(compo.incu$ARO.s),
                                  max(compo.incu$ARO.s),
                                  length = 100),
                          pop = levels(compo.incu$pop))


pred.shannon.egg <- predict(mod.shannon.eggshell,
                            newdat.egg,
                            se.fit = T)

pred.shannon.egg <- data.frame(newdat.egg,
                               shannon = pred.shannon.egg$fit,
                               lowerCI = pred.shannon.egg$fit+(qnorm(0.025)*pred.shannon.egg$se.fit),
                               upperCI = pred.shannon.egg$fit+(qnorm(0.975)*pred.shannon.egg$se.fit))

#range per pop and stage
range(compo.incu$ARO.s)
DMC <- range(compo.incu[compo.incu$pop == "DMuro", "ARO.s"])
EMC <- range(compo.incu[compo.incu$pop == "EMuro", "ARO.s"])
EPC <- range(compo.incu[compo.incu$pop == "EPirio", "ARO.s"])

#Truncate predicted values to avoid extrapolation
pred.shannon.egg <- filter(pred.shannon.egg, !(pop == "DMuro" & ARO.s > max(DMC)))
pred.shannon.egg <- filter(pred.shannon.egg, !(pop == "EMuro" & ARO.s > max(EMC)))

#Unscale ARO.log
pred.shannon.egg$ARO.log <- pred.shannon.egg$ARO.s * sd(compo.incu$ARO.log) + mean(compo.incu$ARO.log)

#order pop
pred.shannon.egg$pop <- factor(pred.shannon.egg$pop, levels=c("DMuro", "EMuro", "EPirio"))
compo.incu$pop <- factor(compo.incu$pop, levels=c("DMuro", "EMuro", "EPirio"))

p.shannon.c.egg <- ggplot(data = pred.shannon.egg, 
                          aes(x = ARO.log,
                              y = shannon,
                              group = pop)) +
  geom_point(data = compo.incu,
             aes(x = ARO.log,
                 y = shannon,
                 color = pop,
                 shape = pop),
             size = 1) +
  
  geom_ribbon(aes(x = ARO.log,
                  y = shannon,
                  ymin = lowerCI,
                  ymax = upperCI,
                  fill = pop),
              alpha = 0.3) +
  geom_line(aes(x = ARO.log,
                y = shannon,
                color = pop,
                linetype = pop),
            size=0.6) +
  
  scale_color_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  scale_fill_manual(name = "Population",
                    labels = c(EMuro = "EMuro",
                               DMuro = "DMuro",
                               EPirio = "EPirio"),
                    values = c("DMuro" = "steelblue1",
                               "EMuro" = "green3",
                               "EPirio" = "forestgreen"))+
  scale_shape_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = 15,
                                "EMuro" = 16,
                                "EPirio" = 17))+
  scale_linetype_manual(name = "Population",
                        labels = c(EMuro = "EMuro",
                                   DMuro = "DMuro",
                                   EPirio = "EPirio"),
                        values = c("DMuro" = 5,
                                   "EMuro" = 5,
                                   "EPirio" = 1))+
  annotate("text",
           x = 2.5,
           y = 4.2,
           label = paste("R2=",
                         round(summary(mod.shannon.eggshell)$adj.r.squared, 2),
                         sep = ""),
           family = "Times New Roman") +
  
  ggtitle("Eggshell microbiota - Incubation") +
  ylab("Bacterial Shannon diversity") + 
  xlab("Quantity of aromatic plants (log(%))") + 
  
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "none",
        legend.text = element_text(size=12), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))

p.shannon.c.egg

# Difference in quantity of aromatic plant between population ---------------

newdat.aro <- expand.grid(pop = c("EPirio",
                                  "EMuro", 
                                  "DMuro"),
                          stade = levels(compo$stade),
                          clutch = 1)

pred.aro <- predict(mod.aro.pop,
                newdat.aro,
                type = "link",
                se.fit = T,
                allow.new.levels =T,
                re.form=NA)

newdat.aro$ARO.log <- pred.aro$fit
newdat.aro$lowerIC <- pred.aro$fit + qnorm(0.025) * pred.aro$se.fit
newdat.aro$upperIC <- pred.aro$fit + qnorm(0.975) * pred.aro$se.fit

newdat.aro$stade <- relevel(newdat.aro$stade, ref = "couvaison")

p.pop.aro <- ggplot(newdat.aro,
                     aes(y = ARO.log,
                         x = pop,
                         group = stade)) +
  # geom_jitter(data = compo.couv,
  #            aes(color = pop),
  #            alpha = 0.5) +
  geom_pointrange(aes(ymin = lowerIC,
                      ymax = upperIC,
                      color = pop,
                      shape = stade), 
                  size=0.6, 
                  position = position_dodge(width = 0.5)) +

  scale_color_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"),
                     breaks = c(DMuro = "DMuro",
                                EMuro = "EMuro",
                                EPirio = "EPirio"))+
  scale_shape_manual(name = "",
                     labels = c(couvaison = "Incubation",
                                elevage = "Post-hatching"),
                     values = c(couvaison = 16,
                                elevage = 17))+
  annotate("text",
           x = 0.8,
           y = 1.8,
           label = paste("R2=",
                         round(rpt.aro.pop$R$Fixed, 2),
                         sep = ""),
           family = "Times New Roman")+
  ylim(0.5, 2)+
  ylab("Quantity of aromatic plants (log(%))")+
  xlab("") +
  coord_flip() +
  theme_minimal(base_size = 12,
                base_family = "Times New Roman")+
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position= "top",
        legend.text = element_text(size=12), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent", 
                                         color = "transparent"),
        plot.background = element_rect(fill = "white", 
                                        color = "white"))+
  guides(colour = "none")


p.pop.aro


# Figure composite -------------------------------------------------------

p <- ggarrange(p.shannon.nest.c, p.shannon.nest.e,
               p.shannon.c.egg, p.pop.aro,
                       ncol = 2, nrow = 2,
                       labels = c("A", "B", "C", "D"))
p

# Export ------------------------------------------------------------------

ggsave(
  "figure/Figure_4.png",
  p,
  width = 18,
  height = 18,
  units = "cm",
  dpi = 300
)
