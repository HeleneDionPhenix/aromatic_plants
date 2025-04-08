# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Diversity analysis - Figure creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-18

# Libraries --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(AICcmodavg)
library(ggpubr)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/2-output_models_diversity.RData")

# Predicted data ----------------------------------------------------------

newdat <- expand.grid(traitement = levels(data.nest$traitement),
                      stage_pheno = levels(data.nest$stage_pheno),
                      nichoir = "1")

pred <- predictSE(mod.trt, 
                  newdat, 
                  allow.new.levels =T, 
                  type = "link", 
                  re.form=NA, 
                  se.fit = T)

newdat=data.frame(newdat, 
                  shannon = pred$fit, 
                  lowerCI=pred$fit+(qnorm(0.025)*pred$se.fit), 
                  upperCI=pred$fit+(qnorm(0.975)*pred$se.fit))

# Plot --------------------------------------------------------------------

# To avoid overlap between stage
data.nest$pos.line <- 0
data.nest[which(data.nest$stage_pheno == "prebrood"), "pos.line"] <- +0.15
data.nest[which(data.nest$stage_pheno == "breed"), "pos.line"] <- -0.15


p.eff <- ggplot(newdat, 
                aes(x = traitement, 
                    y = shannon, 
                    group = stage_pheno))+

  geom_pointrange(aes(ymin = lowerCI, 
                      ymax = upperCI, 
                      color = stage_pheno, 
                      shape = stage_pheno), 
                  size=0.6, 
                  position = position_dodge(width = 0.5)) +
  
  scale_shape_manual(name = "Phenologic stage",
                     labels = c(prebrood = "Pre-incubation",
                                brood = "Incubation",
                                breed = "Post-hatching"),
                     values = c(prebrood = 15,
                                brood = 16,
                                breed = 17))+
  scale_color_manual(name = "Phenologic stage",
                     labels = c(prebrood = "Pre-incubation",
                                brood = "Incubation",
                                breed = "Post-hatching"),
                     values = c(prebrood = "#F2BB05",
                                brood = "#D74E09",
                                breed = "#6E0E0A")) +
  annotate("text",
           x = 9,
           y = 0.1,
           size = 5,
          label = paste("R2=",
                        round(rpt.trt$R$Fixed, 2),
                        sep = ""),
          family = "Times New Roman") +
  
  theme_classic(base_size = 16,
                base_family = "Times New Roman") +
  xlab("") + 
  ylim(0,6) +
  ylab("Bacterial Shannon diversity") +
  theme(legend.title = element_blank(),
        legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p.eff

p.raw <- ggplot(newdat, 
                aes(x=traitement, 
                    y =shannon, 
                    group = stage_pheno)) +
  
  geom_line(data = data.nest,
            aes(x = (as.numeric(traitement) + pos.line),
                group = nichoir,
                color = stage_pheno),
            linetype = 1,
            alpha = 0.5)+
  
  geom_jitter(data=data.nest,
              aes(color = stage_pheno, 
                  shape = stage_pheno),
              alpha = 0.8,
              size = 1.5,
              position = position_dodge(width = 0.5)) + 
  
  scale_shape_manual(name = "Phenologic stage",
                     labels = c(prebrood = "Pre-incubation",
                                brood = "Incubation",
                                breed = "Post-hatching"),
                     limits = c("prebrood", "brood", "breed"),
                     values = c(prebrood = 15,
                                brood = 16,
                                breed = 17))+
  scale_color_manual(name = "Phenologic stage",
                     labels = c(prebrood = "Pre-incubation",
                                brood = "Incubation",
                                breed = "Post-hatching"),
                     limits = c("prebrood", "brood", "breed"),
                     values = c(prebrood = "#F2BB05",
                                brood = "#D74E09",
                                breed = "#6E0E0A")) +
  theme_classic(base_size = 16,
                base_family = "Times New Roman") +
  xlab("") + 
  scale_x_discrete(labels = c(("Control"), 
                              expression(italic("Achillea")), 
                              expression(italic("Pulicaria")), 
                              expression(italic("Mentha")), 
                              expression(italic("Lavandula")), 
                              expression(italic("Helichryum")), 
                              expression(italic("Lavandula - Mentha")), 
                              expression(italic("Lavandula - Helichryum")), 
                              expression(italic("Mentha - Helichryum"))))+
  ylim(0,6)+
  ylab("Bacterial Shannon diversity") +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.20),
        legend.background = element_rect(color = "transparent"),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1))

p.raw

p <- ggarrange(p.eff, p.raw,
               ncol = 1,
               heights = c(5,8),
               labels = c("A", "B"))
p

ggsave("figure/Figure_2.png",
  p,
  width = 8.5,
  height = 9,
  dpi = 300)
