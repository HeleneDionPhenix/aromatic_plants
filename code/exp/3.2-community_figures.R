# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Community analysis - Figure creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-19

# Libraries --------------------------------------------------------------

library(dplyr)
library(vegan)
library(ggplot2)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/exp/3-output_models_community.RData")


# Plant microbiota composition --------------------------------------------

## RDA ---------------------------------------------------------------------

perc <- round(100*(summary(mod.pl)$cont$importance[2,1:2]), 2)

#extract sites scores
sites <- scores(mod.pl, 
                choices = c(1,2), 
                display = "sites",
                scaling = 2)
sites <- data.frame(sites)

#create dataset for plotting
data.pl <- data.plant[rownames(sites),]
data.pl <- cbind(sites, data.pl)


p.pl <- ggplot(data.pl, aes(x = CAP1, y = CAP2)) +
  geom_polygon(aes(color = traitement),
               fill = "transparent") +
  
  geom_point(aes(color = traitement,
                 fill = traitement,
                 shape = traitement),
             alpha = 1) +
  annotate("text",
           x = -1,
           y = 2.3,
           size = 5,
           label = paste("adjusted R2=", 
                         round(RsquareAdj(mod.pl)$adj.r.squared, 2),
                         " p=",
                         round(p.value.pl, 3),
                         sep = ""),
           family = "Times New Roman") +
  scale_color_manual(name = "Species",
                     labels = c(ach = expression(italic("Achillea ligustica")),
                                pul = expression(italic("Pulicaria odorata")), 
                                men = expression(italic("Mentha suaveolens")), 
                                lav = expression(italic("Lavandula stoechas")), 
                                imm = expression(italic("Helichryum italicum"))),
                     values = c(ach = "#044B7F",
                                pul = "#CC3F0C", 
                                men = "#107E7D", 
                                lav = "#706993", 
                                imm = "#C8AD55")) +
  scale_fill_manual(name = "Species",
                    labels = c(ach = expression(italic("Achillea ligustica")),
                               pul = expression(italic("Pulicaria odorata")), 
                               men = expression(italic("Mentha suaveolens")), 
                               lav = expression(italic("Lavandula stoechas")), 
                               imm = expression(italic("Helichryum italicum"))),
                    values = c(ach = "#044B7F",
                               pul = "#CC3F0C", 
                               men = "#107E7D", 
                               lav = "#706993", 
                               imm = "#C8AD55")) +
  scale_shape_manual(name = "Species",
                     labels = c(ach = expression(italic("Achillea ligustica")),
                                pul = expression(italic("Pulicaria odorata")), 
                                men = expression(italic("Mentha suaveolens")), 
                                lav = expression(italic("Lavandula stoechas")), 
                                imm = expression(italic("Helichryum italicum"))),
                     values = c(ach = 21,
                                pul = 22, 
                                men = 23, 
                                lav = 24, 
                                imm = 25)) +
  
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("RDA 2 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position= c(0.8, 0.85), 
        legend.text = element_text(size=14), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))

p.pl

ggsave("figure/Figure_S2.2.png",
       p.pl,
       dpi = 300,
       width = 8,
       height = 7)


# Nest microbiota composition ---------------------------------------------

perc <- round(100*(summary(mod.trt)$cont$importance[2,1:2]), 2)

#extract sites scores
sites <- scores(mod.trt, 
                choices = c(1,2), 
                display = "sites",
                scaling = 2)
sites <- data.frame(sites)

#create dataset for plotting
data.trt <- data.nest[rownames(sites),]
data.trt <- cbind(sites, data.trt)

pg <- ggplot(data.trt, aes(x = CAP1, y = CAP2)) +
  geom_line(aes(group = nichoir,
                color = stage_pheno), 
            lty = 1, 
            alpha = 0.5) +

  geom_point(aes(color = stage_pheno,
                 shape = trtctrl),
             alpha = 1) +
  annotate("text",
           x = -2.2,
           y = 2.3,
           size = 5,
           label = paste("adjusted R2=", 
                         round(RsquareAdj(mod.trt)$adj.r.squared, 2),
                         " p=",
                         round(p.value.trt, 3),
                         sep = ""),
           family = "Times New Roman") +
  
  scale_color_manual(name = "Phenologic stage",
                     breaks = c("prebrood",
                                "brood",
                                "breed"),
                     labels = c(breed = "Post-hatching",
                                brood = "Incubation",
                                prebrood = "Pre-incubation"),
                     values = c(prebrood = "#F2BB05",
                                brood = "#D74E09",
                                breed = "#6E0E0A")) +
  scale_shape_manual(name = "Treatment",
                     labels = c(ctr = "Control",
                                trt = "Aromatic plant"),
                     values = c(ctr = 0,
                                trt = 20)) +
  
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("RDA 2 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position= "left", 
        legend.text = element_text(size=14), 
        legend.title = element_text(size=14), 
        legend.key.height = unit(0.3, "cm"), 
        legend.spacing.y = unit(0, 'cm'),
        legend.background = element_rect(fill = "transparent"))
pg

ps <- ggplot(data.trt, aes(x = CAP1, y = CAP2)) +
  geom_point(aes(color = stage_pheno,
                 shape = trtctrl),
             alpha = 1) +
  geom_line(aes(group = nichoir,
                color = stage_pheno), 
            lty = 1, 
            alpha = 0.5) +
  facet_wrap(~ factor(stage_pheno,
                      levels = c("prebrood",
                                 "brood",
                                 "breed"),
                      labels = c("Pre-incubation",
                                 "Incubation",
                                 "Post-hatching")),
             scales = "free")+
  scale_color_manual(name = "Phenologic stage",
                     labels = c(breed = "Post-hatching",
                                brood = "Incubation",
                                prebrood = "Pre-incubation"),
                     values = c(prebrood = "#F2BB05",
                                brood = "#D74E09",
                                breed = "#6E0E0A")) +
  scale_shape_manual(name = "Treatment",
                     labels = c(ctr = "Control",
                                trt = "Aromatic plant"),
                     values = c(ctr = 0,
                                trt = 20)) +
  
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("RDA 2 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position= "none", 
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size=14))
ps

p <- ggarrange(pg, ps,
               nrow =2,
               labels = c("A", "B"))
p

ggsave("figure/Figure_3.png",
       p,
       width = 8.5,
       height = 7,
       dpi = 300)
