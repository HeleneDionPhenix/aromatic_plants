# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest and eggshell bacterial microbiota
# Community analysis - Figure creation
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03

# Librairies --------------------------------------------------------------

library(dplyr)
library(vegan)
library(ggplot2)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/obs/3-output_models_community.RData")

# Nest --------------------------------------------------------------------

perc <- round(100*(summary(mod.nest)$cont$importance[2,1:2]), 2)

#extract sites scores
sites <- scores(mod.nest, 
                choices = c(1,2), 
                display = "sites",
                scaling = 2)
sites <- data.frame(sites)

#create dataset for plotting
data.nest <- compo.stage[rownames(sites),]
data.nest <- cbind(sites, data.nest)

p.nest <- ggplot(data.nest, aes(x = CAP1, y = CAP2, 
                              fill = pop)) +
  geom_point(aes(size = ARO,
                 fill = pop,
                 shape = stade),
             color = "black",
             alpha = 0.8)+
  annotate("text",
           x = 0.2,
           y = -1.75,
           label = paste("adjusted R2=0.40", 
                         #round(RsquareAdj(mod.nest)$adj.r.squared, 3),
                         " p=",
                         round(p.nest, 3),
                         sep = ""),
           family = "Times New Roman")+
  scale_shape_manual(name = " Phenologic stage",
                     labels = c(elevage = "Post-hatching",
                                couvaison = "Incubation"),
                     values = c(couvaison = 21,
                                elevage = 22)) +
  scale_fill_manual(name = " Population",
                    labels = c(DMuro = "DMuro",
                               EMuro = "EMuro",
                               EPirio = "EPirio"),
                    values = c("DMuro" = "steelblue1",
                               "EMuro" = "green3",
                               "EPirio" = "forestgreen")) +
  scale_size(name = " Aromatic plant quantity (%)")+
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("RDA 2 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position = "left",
        legend.text = element_text(size=14), 
        legend.title = element_text(size=14), 
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0, 'cm'),
        legend.background = element_rect(fill = "transparent"))+
  guides(shape = guide_legend(override.aes = list(size =3)),
         fill = guide_legend(override.aes = list(shape = 21, size =3)))


p.nest


p.pop <- ggplot(data.nest, aes(x = CAP1, y = CAP2)) +
  geom_line(aes(group = nichoir), 
            lty = 1, 
            alpha = 0.5) +
  facet_wrap(~pop)+
  geom_point(aes(size = ARO,
                 fill = pop,
                 shape = stade),
             color = "black",
             alpha = 0.7)+
  scale_shape_manual(name = "Phenologic stage",
                     labels = c(elevage = "Post-hatching",
                                couvaison = "Incubation"),
                     values = c(couvaison = 21,
                                elevage = 22)) +
  scale_fill_manual(name = "Population",
                    labels = c(DMuro = "DMuro",
                               EMuro = "EMuro",
                               EPirio = "EPirio"),
                    values = c("DMuro" = "steelblue1",
                               "EMuro" = "green3",
                               "EPirio" = "forestgreen")) +
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("RDA 2 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position = "none",
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size=14))
p.pop

p.compo <- ggarrange(p.nest, p.pop,
                     nrow = 2,
                     labels = c("A", "B"))
p.compo

ggsave("figure/Figure_5.png",
       p.compo,
       width = 9.5,
       height = 7,
       dpi = 300)
