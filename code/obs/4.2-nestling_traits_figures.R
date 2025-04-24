# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nestling traits
# Nestling traits figures
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-03


# Librairies --------------------------------------------------------------

library(dplyr)
library(vegan)
library(ggplot2)

# Import data -------------------------------------------------------------

rm(list = ls())
load("data/obs/4-output_models_nestling_traits.RData")


# PcoA --------------------------------------------------------------------

# % explained by first two axis
perc <- round(100*(pcoa$CA$eig[1:2]/sum(pcoa$CA$eig)), 2)

#extract nestlings scores
sites <- scores(pcoa, scaling = 2)$sites
sites <- data.frame(sites)

#extract traits scores
traits <- scores(pcoa, scaling = 2)$species
traits <- traits*0.4
traits <- data.frame(traits)
traits$lab <- c("body mass", "tarsus length", "handling aggression", "developmental index")

#create dataset for plotting
rownames(data.cond) <- data.cond$chick_id
data.cond <- data.cond[rownames(sites),]
data.cond <- cbind(sites, data.cond)

p.pcoa <- ggplot(data.cond, aes(x = -MDS1, y = MDS2)) +
  #nestlings
  geom_point(aes(color = pop),
             alpha = 0.5, 
             size = 2.5)+
  
  #traits
  geom_segment(data = traits,
               aes(x = 0,
                   y = 0,
                   xend = -MDS1,
                   yend = MDS2),
               #arrow = arrow(length = unit(0.02, "npc")),
               color = "black",
               alpha = 0.8,
               linetype = 2,
               size = 0.3) +
  geom_text(data = traits,
            aes(x = -MDS1,
                y = MDS2,
                label = lab),
            hjust= c(0.3, 0.4, 0.6, 0.6), 
            vjust= c(-1, 1.4, -0.6, 0.6), 
            size=7,
            color = "black",
            family = "Times New Roman") +
  
  scale_color_manual(name = "",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  annotate("text",
           x = 1.5,
           y = -1.4,
           size = 7,
           label = "n=1149",
           family = "Times New Roman") +

  xlim(-1.5,1.7)+
  xlab(paste("PCoA 1 (", perc[1], "%)", sEM = "")) + 
  ylab(paste("PCoA 2 (", perc[2], "%)", sEM = "")) + 
  theme_bw(base_size = 20,
           base_family = "Times New Roman")+
  theme(legend.position = c(0.90,0.90),
        legend.text = element_text(size=20), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))

p.pcoa

# Model PCoA1 -------------------------------------------------------------

newdat.1 <- expand.grid(ARO.s = seq(min(data.cond$ARO.s),
                                     max(data.cond$ARO.s),
                                     length = 100),
                       pop = levels(data.cond$pop),
                       clutch = 1)

pred.1 <- predict(mod.1, 
                 newdat.1,
                 type = "link",
                 se.fit = T,
                 allow.new.levels =T,
                 re.form=NA)

pred.1 <- data.frame(newdat.1,
                     pcoa1 = pred.1$fit,
                     lowerCI = pred.1$fit+(qnorm(0.025)*pred.1$se.fit),
                     upperCI = pred.1$fit+(qnorm(0.975)*pred.1$se.fit))

#truncate to avoid extrapolation
range(data.cond$ARO.s)
r.DM <- range(data.cond[data.cond$pop=="DMuro","ARO.s"])
r.EM <- range(data.cond[data.cond$pop=="EMuro","ARO.s"])
r.EP <- range(data.cond[data.cond$pop=="EPirio","ARO.s"])

pred.1 <- filter(pred.1, !(pop == "EMuro" & ARO.s > max(r.EM)))
pred.1 <- filter(pred.1, !(pop == "EPirio" & ARO.s > max(r.EP)))

#unscale aromatic quantity
pred.1$ARO.log <- pred.1$ARO.s * sd(data.cond$ARO.log) + mean(data.cond$ARO.log)

p.1 <- ggplot(data = pred.1,
               aes(x=ARO.log, y=pcoa1)) +
  geom_point(data = data.cond,
             aes(x = ARO.log,
                 y = -MDS1,
                 color = pop),
             alpha = 0.5,
             size = 2) +
  
  geom_ribbon(aes(x = ARO.log,
                  y = pcoa1,
                  ymin = lowerCI,
                  ymax = upperCI,
                  fill = pop),
              alpha = 0.3) +
  geom_line(aes(x = ARO.log,
                y = pcoa1,
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
  scale_linetype_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = 1,
                                "EMuro" = 1,
                                "EPirio" = 2)) +
  annotate("text",
           x = 3.3,
           y = -1.4,
           size = 7,
           label = paste("R2=", 
                         round(rpt.1$R$Fixed, 2),
                         sep = ""),
           family = "Times New Roman") +

  ylab("PCoA 1 score") + 
  xlab("Aromatic plant (log(%))") + 
  
  theme_bw(base_size = 20,
           base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        legend.position = "none",
        legend.text = element_text(size=20), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))
p.1  

# create trait arrows 
empty <- data.frame(x = c(0,1),
                    y = c(0,1))
empty$lab <- c("boby mass", "tarsus length")

p.1.leg <- ggplot(empty, aes(x=x, y=y))+
  geom_segment(aes(x = 0.25,
                   y = 0,
                   xend = 0.25,
                   yend = 1),
               arrow = arrow(length = unit(0.2, "npc")),
               color = "black",
               alpha = 1,
               linetype = 1,
               size = 0.3) +
  geom_segment(aes(x = 0.75,
                   y = 0,
                   xend = 0.75,
                   yend = 1),
               arrow = arrow(length = unit(0.2, "npc")),
               color = "black",
               alpha = 1,
               linetype = 1,
               size = 0.3) +
  
  geom_text(aes(x = c(0.15, 0.65),
                y = c(0.5,0.5),
                label = lab),
            size=7,
            color = "black",
            angle = 90,
            family = "Times New Roman") +
  xlim(0,1)+
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))

p.1.leg
  

# Model PCoA 2 ------------------------------------------------------------

newdat.2 <- expand.grid(ARO.s = seq(min(data.cond$ARO.s),
                                      max(data.cond$ARO.s),
                                      length = 100),
                        pop = levels(data.cond$pop),
                        clutch = 1)

pred.2 <- predict(mod.2, 
                  newdat.2,
                  type = "link",
                  se.fit = T,
                  allow.new.levels =T,
                  re.form=NA)

pred.2 <- data.frame(newdat.2,
                     pcoa2 = pred.2$fit,
                     lowerCI = pred.2$fit+(qnorm(0.025)*pred.2$se.fit),
                     upperCI = pred.2$fit+(qnorm(0.975)*pred.2$se.fit))

#truncate to avoid extrapolation
range(data.cond$ARO.s)
r.DM <- range(data.cond[data.cond$pop=="DMuro","ARO.s"])
r.EM <- range(data.cond[data.cond$pop=="EMuro","ARO.s"])
r.EP <- range(data.cond[data.cond$pop=="EPirio","ARO.s"])

pred.2 <- filter(pred.2, !(pop == "EMuro" & ARO.s > max(r.EM)))
pred.2 <- filter(pred.2, !(pop == "EPirio" & ARO.s > max(r.EP)))

#unscale aro
pred.2$ARO.log <- pred.2$ARO.s * sd(data.cond$ARO.log) + mean(data.cond$ARO.log)

p.2 <- ggplot(data = pred.2,
              aes(x=ARO.log, y=pcoa1)) +
  geom_point(data = data.cond,
             aes(x = ARO.log,
                 y = MDS2,
                 color = pop),
             alpha = 0.5,
             size = 2) +
  
  geom_ribbon(aes(x = ARO.log,
                  y = pcoa2,
                  ymin = lowerCI,
                  ymax = upperCI,
                  fill = pop),
              alpha = 0.3) +
  geom_line(aes(x = ARO.log,
                y = pcoa2,
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
  scale_linetype_manual(name = "Population",
                        labels = c(EMuro = "EMuro",
                                   DMuro = "DMuro",
                                   EPirio = "EPirio"),
                        values = c("DMuro" = 2,
                                   "EMuro" = 2,
                                   "EPirio" = 1)) +
  annotate("text",
           x = 3.3,
           y = -0.7,
           size = 7,
           label = paste("R2=", 
                         round(rpt.2$R$Fixed, 2),
                         sep = ""),
           family = "Times New Roman") +
  ylab("PCoA 2 score") + 
  xlab("Aromatic plant (log(%))") + 
  
  theme_bw(base_size = 20,
           base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 20, face = "bold"),
        #legend.position= c(0.20,0.25), 
        legend.position = "none",
        legend.text = element_text(size=20), 
        legend.key.height = unit(0.3, "cm"), 
        legend.background = element_rect(fill = "transparent"))
p.2  

# create trait arrows
empty <- data.frame(x = c(0,1),
                    y = c(0,1))
empty$lab <- c("handling aggression", "developmental index")

p.2.leg <- ggplot(empty, aes(x=x, y=y))+
  geom_segment(aes(x = 0.25,
                   y = 0,
                   xend = 0.25,
                   yend = 1),
               arrow = arrow(length = unit(0.2, "npc")),
               color = "black",
               alpha = 1,
               linetype = 1,
               size = 0.3) +
  geom_segment(aes(x = 0.75,
                   y = 1,
                   xend = 0.75,
                   yend = 0),
               arrow = arrow(length = unit(0.2, "npc")),
               color = "black",
               alpha = 1,
               linetype = 1,
               size = 0.3) +
  
  geom_text(aes(x = c(0.15, 0.65),
                y = c(0.5,0.5),
                label = lab),
            size=7,
            color = "black",
            angle = 90,
            family = "Times New Roman") +
  xlim(0,1)+
  theme_void(base_family = "Times New Roman")+
  theme(panel.background = element_rect(fill = "white", color = "white"))
  
p.2.leg

# Combine -----------------------------------------------------------------

p.1.2 <- ggarrange(p.1.leg, p.1,
                   ncol = 2,
                   widths = c(4,11))
p.1.2

p.2.2 <- ggarrange(p.2.leg, p.2,
                   ncol = 2,
                   widths = c(4,11))
p.2.2

p.eff <- ggarrange(p.1.2, p.2.2,
                   nrow=2,
                   labels = c("B", "C"))
p.eff

p <- ggarrange(p.pcoa, p.eff,
               ncol = 2,
               labels = c("A", ""),
               widths = c(4,3))
p


# Export ------------------------------------------------------------------

ggsave("figure/Figure_6.png",
       p,
       width = 18,
       height = 12,
       dpi = 400)

