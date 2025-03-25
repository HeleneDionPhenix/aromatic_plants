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


# Figure S4.1. Developmental index ~ hatching date ------------------------

newdat <- expand.grid(date_eclo.s = seq(min(data.cond$date_eclo.s, na.rm =T),
                                        max(data.cond$date_eclo.s, na.rm = T),
                                        length = 100),
                      pop = levels(data.cond$pop),
                      nichoir = 1)

pred <- predict(mod.devel,
                newdat,
                se.fit = T,
                allow.new.levels = T,
                re.form = NA)

newdat$devel <- pred$fit
newdat$lowci <- pred$fit - 1.96*pred$se.fit
newdat$upperci <- pred$fit + 1.96*pred$se.fit

rp <- range(data.cond[which(data.cond$pop == "EPirio"), "date_eclo.s"], na.rm=T)
rm <- range(data.cond[which(data.cond$pop == "DMuro"), "date_eclo.s"], na.rm=T)
re <- range(data.cond[which(data.cond$pop == "EMuro"), "date_eclo.s"], na.rm=T)

# Truncate and unscale by pop
newdat.EP <- filter(newdat, pop == "EPirio")
newdat.EP <- filter(newdat.EP, !(date_eclo.s < min(rp) | date_eclo.s > max(rp)))
data.cond.EP <- filter(data.cond, pop == "EPirio") 
newdat.EP$date_eclo <- newdat.EP$date_eclo.s * sd(data.cond.EP$date_eclo, na.rm=T) + mean(data.cond.EP$date_eclo, na.rm=T)

newdat.DM <- filter(newdat, pop == "DMuro")
newdat.DM <- filter(newdat.DM, !(date_eclo.s < min(rm) | date_eclo.s > max(rm)))
data.cond.DM <- filter(data.cond, pop == "DMuro") 
newdat.DM$date_eclo <- newdat.DM$date_eclo.s * sd(data.cond.DM$date_eclo, na.rm=T) + mean(data.cond.DM$date_eclo, na.rm=T)

newdat.EM <- filter(newdat, pop == "EMuro")
newdat.EM <- filter(newdat.EM, !(date_eclo.s < min(re) | date_eclo.s > max(re)))
data.cond.EM <- filter(data.cond, pop == "EMuro") 
newdat.EM$date_eclo <- newdat.EM$date_eclo.s * sd(data.cond.EM$date_eclo, na.rm=T) + mean(data.cond.EM$date_eclo, na.rm=T)

newdat <- rbind(newdat.EP, newdat.DM, newdat.EM)

newdat$date <- lubridate::as_date(newdat$date_eclo)
data.cond$date <- lubridate::as_date(data.cond$date_eclo)

p.devel <- ggplot(data = newdat,
                  aes(x = date,
                      y = devel,
                      group = pop)) +
  geom_point(data = data.cond,
             aes(color = pop),
             alpha = 0.5) +
  geom_ribbon(aes(ymin = lowci,
                  ymax = upperci,
                  fill = pop),
              alpha = 0.5)+
  geom_line(aes(color = pop)) +
  annotate("text",
           x = lubridate::as_date(155),
           y = -2.5,
           label = paste("R2=",
                         round(rpt.devel$R$Fixed, 2),
                         sep = ""),
           family = "Times New Roman")+
  
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
  xlab("Hatching date") + 
  ylab("Developmental index")+
  theme_bw(base_size = 16,
           base_family = "Times New Roman")

p.devel

# difference between pop
newdat <- expand.grid(date_eclo.s = mean(data.cond$date_eclo.s, na.rm =T),
                      pop = levels(data.cond$pop),
                      nichoir = 1)

pred <- predict(mod.devel,
                newdat,
                se.fit = T,
                allow.new.levels = T,
                re.form = NA)

newdat$devel <- pred$fit
newdat$lowci <- pred$fit - 1.96*pred$se.fit
newdat$upperci <- pred$fit + 1.96*pred$se.fit

p.pop <- ggplot(data = newdat,
                aes(x = pop,
                    y = devel)) +
  geom_jitter(data = data.cond,
              aes(color = pop),
              alpha = 0.5,
              width = 0.15) +
  geom_point()+
  geom_linerange(aes(ymin = lowci,
                     ymax = upperci)) +
  scale_color_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  xlab("") + 
  ylab("Developmental index")+
  theme_bw(base_size = 16,
           base_family = "Times New Roman") +
  theme(legend.position = "none")


p.pop

p1 <- ggarrange(p.pop, p.devel,
                widths = c(8,15),
                labels = c("A", "B"))
p1

ggsave("figure/Figure_S5.1.png",
       p1,
       width = 23,
       height = 10,
       units = "cm",
       dpi = 300)

# Figure S4.2. Aromatic plant quantity ~ hatching date ------------------------

newdat <- expand.grid(date_eclo.s = seq(min(data.nest$date_eclo.s, na.rm =T),
                                        max(data.nest$date_eclo.s, na.rm = T),
                                        length = 100),
                      pop = levels(data.nest$pop))

pred <- predict(mod.aro,
                newdat,
                se.fit = T)

newdat$ARO.log <- pred$fit
newdat$lowci <- pred$fit - 1.96*pred$se.fit
newdat$upperci <- pred$fit + 1.96*pred$se.fit

rp <- range(data.nest[which(data.nest$pop == "EPirio"), "date_eclo.s"], na.rm=T)
rm <- range(data.nest[which(data.nest$pop == "DMuro"), "date_eclo.s"], na.rm=T)
re <- range(data.nest[which(data.nest$pop == "EMuro"), "date_eclo.s"], na.rm=T)


# Truncate and unscale by pop
newdat.EP <- filter(newdat, pop == "EPirio")
newdat.EP <- filter(newdat.EP, !(date_eclo.s < min(rp) | date_eclo.s > max(rp)))
data.nest.EP <- filter(data.nest, pop == "EPirio") 
newdat.EP$date_eclo <- newdat.EP$date_eclo.s * sd(data.nest.EP$date_eclo, na.rm=T) + mean(data.nest.EP$date_eclo, na.rm=T)

newdat.DM <- filter(newdat, pop == "DMuro")
newdat.DM <- filter(newdat.DM, !(date_eclo.s < min(rm) | date_eclo.s > max(rm)))
data.nest.DM <- filter(data.nest, pop == "DMuro") 
newdat.DM$date_eclo <- newdat.DM$date_eclo.s * sd(data.nest.DM$date_eclo, na.rm=T) + mean(data.nest.DM$date_eclo, na.rm=T)


newdat.EM <- filter(newdat, pop == "EMuro")
newdat.EM <- filter(newdat.EM, !(date_eclo.s < min(re) | date_eclo.s > max(re)))
data.nest.EM <- filter(data.nest, pop == "EMuro") 
newdat.EM$date_eclo <- newdat.EM$date_eclo.s * sd(data.nest.EM$date_eclo, na.rm=T) + mean(data.nest.EM$date_eclo, na.rm=T)

newdat <- rbind(newdat.EP, newdat.DM, newdat.EM)

newdat$date <- lubridate::as_date(newdat$date_eclo)
data.nest$date <- lubridate::as_date(data.nest$date_eclo)


p.aro <- ggplot(data = newdat,
                aes(x = date,
                    y = ARO.log,
                    group = pop)) +
  geom_point(data = data.nest,
             aes(color = pop),
             alpha = 0.5) +
  geom_ribbon(aes(ymin = lowci,
                  ymax = upperci,
                  fill = pop),
              alpha = 0.5)+
  geom_line(aes(color = pop)) +
  annotate("text",
           x = lubridate::as_date(155),
           y = 3.5,
           label = paste("R2=",
                         round(summary(mod.aro)$adj.r.squared, 2),
                         sep = ""),
           family = "Times New Roman")+
  
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
  xlab("Hatching date") + 
  ylab("Aromatic plant quantity (log(%))")+
  theme_bw(base_size = 14,
           base_family = "Times New Roman")

p.aro

ggsave("figure/Figure_S5.2.png",
       p.aro,
       width = 15,
       height = 10,
       units = "cm",
       dpi = 300)


# Figure S4.3. Fledging success over the years ----------------------------

# by year

newdat <- expand.grid(pop = levels(data.fledg$pop),
                      an = levels(data.fledg$an))

pred <- predict(mod.fledg,
                newdat,
                se.fit = TRUE,
                type = "link")

inv_link <- family(mod.fledg)$linkinv

newdat$fledg <- inv_link(pred$fit)
newdat$lowerCI <- inv_link(pred$fit - 1.96*pred$se.fit)
newdat$upperCI <- inv_link(pred$fit + 1.96*pred$se.fit)

y2023 <- filter(newdat, an == "2023")
newdat <- filter(newdat, an %in% as.character(2003:2023))

fledg.tab <- newdat

tmp <- select(fledg.tab,
              pop, an)
tmp <- left_join(tmp, mean.fledg.tab, by = "pop")

p.fledg <- ggplot(newdat,
                  aes(y = fledg,
                      x = an,
                      group = pop)) +
  geom_line(data = tmp,
            aes(color = pop),
            linetype = 5,
            alpha = 0.7) +
  geom_line(aes(color = pop),
            position = position_dodge(width = 0.5))+
  geom_point(aes(color = pop),
              position = position_dodge(width = 0.5))+
  geom_linerange(aes(ymin = lowerCI,
                     ymax = upperCI,
                     color = pop),
             position = position_dodge(width = 0.5))+
  geom_point(data = y2023,
              position = position_dodge(width = 0.5),
              color = "red")+

  scale_color_manual(name = "Population",
                     labels = c(EMuro = "EMuro",
                                DMuro = "DMuro",
                                EPirio = "EPirio"),
                     values = c("DMuro" = "steelblue1",
                                "EMuro" = "green3",
                                "EPirio" = "forestgreen"))+
  xlab("Year")+
  ylab("Fledging success (%)")+
  theme_pubclean(base_size = 16,
                 base_family = "Times New Roman")+
  theme(legend.position = "none")


p.fledg  

ggsave(
  "figure/Figure_S5.3.png",
  p.fledg,
  width = 25,
  height = 10,
  units = "cm",
  dpi = 300
)

