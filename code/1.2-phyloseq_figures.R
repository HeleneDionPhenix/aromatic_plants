# Description -------------------------------------------------------------
# Experimental & observational approach
# Create rarefaction curves and control composition bar plots
# Supplementary Material
# Sequencing: Illumina mi-seq June 2024
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-18

# Libraries --------------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(fantaxtic)
library(scales)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/1-phyloseq_objects.RData")

# change object names
comm.exp <- comm
data.exp <- data
ps.exp <- ps # not rarefied
ps.exp.nest <- ps.nest
data.exp.nest <- data.nest


load("data/obs/1-phyloseq_objects_nest.RData")
load("data/obs/1-phyloseq_objects_eggshell.RData")

set.seed(11)


# SECTION 2 ---------------------------------------------------------------

# Rarefaction figures ---------------------------------------------------------

## Fig. S2.1a. Experimental - Nest and Plant ------------------------------------

#create ggplot2 object for rarecurves
rc.exp <- rarecurve(comm.exp, 
                    step=100,
                    label=FALSE, 
                    xlim=c(0,15000), 
                    tidy = T)
names(rc.exp)[1] <- "sample.names"

#add data
rc.exp <- left_join(rc.exp, data.exp, by = "sample.names")

rc.exp$rareGroup <- "sample"
rc.exp[which(rc.exp$nseq < 14027), 
       "rareGroup"] <- "excluded"
rc.exp[which(rc.exp$type %in% c("ctrl_ext", "ctrl_PCR_neg")), 
       "rareGroup"] <- "ctrl_neg"

p.exp <- ggplot(rc.exp,
                aes(x = Sample,
                    y = Species,
                    group = sample.names)) +
  geom_line(aes(color = rareGroup,
                size = rareGroup)) +
  scale_color_manual(name = "",
                     values = c("ctrl_neg" = "blue",
                                "sample" = "black",
                                "excluded" = "red"),
                     labels = c("ctrl_neg" = "Negative control",
                                "sample" = "Sample",
                                "excluded" = "Excluded")) +
  scale_size_manual(name = "",
                    values = c("ctrl_neg" = 1.0,
                               "sample" = 0.3,
                               "excluded" = 1),
                    labels = c("ctrl_neg" = "Negative control",
                               "sample" = "Sample",
                               "excluded" = "Excluded")) +
  geom_vline(xintercept = 14027,
             linetype = 5) +
  xlim(1,15000) +
  xlab("Number of sequences") +
  ylab("Number of species") +
  ggtitle("Experimental Approach - Nest Samples") +
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))
p.exp

## Fig. S2.1b. Observational - Nest --------------------------------------------

#create ggplot2 object for rarecurves
rc.nest <- rarecurve(comm.nest, 
                    step=100,
                    label=FALSE, 
                    xlim=c(0,8000), 
                    tidy = T)
names(rc.nest)[1] <- "id"

#add data
rc.nest <- left_join(rc.nest, mdata.nest, by = "id")

rc.nest$rareGroup <- "sample"
rc.nest[which(rc.nest$nseq < 7363), 
       "rareGroup"] <- "excluded"
rc.nest[which(rc.nest$categorie == "ctrl"), 
       "rareGroup"] <- "ctrl_neg"
rc.nest[which(rc.nest$id == "ctrl-PCR-pos-N1"),
        "rareGroup"] <- "sample"

rc.nest <- rc.nest[order(rc.nest$nseq),]

p.nest <- ggplot(rc.nest,
                aes(x = Sample,
                    y = Species,
                    group = id)) +
  geom_line(aes(color = rareGroup,
                size = rareGroup)) +
  scale_color_manual(name = "",
                     values = c("ctrl_neg" = "blue",
                                "sample" = "black",
                                "excluded" = "red"),
                     labels = c("ctrl_neg" = "Negative control",
                                "sample" = "Sample",
                                "excluded" = "Excluded")) +
  scale_size_manual(name = "",
                     values = c("ctrl_neg" = 1,
                                "sample" = 0.3,
                                "excluded" = 1),
                     labels = c("ctrl_neg" = "Negative control",
                                "sample" = "Sample",
                                "excluded" = "Excluded")) +
  geom_vline(xintercept = 7363,
             linetype = 5) +
  xlim(1,8000) +
  xlab("Number of sequences") +
  ylab("Number of species") +
  ggtitle("Observational Approach - Nest Samples") +
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, face = "bold"))

p.nest

## Fig. S2.1c. Observational - Eggshell ---------------------------------------

#create ggplot2 object for rarecurves
rc.eggshell <- rarecurve(comm.eggshell, 
                     step=100,
                     label=FALSE, 
                     xlim=c(0,2000), 
                     tidy = T)
names(rc.eggshell)[1] <- "id"

#add data
rc.eggshell <- left_join(rc.eggshell, mdata.eggshell, by = "id")

rc.eggshell$rareGroup <- "sample"
rc.eggshell[which(rc.eggshell$nseq < 1464), 
        "rareGroup"] <- "excluded"
rc.eggshell[which(rc.eggshell$categorie == "ctrl"), 
        "rareGroup"] <- "ctrl_neg"
rc.eggshell[which(rc.eggshell$id == "ctrl-PCR-pos-N1"),
        "rareGroup"] <- "sample"

rc.eggshell <- rc.eggshell[order(rc.eggshell$nseq),]

p.eggshell <- ggplot(rc.eggshell,
                 aes(x = Sample,
                     y = Species,
                     group = id)) +
  geom_line(aes(color = rareGroup,
                size = rareGroup)) +
  scale_color_manual(name = "",
                     values = c("ctrl_neg" = "blue",
                                "sample" = "black",
                                "excluded" = "red"),
                     labels = c("ctrl_neg" = "Negative control",
                                "sample" = "Sample",
                                "excluded" = "Excluded sample")) +
  scale_size_manual(name = "",
                    values = c("ctrl_neg" = 1,
                               "sample" = 0.3,
                               "excluded" = 1),
                    labels = c("ctrl_neg" = "Negative control",
                               "sample" = "Sample",
                               "excluded" = "Excluded sample")) +
  geom_vline(xintercept = 1464,
             linetype = 5) +
  xlim(1,2000) +
  xlab("Number of sequences") +
  ylab("Number of species") +
  ggtitle("Observational Approach - Eggshell Samples") +
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(legend.text = element_text(size=16))

p.eggshell


## Legend ------------------------------------------------------------------


legend <- get_legend(p.eggshell, position = "right")

p.eggshell <- ggplot(rc.eggshell,
                     aes(x = Sample,
                         y = Species,
                         group = id)) +
  geom_line(aes(color = rareGroup,
                size = rareGroup)) +
  scale_color_manual(name = "",
                     values = c("ctrl_neg" = "blue",
                                "sample" = "black",
                                "excluded" = "red"),
                     labels = c("ctrl_neg" = "Negative control",
                                "sample" = "Sample",
                                "excluded" = "Excluded sample")) +
  scale_size_manual(name = "",
                    values = c("ctrl_neg" = 1,
                               "sample" = 0.3,
                               "excluded" = 1),
                    labels = c("ctrl_neg" = "Negative control",
                               "sample" = "Sample",
                               "excluded" = "Excluded sample")) +
  geom_vline(xintercept = 1464,
             linetype = 5) +
  xlim(1,2000) +
  ylim(0, 160) +
  xlab("Number of sequences") +
  ylab("Number of species") +
  ggtitle("Observational Approach - Eggshell Samples") +
  theme_bw(base_size = 12,
           base_family = "Times New Roman") +
  theme(legend.position = "none",
        legend.text = element_text(size=16),
        plot.title = element_text(size = 11, face = "bold"))
p.eggshell

## Fig. S2.1. Composite figure --------------------------------------------------------

p.rc <- ggarrange(p.exp, legend, p.nest, p.eggshell,
                  nrow = 2, ncol = 2,
                  labels = c("a", "", "b", "c"))

p.rc

ggsave(
  "figure/Figure_S2.1.png",
  p.rc,
  width = 18,
  height = 18,
  units = "cm",
  dpi = 300
)


# Extraction and field controls compare to samples -----------------------------

## Fig. S2.2. Experimental approach --------------------------------------------

# select only nest and plant samples and calculate the mean composition
ps.exp.comp <- subset_samples(ps.exp, 
                              sample_data(ps.exp)$type %in% c("ctrl_ext",
                                                              "nest",
                                                              "plant"))
table(sample_data(ps.exp.comp)$type)
# Get the most abundant phyla and the most abundant families within those phyla
top_nested.exp <- nested_top_taxa(ps.exp.comp,
                                  top_tax_level = "Phylum",
                                  nested_tax_level = "Family",
                                  n_top_taxa = 3, 
                                  n_nested_taxa = 3,
                                  by_proportion = FALSE)

pal.exp <- taxon_colours(top_nested.exp$ps_obj, tax_level = "Phylum",
                         palette = c(Actinobacteriota = "#345995",
                                     Firmicutes = "#EAC435", 
                                     Proteobacteria = "#FB4D3D",
                                     Other = "#7C8483"))
scales::show_col(pal.exp)

### Relative ----------------------------------------------------------------

# Plot the relative abundances at two levels.
p.bp.exp.comp1 <- plot_nested_bar(ps_obj = top_nested.exp$ps_obj,
                                 top_level = "Phylum",
                                 nested_level = "Family",
                                 palette = pal.exp,
                                 relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ factor(type,
                      labels = c("ctrl_ext" = "Extraction control",
                                 "nest" = "Nest sample",
                                 "plant" = "Plant sample")),
             nrow = 3, ncol = 1,
             scales = "free",
             switch = "y")+
  ylab('Relative abundance') +
  xlab("")+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  scale_y_continuous(labels = label_percent(suffix = "%"))+
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 16),
        plot.margin = margin(t=0,r=0,l=0,b=0, unit = "pt"))

p.bp.exp.comp1

### Absolute ----------------------------------------------------------------

# Plot the absolute abundances at two levels.
p.bp.exp.comp2 <- plot_nested_bar(ps_obj = top_nested.exp$ps_obj,
                              top_level = "Phylum",
                              nested_level = "Family",
                              palette = pal.exp,
                              relative_abundances = FALSE) +
  coord_flip() +
  facet_wrap(~ factor(type,
                      labels = c("ctrl_ext" = "Extraction control",
                                 "nest" = "Nest sample",
                                 "plant" = "Plant sample")),
             nrow = 3, ncol = 1,
             scales = "free",
             switch = "y")+
  ylab('Number of sequences') +
  xlab("")+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  scale_y_continuous(limits = c(0, 70000),
                     labels = label_number(suffix = "K", scale = 1e-3)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=0,l=0,b=0, unit = "pt"))

p.bp.exp.comp2

### Composite ---------------------------------------------------------------

p.ext.comp <- ggarrange(p.bp.exp.comp1, p.bp.exp.comp2,
                   nrow = 1, ncol =2,
                   widths = c(15,20))
p.ext.comp

ggsave("figure/Figure_S2.2.png",
       p.ext.comp,
       dpi = 300,
       width = 10,
       height = 12)

## Fig. S2.3. Observational - Nest --------------------------------------------

table(mdata.nest$type)
table(mdata.nest$categorie)

# select samples
ps.sub.nest <- subset_samples(ps.nest, sample_data(ps.nest)$type != "PCR_ctrl")
mdata.nest.sub <- mdata.nest[sample_names(ps.sub.nest),]
mdata.nest.sub$type <- as.character(mdata.nest.sub$type)
mdata.nest.sub[which(mdata.nest.sub$categorie == "nid"), "type"] <- "nest"
mdata.nest.sub$type <- as.factor(mdata.nest.sub$type)
sample_data(ps.sub.nest) <- mdata.nest.sub

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.nest <- nested_top_taxa(ps.sub.nest,
                                   top_tax_level = "Phylum",
                                   nested_tax_level = "Family",
                                   n_top_taxa = 3, 
                                   n_nested_taxa = 3,
                                   by_proportion = FALSE)

pal.nest <- taxon_colours(top_nested.nest$ps_obj, tax_level = "Phylum",
                          palette = c(Actinobacteriota = "#345995",
                                      Firmicutes = "#EAC435", 
                                      Proteobacteria = "#FB4D3D",
                                      Other = "#7C8483"))
scales::show_col(pal.nest)

### Relative ----------------------------------------------------------------

# Plot the relative abundances at two levels.
p.bp.nest.1 <- plot_nested_bar(ps_obj = top_nested.nest$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.nest,
                               relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ factor(type,
                      levels = c("field_ctrl",
                                 "ext_ctrl",
                                 "nest"),
                      labels = c("Field control",
                                 "Extraction control",
                                 "Nest sample")),
             nrow = 3, ncol = 1,
             scales = "free",
             switch = "y")+
  ylab('Relative abundance') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = "%"))+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=10,l=0,b=10, unit = "pt"))

p.bp.nest.1

### Absolute ----------------------------------------------------------------

# Plot the absolute abundances at two levels.
p.bp.nest.2 <- plot_nested_bar(ps_obj = top_nested.nest$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.nest,
                               relative_abundances = FALSE) +
  coord_flip() +
  facet_wrap(~type,
             nrow = 3, ncol = 1,
             scales = "free")+
  ylab('Number of sequences') +
  xlab("")+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  scale_y_continuous(limits = c(0, 20000),
                     labels = label_number(suffix = "K", scale = 1e-3)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=0,l=0,b=10, unit = "pt"))

p.bp.nest.2

### Composite ---------------------------------------------------------------

p.ext <- ggarrange(p.bp.nest.1, p.bp.nest.2,
                   nrow = 1, ncol =2,
                   widths = c(15,20))
p.ext

ggsave("figure/Figure_S2.3.png",
       p.ext,
       dpi = 300,
       width = 10,
       height = 10)

## Fig. S2.4. Observational - Eggshell ----------------------------------------

table(mdata.eggshell$type)

ps.sub.egg <- subset_samples(ps.eggshell, 
                             sample_data(ps.eggshell)$type != "PCR_ctrl")

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.egg <- nested_top_taxa(ps.sub.egg,
                                   top_tax_level = "Phylum",
                                   nested_tax_level = "Family",
                                   n_top_taxa = 3, 
                                   n_nested_taxa = 3,
                                   by_proportion = FALSE)

pal.egg <- taxon_colours(top_nested.egg$ps_obj, tax_level = "Phylum",
                          palette = c(Actinobacteriota = "#345995",
                                      Bacteroidota = "purple3", 
                                      Proteobacteria = "#FB4D3D",
                                      Other = "#7C8483"))
scales::show_col(pal.egg)

### Relative ----------------------------------------------------------------

# Plot the relative abundances at two levels.
p.bp.egg.1 <- plot_nested_bar(ps_obj = top_nested.egg$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.egg,
                               relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ factor(type,
                      levels = c("field_ctrl",
                                 "ext_ctrl",
                                 "swab_w"),
                      labels = c("Field control",
                                 "Extraction control",
                                 "Eggshell swab sample")),
             nrow = 3, ncol = 1,
             scales = "free",
             switch = "y")+
  ylab('Relative abundance') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = "%"))+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=10,l=0,b=10, unit = "pt"))

p.bp.egg.1

### Absolute ----------------------------------------------------------------

# Plot the absolute abundances at two levels.
p.bp.egg.2 <- plot_nested_bar(ps_obj = top_nested.egg$ps_obj,
                               top_level = "Phylum",
                               nested_level = "Family",
                               palette = pal.egg,
                               relative_abundances = FALSE) +
  coord_flip() +
  facet_wrap(~type,
             nrow = 3, ncol = 1,
             scales = "free")+
  ylab('Number of sequences') +
  xlab("")+
  theme_classic(base_size=16,
                base_family = "Times New Roman") +
  scale_y_continuous(limits = c(0, 12000),
                     labels = label_number(suffix = "K", scale = 1e-3)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=0,r=0,l=0,b=10, unit = "pt"))

p.bp.egg.2


### Composite ---------------------------------------------------------------

p.egg <- ggarrange(p.bp.egg.1, p.bp.egg.2,
                   nrow = 1, ncol =2,
                   widths = c(15,25))
p.egg

ggsave("figure/Figure_S2.4.png",
       p.egg,
       dpi = 300,
       width = 10,
       height = 10)


## Fig. S2.5. Eggshell composition comparison ---------------------------------

mdata.egg.neg <- filter(mdata.eggshell, type != "PCR_ctrl")
comm.egg.neg <- comm.eggshell[rownames(mdata.egg.neg),]
comm.egg.neg <- comm.egg.neg[rowSums(comm.egg.neg)>=1464,]
mdata.egg.neg <- mdata.egg.neg[rownames(comm.egg.neg),]

dist.egg.neg <- vegdist(comm.egg.neg, method = "hellinger")

disper.egg.neg <- betadisper(dist.egg.neg, mdata.egg.neg$categorie)
anova(disper.egg.neg)
plot(disper.egg.neg, label = FALSE)

mod.egg.neg <- capscale(comm.egg.neg ~ categorie,
                        data = mdata.egg.neg,
                        distance = "hellinger")
anova(mod.egg.neg, by = "term")
RsquareAdj(mod.egg.neg)
anova(mod.egg.neg)

# Figure

perc <- round(100*(summary(mod.egg.neg)$cont$importance[2,1:2]), 2)

#extract sites scores
sites <- scores(mod.egg.neg, 
                choices = c(1,2), 
                display = "sites",
                scaling = 2)
sites <- data.frame(sites)

#create dataset for plotting
data.egg.neg <- mdata.egg.neg[rownames(sites),]
data.egg.neg <- cbind(sites, data.egg.neg)

p <- ggplot(data.egg.neg, aes(x = CAP1, y = MDS1, 
                              fill = categorie)) +
  geom_point(aes(fill = categorie,
                 shape = categorie),
             color = "black",
             size = 2.5)+
  annotate("text",
           x = 2.3,
           y = 1.9,
           label = paste("adjusted R2=0.01", 
                         #round(RsquareAdj(mod.aro)$adj.r.squared, 3),
                         " p=0.02",
                         sep = ""),
           family = "Times New Roman",
           size = 5)+
  scale_shape_manual(name = "",
                     labels = c(swab = "Eggshell sample",
                                ctrl = "Negative control"),
                     values = c(swab = 21,
                                ctrl = 22)) +
  scale_fill_manual(name = "",
                    labels = c(swab = "Eggshell sample",
                               ctrl = "Negative control"),
                    values = c(swab = "black",
                               ctrl = "red")) +
  xlab(paste("RDA 1 (", perc[1], "%)", sep = "")) + 
  ylab(paste("MDS 1 (", perc[2], "%)", sep = "")) + 
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(legend.position = c(0.15, 0.15),
        legend.title = element_text(size=14), 
        legend.text = element_text(size=14), 
        legend.key.height = unit(0.3, "cm"),
        legend.spacing.y = unit(0, 'cm'),
        legend.background = element_rect(fill = "transparent"))+
  guides(shape = guide_legend(override.aes = list(size =3)),
         fill = guide_legend(override.aes = list(shape = 21, size =3)))


p

ggsave("figure/Figure_S2.5.png",
       p,
       width = 9.5,
       height = 7,
       dpi = 300)

# PCR control - bar plots ------------------------------------------------

# Mock community composition 
# https://zymoresearch.eu/products/zymobiomics-microbial-community-dna-standard
# Bacillus subtilis - 12%
# Enterococcus faecalis - 12% 
# Escherichia coli - 12%
# Lactobacillus fermentum - 12%
# Listeria monocytogenes - 12%
# Pseudomonas aeruginosa - 12%
# Salmonella enterica - 12%
# Staphylococcus aureus - 12%

## Fig. S2.6. Experimental --------------------------------------------------

### c. Absolute ----------------------------------------------------------------

table(sample_data(ps.exp)$type)

# subset PCR controls
ps.exp.pcr <- subset_samples(ps.exp, sample_data(ps.exp)$type %in% c("ctrl_PCR_neg",
                                                                     "ctrl_PCR_pos"))
ps.exp.pcr <- prune_taxa(taxa_sums(ps.exp.pcr) > 0, ps.exp.pcr)

# extract legend labels and colors
p.pcr <- plot_bar(ps.exp.pcr, fill="Genus")
g <- ggplot_build(p.pcr)
pal <- data.frame(colours = unique(g$data[[1]]["fill"]), 
                  label = unique(g$plot$data[, g$plot$labels$fill]))
pal <- setNames(pal$fill, pal$label)

# plot number of sequences
p.pcr <- plot_bar(ps.exp.pcr, fill="Genus") +
  scale_x_discrete(name ="",
                   labels=c("N1",
                            "S1")) +
  facet_wrap(~type,
             scales = "free_x",
             labeller = as_labeller(c("ctrl_PCR_neg" = "PCR negative control",
                                      "ctrl_PCR_pos" = "PCR positive control")))+
  ylab("Number of bacterial sequences") +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.pcr

### b. Relative - positive ---------------------------------------------------

ps.exp.pcr.pos <- subset_samples(ps.exp.pcr, 
                                 sample_data(ps.exp.pcr)$type == "ctrl_PCR_pos")

# Transform in relative abundance
pos.taxo <- names(sort(taxa_sums(ps.exp.pcr.pos), decreasing=TRUE))
pos.taxo.ps <- transform_sample_counts(ps.exp.pcr.pos, function(OTU) OTU/sum(OTU))
pos.taxo.ps<- prune_taxa(pos.taxo, pos.taxo.ps)

# Extract legend elements
p.pos <- plot_bar(pos.taxo.ps, fill="Genus")
g <- ggplot_build(p.pos)

# Change colors to fit with first plot
pal.pos <- pal
pal.pos <- pal.pos[unique(g$plot$data[, g$plot$labels$fill])]
pal.pos[which(is.na(pal.pos))] <- "grey50"

# Plot relative abundances
p.pos <- plot_bar(pos.taxo.ps, fill="Genus") +
  scale_fill_manual(values = pal.pos,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("ctrl_PCR_pos" = "PCR positive control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.pos

### a. Relative - negative ----------------------------------------------------

ps.exp.pcr.neg <- subset_samples(ps.exp.pcr, 
                                 sample_data(ps.exp.pcr)$type == "ctrl_PCR_neg")

# Transform in relative abundance
neg.taxo <- names(sort(taxa_sums(ps.exp.pcr.neg), decreasing=TRUE))
neg.taxo.ps <- transform_sample_counts(ps.exp.pcr.neg, function(OTU) OTU/sum(OTU))
neg.taxo.ps<- prune_taxa(neg.taxo, neg.taxo.ps)

# Extract legend elements
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") 
g <- ggplot_build(p.neg)

# Change colors to fit with first plot
pal.neg <- pal
pal.neg <- pal.neg[unique(g$plot$data[, g$plot$labels$fill])]
pal.neg[which(is.na(pal.neg))] <- "grey50"

# Plot relative abundances
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") +
  scale_fill_manual(values = pal.neg,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("ctrl_PCR_neg" = "PCR negative control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.neg

### Composite ---------------------------------------------------------------

p.pcr.all.exp <- ggarrange(p.neg, p.pos, p.pcr,
                           nrow = 1, ncol = 3,
                           widths = c(2,2,6),
                           labels = c("a", "b", "c"))

## Fig. S2.6. Observational --------------------------------------------------

### f. Absolute ----------------------------------------------------------------

# subset PCR controls
ps.nest.pcr <- subset_samples(ps.nest, sample_data(ps.nest)$type == "PCR_ctrl")
ps.nest.pcr <- prune_taxa(taxa_sums(ps.nest.pcr) > 0, ps.nest.pcr)

ps.eggshell.pcr <- subset_samples(ps.eggshell, sample_data(ps.eggshell)$type == "PCR_ctrl")
ps.eggshell.pcr <- prune_taxa(taxa_sums(ps.eggshell.pcr) > 0, ps.eggshell.pcr)

# merge pcr controls
ps.pcr <- merge_phyloseq(ps.nest.pcr, ps.eggshell.pcr)

#add positive and negative
mdata.pcr <- sample_data(ps.pcr)
mdata.pcr$pcr <- "positive"
mdata.pcr[which(mdata.pcr$id %in% c("ctrl-PCR-neg-N1",
                                    "ctrl-PCR-neg-S1",
                                    "ctrl-PCR-neg-S2")), "pcr"] <- "negative"
sample_data(ps.pcr) <- mdata.pcr

# extract legend labels and colors
p.pcr <- plot_bar(ps.pcr, fill="Genus")
g <- ggplot_build(p.pcr)
pal <- data.frame(colours = unique(g$data[[1]]["fill"]), 
                  label = unique(g$plot$data[, g$plot$labels$fill]))
pal <- setNames(pal$fill, pal$label)

# plot number of sequences
p.pcr <- plot_bar(ps.pcr, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("N1",
                            "S1",
                            "S2")) +
  facet_wrap(~pcr,
             scales = "free_x",
             labeller = as_labeller(c("negative" = "PCR negative control",
                                      "positive" = "PCR positive control")))+
  ylab("Number of bacterial sequences") +
  scale_fill_manual(values = pal) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.pcr

### e. Relative - positive ---------------------------------------------------

ps.pcr.pos <- subset_samples(ps.pcr, sample_data(ps.pcr)$pcr == "positive")

# Transform in relative abundance
pos.taxo <- names(sort(taxa_sums(ps.pcr.pos), decreasing=TRUE))
pos.taxo.ps <- transform_sample_counts(ps.pcr.pos, function(OTU) OTU/sum(OTU))
pos.taxo.ps<- prune_taxa(pos.taxo, pos.taxo.ps)

# Extract legend elements
p.pos <- plot_bar(pos.taxo.ps, fill="Genus")
g <- ggplot_build(p.pos)

# Change colors to fit with first plot
pal.pos <- pal
pal.pos <- pal.pos[unique(g$plot$data[, g$plot$labels$fill])]
pal.pos[which(is.na(pal.pos))] <- "grey50"

# Plot relative abundances
p.pos <- plot_bar(pos.taxo.ps, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("N1",
                            "S1",
                            "S2")) +
  scale_fill_manual(values = pal.pos,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("PCR_ctrl" = "PCR positive control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.pos

### d. Relative - negative ----------------------------------------------------

ps.pcr.neg <- subset_samples(ps.pcr, sample_data(ps.pcr)$pcr == "negative")

# Transform in relative abundance
neg.taxo <- names(sort(taxa_sums(ps.pcr.neg), decreasing=TRUE))
neg.taxo.ps <- transform_sample_counts(ps.pcr.neg, function(OTU) OTU/sum(OTU))
neg.taxo.ps<- prune_taxa(neg.taxo, neg.taxo.ps)

# Extract legend elements
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") 
g <- ggplot_build(p.neg)

# Change colors to fit with first plot
pal.neg <- pal
pal.neg <- pal.neg[unique(g$plot$data[, g$plot$labels$fill])]
pal.neg[which(is.na(pal.neg))] <- "grey50"

# Plot relative abundances
p.neg <- plot_bar(neg.taxo.ps, fill="Genus") +
  scale_x_discrete(name ="", 
                   labels=c("N1",
                            "S1",
                            "S2")) +
  scale_fill_manual(values = pal.neg,
                    guide = "none") +
  facet_wrap(~type,
             labeller = as_labeller(c("PCR_ctrl" = "PCR negative control")))+
  ylab("Bacterial relative abundance") +
  theme_bw(base_size = 16,
           base_family = "Times New Roman")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p.neg

### Composite ---------------------------------------------------------------

p.pcr.all.obs <- ggarrange(p.neg, p.pos, p.pcr,
                           nrow = 1, ncol = 3,
                           widths = c(2,2,6),
                           labels = c("d", "e", "f"))

p.pcr.all.obs

## Fig. S2.6. Composite figure --------------------------------------

# Create treatment labels
empty <- data.frame(x = 0,
                    y = 0.5)

empty$lab <- c("Experimental approach")

lab.exp <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(label = lab),
            size=7,
            color = "black",
            hjust = 0,
            family = "Times New Roman") +
  xlim(0, 15) +
  ylim(0,1) +
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))
lab.exp

empty$lab <- c("Observational approach")

lab.obs <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(label = lab),
            size=7,
            color = "black",
            hjust = 0,
            family = "Times New Roman") +
  xlim(0, 15) +
  ylim(0,1) +
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))
lab.obs

p.pcr.all <- ggarrange(lab.exp,
                       p.pcr.all.exp, 
                       lab.obs,
                       p.pcr.all.obs,
                       nrow = 4, ncol = 1,
                       heights = c(1,10,1,10))
p.pcr.all

ggsave("figure/Figure_S2.6.png",
       p.pcr.all,
       dpi = 300,
       width = 15,
       height = 10)


# SECTION 3 ---------------------------------------------------------------

# Composition - comparison between groups ---------------------------------

## Fig. S3.1. Experimental approach ----------------------------------------

# Combine nest samples and plant samples (rarefied)
ps.exp.sub <- merge_phyloseq(ps.exp.nest, ps.plant)

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.exp <- nested_top_taxa(ps.exp.sub,
                                  top_tax_level = "Phylum",
                                  nested_tax_level = "Family",
                                  n_top_taxa = 3, 
                                  n_nested_taxa = 3,
                                  by_proportion = FALSE)

pal.exp <- taxon_colours(top_nested.exp$ps_obj, tax_level = "Phylum",
                         palette = c(Actinobacteriota = "#345995",
                                     Firmicutes = "#EAC435", 
                                     Proteobacteria = "#FB4D3D",
                                     Other = "#7C8483"))
scales::show_col(pal.exp)


# Plot the relative abundances by treatment for nests and plants
bp.sample.exp.1 <- plot_nested_bar(ps_obj = top_nested.exp$ps_obj,
                              top_level = "Phylum",
                              nested_level = "Family",
                              palette = pal.exp,
                              relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ factor(stage_pheno,
                      levels = c("prebrood", "brood", "breed")) + 
               type + 
               factor(traitement,
                      levels = c("ach", "pul", "men", "lav", "imm",
                                 "lavmen", "lavimm", "menimm", "ctr")),
             ncol = 4, nrow = 9,
             scales = "free",
             dir = "v",
             switch = "y")+
  ylab('Relative abundance (%)') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = ""))+
  theme_classic(base_size=14,
                base_family = "Times New Roman") +
  ggtitle("Nest in pre-incubation    Nest in incubation          Nest in post-hatching      Plant sample")+
  theme(plot.title = element_text(size = 14, hjust = 0),
        legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=10,r=10,l=0,b=10, unit = "pt"))

bp.sample.exp.1 

# Create treatment labels
empty <- data.frame(x = rep(2.5, 9),
                    y = c(0.5, 1.5, 2.5, 
                          3.5, 4.5, 5.5, 
                          6.5, 7.5, 8.5))

empty$lab <- c("Control",
               "Mentha - Helichryum",
               "Lavandula - Helichryum",
               "Lavandula - Mentha",
               "Helichryum", "Lavandula", "Mentha", "Pulicaria", "Achillea")

p.lab <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(x = Inf,
                label = lab),
            size=5,
            color = "black",
            hjust = 1,
            family = "Times New Roman") +
  xlim(0, 5) +
  ylim(-0.1,8.9) +
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))
p.lab

# Combine treatment and plot
bp.sample.exp <- ggarrange(p.lab, bp.sample.exp.1 ,
                    nrow = 1, ncol =2,
                    widths = c(4,20))
bp.sample.exp

ggsave(
  "figure/Figure_S3.1.png",
  bp.sample.exp,
  width = 32,
  height = 18,
  units = "cm",
  dpi = 500
)


## Fig. S3.3a. Observational - Nest ----------------------------------------

# load metadata
load("data/obs/0-data.RData")
rm(ad, demo, pou, compo.nest, mdata.e, mdata.n, compo)

# add data to phyloseq object
data.nest.r <- filter(compo.stage, id %in% sample_names(ps.nest.r))
rownames(data.nest.r) <- data.nest.r$id

# To order by quantity of aromatic plant
data.nest.r <- data.nest.r[order(data.nest.r$ARO),]
data.nest.r$ordre <- 1:nrow(data.nest.r)

data.nest.r <- data.nest.r[sample_names(ps.nest.r),]
sample_data(ps.nest.r) <- data.nest.r
sample_names(ps.nest.r) <- sample_data(ps.nest.r)$ordre

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.obs.nest <- nested_top_taxa(ps.nest.r,
                                  top_tax_level = "Phylum",
                                  nested_tax_level = "Family",
                                  n_top_taxa = 3, 
                                  n_nested_taxa = 3,
                                  by_proportion = FALSE)

pal.obs.nest <- taxon_colours(top_nested.obs.nest$ps_obj,
                              tax_level = "Phylum",
                         palette = c(Actinobacteriota = "#345995",
                                     Firmicutes = "#EAC435", 
                                     Proteobacteria = "#FB4D3D",
                                     Other = "#7C8483"))
scales::show_col(pal.obs.nest)

# Plot the relative abundances by treatment for nests and plants
bp.obs.nest.1 <- plot_nested_bar(ps_obj = top_nested.obs.nest$ps_obj,
                                   top_level = "Phylum",
                                   nested_level = "Family",
                                   palette = pal.obs.nest,
                                   relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ stade + pop,
             nrow = 3, ncol = 2,
             scales = "free",
             dir = "v",
             switch = "y")+
  ylab('Relative abundance') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = "%"))+
  theme_classic(base_size=14,
                base_family = "Times New Roman") +
  ggtitle("Nest in incubation                        Nest in post-hatching")+
  theme(plot.title = element_text(size = 14, hjust = 0),
        legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=10,r=10,l=0,b=10, unit = "pt"))

bp.obs.nest.1

# Create treatment labels
empty <- data.frame(x = rep(2.5, 3),
                    y = c(0.5, 1.5, 2.5))

empty$lab <- c("EPirio",
               "EMuro",
               "DMuro")

p.lab <- ggplot(empty, aes(x=x, y=y))+
  
  geom_text(aes(label = lab),
            size=5,
            color = "black",
            angle = 90,
            family = "Times New Roman") +
  xlim(0, 3.5) +
  ylim(-0.1,2.9) +
  theme_void()+
  theme(panel.background = element_rect(fill = "white", color = "white"))
p.lab

# Combine treatment and plot
bp.obs.nest <- ggarrange(p.lab, bp.obs.nest.1 ,
                           nrow = 1, ncol =2,
                           widths = c(1,30))
bp.obs.nest

## Fig. S3.3b. Observational - Eggshell ---------------------------------------

# add data to phyloseq object
data.eggshell.r <- filter(compo.stage, id %in% sample_names(ps.eggshell.r))
rownames(data.eggshell.r) <- data.eggshell.r$id

# To order by quantity of aromatic plant
data.eggshell.r <- data.eggshell.r[order(data.eggshell.r$ARO),]
data.eggshell.r$ordre <- 1:nrow(data.eggshell.r)

data.eggshell.r <- data.eggshell.r[sample_names(ps.eggshell.r),]
sample_data(ps.eggshell.r) <- data.eggshell.r
sample_names(ps.eggshell.r) <- sample_data(ps.eggshell.r)$ordre

# Get the most abundant phyla and the most abundant families within those phyla
top_nested.obs.egg <- nested_top_taxa(ps.eggshell.r,
                                       top_tax_level = "Phylum",
                                       nested_tax_level = "Family",
                                       n_top_taxa = 3, 
                                       n_nested_taxa = 3,
                                       by_proportion = FALSE)

pal.obs.egg <- taxon_colours(top_nested.obs.egg$ps_obj,
                              tax_level = "Phylum",
                              palette = c(Actinobacteriota = "#345995",
                                          Bacteroidota = "purple3", 
                                          Proteobacteria = "#FB4D3D",
                                          Other = "#7C8483"))
scales::show_col(pal.obs.egg)

# Plot the relative abundances by treatment for nests and plants
bp.obs.egg.1 <- plot_nested_bar(ps_obj = top_nested.obs.egg$ps_obj,
                                 top_level = "Phylum",
                                 nested_level = "Family",
                                 palette = pal.obs.egg,
                                 relative_abundances = TRUE) +
  coord_flip() +
  facet_wrap(~ pop,
             nrow = 3, ncol = 1,
             scales = "free",
             dir = "v",
             switch = "y")+
  ylab('Relative abundance') +
  xlab("")+
  scale_y_continuous(labels = label_percent(suffix = "%"))+
  theme_classic(base_size=14,
                base_family = "Times New Roman") +
  ggtitle("Eggshell in incubation")+
  theme(plot.title = element_text(size = 14, hjust = 0),
        legend.position = "right",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = margin(t=10,r=10,l=0,b=10, unit = "pt"))

bp.obs.egg.1

# Combine treatment and plot
bp.obs.egg <- ggarrange(p.lab, bp.obs.egg.1 ,
                         nrow = 1, ncol =2,
                         widths = c(1,10))
bp.obs.egg

bp.obs <- ggarrange(bp.obs.nest, bp.obs.egg,
                    nrow = 1, ncol =2,
                    widths = c(30,20),
                    labels = "auto")
bp.obs

ggsave(
  "figure/Figure_S3.3.png",
  bp.obs,
  width = 35,
  height = 18,
  units = "cm",
  dpi = 500
)
