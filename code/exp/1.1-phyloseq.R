# Description -------------------------------------------------------------
# Experimental approach
# Effect of aromatic plants on nest bacterial microbiota
# Rarefaction and creation of phyloseq objects
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-03-05

# Libraries --------------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(picante); packageVersion("picante")
library(dplyr)

# Import data -------------------------------------------------------------

rm(list = ls())

load("data/exp/bioinfo/0-dada2.RData")
load("data/exp/0-data.RData")

# Phyloseq ----------------------------------------------------------------

#### create community, taxonomy, metadata files and match them by sample name
comm <- seqtab.nochim
#### load metadata
all(rownames(comm) %in% rownames(data))

# create phyloseq objects --------------------------------------------------

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               sample_data(data))


# Remove some taxa --------------------------------------------------------

#How many sequences are non-bacterial
(chloro = subset_taxa(ps, Order=="Chloroplast"))#3
sum(chloro@otu_table@.Data)#(237 sequences)
#(euk = subset_taxa(ps, Kingdom=="Eukaryota"))#0
(arch = subset_taxa(ps, Kingdom=="Archaea"))#1
sum(arch@otu_table@.Data)#6
(mito = subset_taxa(ps, Family=='Mitochondria'))#26
sum(mito@otu_table@.Data)#5012
sum(comm)
(237+6+5012)/sum(comm)*100

#How many sequences are not identified
(na.o = subset_taxa(ps, is.na(Order)))#753
(na.c = subset_taxa(ps, is.na(Class)))#396
(na.p = subset_taxa(ps, is.na(Phylum)))#259
sum(na.p@otu_table@.Data)#22002
#(na.k = subset_taxa(ps, is.na(Kingdom)))#0
sum(na.p@otu_table@.Data)/sum(comm)*100

#Remove eucaryote, archea, Chloroplast, mitochondria and NA at phylum level
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(chloro))))
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(euk))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(arch))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(mito))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(na.p))))

rm(list='chloro','arch','mito','na.o',"na.c","na.p")

#Extract picante/vegan format objects
comm <- otu_table(ps)
comm <- comm@.Data
taxo <- tax_table(ps)
metadata<- sample_data(ps)

all(rownames(comm) %in% rownames(data))

data <- data[rownames(comm),]

data$nseq <- rowSums(comm)
data$rich <- specnumber(comm)

# Rarefaction ------------------------------------------------------------

rarecurve(comm, step=100, label=FALSE)
rarecurve(comm, step=100, label=TRUE, xlim = c(0, 10000))

data[order(data$nseq),c("sample.names","type","nseq")]
#14027 sequences for the first sample

## Rarefaction threshold ---------------------------------------------------

#Comparison of rarefaction threshold
S <- specnumber(comm)

Srare <- rarefy(comm, 14027)
p1<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "14027 (0 sample, all negative ctrls)")
abline(0, 1)


rarecurve(comm, step=100,label=FALSE, xlim=c(0,15000))
abline(v = 14027)

# Rarefaction -------------------------------------------------------------

set.seed(1)
comm.r <- rrarefy(comm[which(apply(comm,1,sum)>=14027),], sample=14027)
comm.r <- comm.r[,apply(comm.r,2,sum)>0]
taxo.r <- taxo[colnames(comm.r),]
data.r <- data[rownames(comm.r),]
ps.r <- phyloseq(otu_table(comm.r, taxa_are_rows=FALSE),
                 tax_table(taxo.r),
                 sample_data(data.r))
rarecurve(comm.r, step=100,label=TRUE)


# Create phyloseq objects -------------------------------------------------

table(data.r$type)
#still 2 ctrls (2 PCR positive ctrl)

#Remove controls
ps.r <- subset_samples(ps.r, sample_data(ps.r)$category != "control")
comm.r <-otu_table(ps.r)
comm.r <-  comm.r@.Data
taxo.r <- tax_table(ps.r)
data.r <- data.r[rownames(comm.r),]

#Create obejct with only nest samples
ps.nest <- subset_samples(ps.r, sample_data(ps.r)$type == "nest")
comm.nest <-otu_table(ps.nest)
comm.nest <-  comm.nest@.Data
taxo.nest <- tax_table(ps.nest)
data.nest <- data.r[rownames(comm.nest),]

#Create obejct with only plant samples
ps.plant <- subset_samples(ps.r, sample_data(ps.r)$type == "plant")
comm.plant <-otu_table(ps.plant)
comm.plant <-  comm.plant@.Data
taxo.plant <- tax_table(ps.plant)
data.plant <- data.r[rownames(comm.plant),]

# Clean and save environment ----------------------------------------------

rm(metadata, seqtab.nochim, taxa, S, p1,
   sample.names, Srare)

save.image("data/exp/1-phyloseq_objects.RData")

