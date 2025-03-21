# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on eggshell bacterial microbiota
# Create rarefied phyloseq objects
# Sequencing: Illumina mi-seq June 2024
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-02-28

# Libraries --------------------------------------------------------------

library(phyloseq); packageVersion("phyloseq")
library(picante); packageVersion("picante")
library(dplyr)
library(lubridate)

# Import data -------------------------------------------------------------

rm(list = ls())

# ASV table
load("data/obs/bioinfo/0-dada2_eggshell.RData")

#metadata
mdata <- read.csv2("data/obs/0-data_eggshell_sample.csv",
                   row.names = 1)


# Format variables --------------------------------------------------------

mdata$categorie <- as.factor(mdata$categorie)
mdata$nichoir <- as.factor(mdata$nichoir)
mdata$date <- ymd(mdata$date)
mdata$type <- as.factor(mdata$type)

# Phyloseq ----------------------------------------------------------------

#### create community, taxonomy, metadata files and match them by sample name
comm <- seqtab.nochim
taxo <- taxa
#### load metadata
all(rownames(comm) %in% rownames(mdata))
mdata <- mdata[rownames(comm),]

# create phyloseq objects --------------------------------------------------

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),
               sample_data(mdata))

sample_names(ps) <- sample_data(ps)$id
rownames(mdata) <- mdata$id

# Remove some taxa --------------------------------------------------------

#How many sequences are non-bacterial
#(chloro = subset_taxa(ps, Order=="Chloroplast"))#0
#sum(chloro@otu_table@.Data)
#(euk = subset_taxa(ps, Kingdom=="Eukaryota"))#0
#(arch = subset_taxa(ps, Kingdom=="Archaea"))#0
(mito = subset_taxa(ps, Family=='Mitochondria'))#3
sum(mito@otu_table@.Data)/sum(comm)*100
#total: 0.04%

#How many sequences are not identified
(na.o = subset_taxa(ps, is.na(Order)))#100
(na.c = subset_taxa(ps, is.na(Class)))#30
(na.p = subset_taxa(ps, is.na(Phylum)))#15
sum(na.p@otu_table@.Data)#5313
#(na.k = subset_taxa(ps, is.na(Kingdom)))#0
sum(na.p@otu_table@.Data)/sum(comm)*100
#0.5%

#Remove eucaryote, archea, Chloroplast, mitochondria and NA at phylum level
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(chloro))))
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(euk))))
#ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(arch))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(mito))))
ps = subset_taxa(ps, !(rownames(tax_table(ps)) %in% rownames(tax_table(na.p))))

rm(list='mito','na.o',"na.c","na.p")

#Extract picante/vegan format objects
comm <- otu_table(ps)
comm <- comm@.Data
taxo <- tax_table(ps)
metadata<- sample_data(ps)

mdata$nseq <- rowSums(comm)
mdata$rich <- specnumber(comm)

hist(mdata$nseq, breaks = 20)
hist(mdata$rich)

# Rarefaction ------------------------------------------------------------

rarecurve(comm, step=100, label=FALSE)
rarecurve(comm, step=100, label=TRUE, xlim = c(0, 10000))

mdata[order(mdata$nseq),c("nom_cermo","nseq")]
hist(mdata$nseq, breaks = 50)

## Rarefaction threshold ---------------------------------------------------

#Comparison of rarefaction threshold
S <- specnumber(comm)

Srare <- rarefy(comm, 934)
p1<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "934 (0 samples, 3/7 ctrl)")
abline(0, 1)

Srare <- rarefy(comm, 1464)
p2<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "1464 (3 samples, 4/7 ctrl)")
abline(0, 1)

Srare <- rarefy(comm, 1709)
p3<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "1709 (4 samples, 5/7 ctrl)")
abline(0, 1)

Srare <- rarefy(comm, 3436)
p4<- plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main= "3436 (6 samples, 7/7 ctrl)")
abline(0, 1)

# Rarefaction -------------------------------------------------------------

set.seed(1)
comm.r <- rrarefy(comm[which(apply(comm,1,sum)>=1464),], sample=1464)
comm.r <- comm.r[,apply(comm.r,2,sum)>0]
taxo.r <- taxo[colnames(comm.r),]
mdata.r <- mdata[rownames(comm.r),]
ps.r <- phyloseq(otu_table(comm.r, taxa_are_rows=FALSE),
                 tax_table(taxo.r),
                 sample_data(mdata.r))
rarecurve(comm.r, step=100,label=TRUE)


# Create phyloseq objects -------------------------------------------------

#Remove controls
ps.r <- subset_samples(ps.r, sample_data(ps.r)$categorie != "ctrl")
comm.r <-otu_table(ps.r)
comm.r <-  comm.r@.Data
taxo.r <- tax_table(ps.r)
mdata.r <- mdata.r[rownames(comm.r),]

# Clean and save environment ----------------------------------------------

# rename objects
comm.eggshell <- comm
comm.eggshell.r <- comm.r
mdata.eggshell <- mdata
mdata.eggshell.r <- mdata.r
ps.eggshell <- ps
ps.eggshell.r <- ps.r
taxo.eggshell <- taxo
taxo.eggshell.r <- taxo.r

rm(metadata, seqtab.nochim, taxa, S,
   sample.names, Srare, p1, tmp,
   comm, comm.r, mdata, mdata.r, ps, ps.r, taxo, taxo.r)

save.image("data/obs/1-phyloseq_objects_eggshell.RData")
