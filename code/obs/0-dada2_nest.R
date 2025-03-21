# Description -------------------------------------------------------------
# Observational approach
# Effect of aromatic plants on nest bacterial microbiota
# DADA2 analysis
# Sequencing: Illumina mi-seq June 2024
# Authors: Hélène Dion-Phénix & Gabrielle Gingras
# Last edition: 2025-02-28

# Libraries --------------------------------------------------------------

library(dada2); packageVersion("dada2")

# DADA2 -------------------------------------------------------------------

#Set file paths which contain unzip fastq files and extract sample names

#Set file paths
path <- "data/obs/bioinfo/fastq/nest"
# Forward and reverse fastq filenames have format: SAMPLENAME_L001_R1_001.fastq and SAMPLENAME_L001_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R2_", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R1_", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names

# Plot quality profile of forward reads
plotQualityProfile(fnFs[1:9])
# Plot quality profile of reverse reads
plotQualityProfile(fnRs[1:9])

#Trimming
#Primer: 
#AACMGGATTAGATACCCKG  ->    799F
#AGGGTTGCGCTCGTTG   ->   1115R
#need to trim 19 from start of foward reads and 16 from start of reverse reads (gets rid of primer)
#Quality trimming:
#Forward reads quality crashes around 160
#Reverse reads quality crashes around 200


fowardCutOff <- 160
reverseCutOff <- 200
# set filtered file folder path
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(19,16), 
                     truncLen = c(fowardCutOff, reverseCutOff),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE)

print(out[1:10])

# Learn error rates
set.seed(100)
errF <- learnErrors(filtFs)
errR <- learnErrors(filtRs)

#sanity check - visualize error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Dereplicate the filtered fastq files
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF)
dadaRs <- dada(derepRs, err=errR)

# e.g. inspect results
dadaFs[[1]]
dadaRs[[1]]

#Merge the denoised forward and reverse reads:
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                      returnRejects = TRUE, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[10]])

#construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab) 
table(nchar(getSequences(seqtab)))

#save sequence table
saveRDS(seqtab,
        "data/obs/bioinfo/seqtab_nest.rds")

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Remove Name. and .fastq.gz from community names so they match
rownames(seqtab.nochim) <- gsub(".fastq", "", rownames(seqtab.nochim))
rownames(seqtab.nochim)

#Summary - track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab), 
               rowSums(seqtab.nochim))

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#To get the number of sequences 
track <- as.data.frame(track)
sum(track$input) 
sum(track$filtered) 
sum(track$nonchim) 

# identify taxonomy
taxa <- assignTaxonomy(seqtab.nochim, 
                       "data/SILVA138.1/silva_nr99_v138.1_train_set.fa.gz",
                       multithread=TRUE, tryRC=TRUE)

# exact species matching (won't work if concatenating sequences)
taxa <- addSpecies(taxa,  
                   "data/SILVA138.1/silva_species_assignment_v138.1.fa.gz",
                   allowMultiple = TRUE, tryRC = TRUE)

#Inspect the taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
capture.output(taxa.print, file ="data/obs/bioinfo/taxa_nest_dada2.txt")

# Clean and save environment ----------------------------------------------

#Save environment
save.image("data/obs/bioinfo/0-dada2_nest_complete.RData")

#Save output - only seqtab.nochim, taxa, track and sample.names
rm(list=c('dadaFs',"dadaRs","derepFs","derepRs","errF","errR","mergers",
          "out","seqtab","taxa.print","filt_path","filtFs","filtRs",
          "fnFs","fnRs","fowardCutOff","path","reverseCutOff","getN",
          "track"))

save.image("data/obs/bioinfo/0-dada2_nest.RData")
