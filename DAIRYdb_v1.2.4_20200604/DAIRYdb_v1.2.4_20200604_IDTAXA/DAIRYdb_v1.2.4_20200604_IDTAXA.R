# Install DECIPHER
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

# load the DECIPHER library in R
library(DECIPHER)
library(seqinr)
library(XVector)
library(Biostrings)
library(doParallel)

# Import database
setwd("~/DAIRYdb_v1.2.4_20200604/DAIRYdb_v1.2.4_20200604_IDTAXA/")
# import training sequences
dna <- readDNAStringSet("DAIRYdb_v1.2.4_20200604_IDTAXA.fasta")

# parse the headers to obtain a taxonomy
s <- strsplit(names(dna), ";")
kingdom <- sapply(s, `[`, 1)
phylum <- sapply(s, `[`, 2)
class <- sapply(s, `[`, 3)
order <- sapply(s, `[`, 4)
family <- sapply(s, `[`, 5)
genus <- sapply(s, `[`, 6)
species <- sapply(s, `[`, 7)
taxonomy <- paste("Root", kingdom, phylum, class, order, family, genus, species, sep="; ")
head(taxonomy)

trainingSet <- LearnTaxa(dna,
          taxonomy,
          rank = NULL,
          K = floor(log(100*quantile(width(dna), 0.99), 4)),
          minFraction = 0.01,
          maxFraction = 0.06,
          maxIterations = 100000,
          multiplier = 100,
          maxChildren = 200,
          verbose = TRUE)

# view information about the classifier
plot(trainingSet)

trainingSet$problemGroups
trainingSet$problemSequences

# Analysis

# specify the path to the FASTA file (in quotes)
fas <- "FAI_all.fasta"
fas <- "priors.fa"

# load the sequences from the file
seqs <- readDNAStringSet(fas) # or readRNAStringSet

# remove any gaps (if needed)
seqs <- RemoveGaps(seqs)

# for help, see the IdTaxa help page (optional)
?IdTaxa

# load a training set object (trainingSet)
# see http://DECIPHER.codes/Downloads.html
load("<<REPLACE WITH PATH TO RData file>>")

# classify the sequences
ids <- IdTaxa(seqs,
              trainingSet,
              type = "extended",
              strand="both", # or "top" if same as trainingSet
              threshold=10, # 60 (very high) or 50 (high)
              bootstraps = 100,
              samples = L^0.47,
              minDescend = 0.98,
              processors=NULL, # use all available processors
              verbose = TRUE)

# look at the results
print(ids)
plot(ids)

class(ids)
ids.m <- as.matrix(ids)

getwd()
write.table(ids.m, file="taxonomy.IDTAXA")
#write.table(dna, file="dna")
#write.table(taxonomy, file="taxonomy")
#write.table(fas, file="fas")

# Saving workspace image
wd = "DAIRYdb_v1.2.4_20200604_IDTAXA"
filename = paste(wd,".RData",sep = "")
save(trainingSet, file = filename)

# Loading workspace image
#wd = "DAIRYdb_v1.1.2_20180914_IDTAXA" 
#filename = paste(wd,".RData",sep  =  "")
#load(filename)
