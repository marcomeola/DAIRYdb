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

# start doParallel
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
m <- matrix(rnorm(9), 3, 3)
foreach(i=1:nrow(m), .combine=rbind)

# Import database
setwd("~path/to/Database")
# import training sequences
dna <- readDNAStringSet("DAIRYdb_v1.1.2_10290_20180914_IDTAXA.fasta")

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

# train the classifier
trainingSet <- LearnTaxa(dna, taxonomy)
trainingSet

trainingSet2 <- LearnTaxa(dna,
          taxonomy,
          rank = NULL,
          K = floor(log(100*quantile(width(train), 0.99), 4)),
          minFraction = 0.01,
          maxFraction = 0.06,
          maxIterations = 10,
          multiplier = 100,
          maxChildren = 200,
          verbose = TRUE)

# view information about the classifier
plot(trainingSet)

# Analysis

# specify the path to the FASTA file (in quotes)
fas <- "/the/path/to/otus.fa"

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
              strand="both", # or "top" if same as trainingSet
              threshold=50, # 60 (very high) or 50 (high)
              processors=NULL) # use all available processors

# stop doParallel
stopCluster(cl)

# look at the results
print(ids)
plot(ids)

class(ids)
ids.m <- as.matrix(ids)

write.table(ids.m, file="taxonomy.IDTAXA")
#write.table(dna, file="dna")
#write.table(taxonomy, file="taxonomy")
#write.table(fas, file="fas")

# Saving workspace image
wd = "DAIRYdb_v1.1.2_20180914_IDTAXA"
filename = paste(wd,".RData",sep = "")
save.image(file = filename)

# Loading workspace image
#wd = "DAIRYdb_v1.1.2_20180914_IDTAXA" 
#filename = paste(wd,".RData",sep  =  "")
#load(filename)
