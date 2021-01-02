# for ICU ITS 
# dada2 for ITS sequencing 
# name : cuiqin 
# aim : dada2 protocol for ITS sequencing data and produce species table 

library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions # change the ref argument to get other versions
library(dada2)
packageVersion("dada2") #1.14
library(ShortRead)
packageVersion("ShortRead") #‘1.44.3’
library(Biostrings)
packageVersion("Biostrings") #‘2.54.0’
library(magrittr)
library(tidyverse)

load("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/dada.RData")

load("/media/lu/Lucy//documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/dada.RData")

# set WD
#setwd("/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS")
setwd("/home/lu/Downloads/ICU_data/result/")


# # data directory 
path <- "/home/lu/Downloads/ICU_data/ITS/"
#path <- "/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/data"
list.files(path) 

# GETTING READY -----------------------------------------------------------

# read names of fastq file names 
# Forward and reverse fastq filenames have format: SAMPLENAME_R1.raw.fastq.gz and SAMPLENAME_R2.raw.fastq.gz
fnFs <- sort(list.files(path, pattern="_1.rev.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.rev.fastq.gz", full.names = TRUE))

# IDENTIFY PRIMERS --------------------------------------------------------

FWD <- "GCATCGATGAAGAACGCAGC"  ## CHANGE ME to your forward primer sequence
REV <- "TCCTCCGCTTATTGATATGC"  ## CHANGE ME...

# fixing direction 
# If, for example, the REV primer is matching the Reverse reads in its RevComp
# orientation, then replace REV with its reverse-complement orientation (REV <- REV.orient[["RevComp"]]) before proceeding
# REV <- "GCATCGATGAAGAACGCAGCA"

# verification of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# pre filter the sequences just to remoce those with Ns
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# count the number of times the primers appear 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


# REMOVE PRIMERS ----------------------------------------------------------

#cutadapt <- '/Users/zilu/miniconda3/envs/cutadaptenv/bin/cutadapt' # CHANGE ME to the cutadapt path on your machine
cutadapt <- '/home/lu/conda/bin/cutadapt'
system2(cutadapt, args = "--version") # Run shell commands from R 1.18

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads) 
# I have an error message in later 

R1.flags <- paste("-g", FWD, "-a", REV.RC, "-m 1") 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "-m 1") 


# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:
# otuput should have all 0 values
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1.rev.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.rev.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)



# INSPECT READ QUALITY PROFILES -------------------------------------------

# forward reads
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutFs[1])

# reverse reads
plotQualityProfile(cutRs[1:2])
plotQualityProfile(cutRs[1])

# FILTER AND TRIM ---------------------------------------------------------

# Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE
head(out)

# LEARN THE ERROR RATES ---------------------------------------------------

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

# Dereplicate identical reads ---------------------------------------------
# Sample inference 
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)


# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# SAMPLE INFERENCE --------------------------------------------------------
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)


# MERGE PEIRED READS ------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# CONSTRUCT SEQUENCE TABLE ------------------------------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab) #74 9690


# REMOVE QUIMERAS ---------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE) #Identified 4684 bimeras out of 9690 input sequences.
#Inspect distribution of sequence lengths:
table(nchar(getSequences(seqtab.nochim)))

# seqtab.nochim is the final table, but colnames are sequences, we do the taxonomy step to have the taxonomy that correspond to the column seq

# TRACK READS THROUGH THE PIPELINE ----------------------------------------
# We now inspect the the number of reads that made it through each step in the pipeline to verify everything worked as expected.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

# ASSIGN TAXONOMY ---------------------------------------------------------
unite.ref <- "/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/UNITE_retrained_28062017/sh_general_release_dynamic_s_28.06.2017.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
# in taxa You have the seq and the taxonomy that correspond
head(taxa) 

save.image("/home/lu/Downloads/ICU_data/ITS/dada.RData")
load("/home/lu/Downloads/ICU_data/ITS/dada.RData")
# Inspecting the taxonomic assignments:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


# manipulate with plantae kingdom(remove them) 

OTU <- t(seqtab.nochim)
OTU_fungi <- cbind(OTU,taxa) %>% 
  tidyr::as_tibble() %>% 
  add_column(., id = rownames(OTU)) %>% 
  filter(., Kingdom == "k__Fungi") %>% 
  column_to_rownames(., var = "id") %>% 
  .[,c(1:74)] %>%
  type_convert(.)
   
#taxa_fungi <- taxa %>% 
#  as_tibble() %>% 
#  add_column(., id = rownames(OTU)) %>%
#  filter(., Kingdom == "k__Fungi") %>%
#  column_to_rownames(., var = "id") 



##### SPECIES #######
OTU_colnames_seq <- rownames(OTU_fungi)
species <- data.frame( taxa[OTU_colnames_seq,'Species'] )
genus <- data.frame( taxa[OTU_colnames_seq,'Genus'] )

species_full_name <- data.frame(paste(gsub('g__', '', genus$taxa.OTU_colnames_seq...Genus..) ,gsub('s_', '', species$taxa.OTU_colnames_seq...Species..), sep = '' ))

colnames(species_full_name) <- 'species'
OTU_with_species <- data.frame( species_full_name, OTU_fungi) 
OTU_with_species <- OTU_with_species %>% na.omit()
test1<- OTU_with_species


OTU <- t(seqtab.nochim)
OTU_colnames_seq <- rownames(OTU)
species <- data.frame( taxa[OTU_colnames_seq,'Species'] )
genus <- data.frame( taxa[OTU_colnames_seq,'Genus'] )
species_full_name <- data.frame(paste(gsub('g__', '', genus$taxa.OTU_colnames_seq...Genus..) ,gsub('s_', '', species$taxa.OTU_colnames_seq...Species..), sep = '' ))
colnames(species_full_name) <- 'species'
OTU_with_species <- data.frame( species_full_name, OTU) 

OTU_with_species <- OTU_with_species %>% 
  na.omit() %>%
  rownames_to_column(., var = "id") %>% 
  filter(., species != "NANA") %>%
  column_to_rownames(., var = "id")  

test2<-OTU_with_species

# these are the species not in fungi kindom manually remove them 
non_fungi <- setdiff(test2$species,test1$species)
#[1] "Trebouxia_decolorans"           "TrebouxiaNA"                    "QuercusNA"                      "Paramicrosporidium_saccamoebae"
#[5] "Trebouxia_arboricola"

setdiff(test1$species,test2$species)


# I want to sum the species repeated
group_OTU_table_tax <- group_by(OTU_with_species, species) 
OTU_subset_species <- summarise_all(group_OTU_table_tax, sum)
table_species <- data.frame(OTU_subset_species)
rownames(table_species) <- table_species$species
table_species_DADA2 <- table_species[,-1]
table_species_DADA2 %>% dim()
colnames(table_species_DADA2) <- gsub('_1.rev.fastq.gz', '', colnames(table_species_DADA2))

# manually remove those are not in fungi kingdom 

test3 <- table_species_DADA2 %>% rownames_to_column(.,var = "species") %>%
  filter(.,!species %in% non_fungi) 
  

table_species_DADA2 <- test3 
# save file
write.csv2(table_species_DADA2, file = '/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/table_species_DADA2_human.csv')

phylum <- data.frame( taxa[OTU_colnames_seq,'Phylum'] )
species_with_matched_phylum <- data.frame(species = species_full_name,phylum = phylum$taxa.OTU_colnames_seq...Phylum..)
write.csv2(species_with_matched_phylum, file = '/Volumes/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/species_matched_phylum.csv')


save.image("dada_output.RData")
load("/media/lu/Lucy/documents/201910_ICU_ITS/74samples_ITS/dada/dada_ITS/dada_output.RData")
