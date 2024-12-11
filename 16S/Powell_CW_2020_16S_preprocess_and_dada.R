# 16S Analysis of Curacao 2020 and 2021 data
#Maya Powell
#March 2023 and ongoing
###Based on DADA2 Pipeline 1.16, Ana Dulskiy's oculina 16S analysis, and Nicola Kriefall's Moorea holobiont analysis

#~########################~#
##### PRE-PROCESSING #######
#~########################~#
#all pre-processing done on longleaf ondemand which is the UNC computing cluster

#make fasta files
#my adapters for 16S, which I saved as "adapters.fasta"
>forward
AATGATACGGCGACCAC
>forwardrc
GTGGTCGCCGTATCATT
>reverse
CAAGCAGAAGACGGCATAC
>reverserc
GTATGCCGTCTTCTGCTTG

#primers for 16S, which I saved as "primers.fasta": 
>forward
GTGYCAGCMGCCGCGGTA
>reverse
GGACTACHVGGGTWTCTAAT

#Open terminal & sign in to longleaf (either online on ondemand or log in to longleaf from home terminal)
#Navigate to your folder where you want to do your preprocessing
#put all files from sequencing facility into this folder
#these should be in the fastq.gz format with R1 = forward reads and R2 = reverse reads

#load bbmap program
load bbmap

#In terminal - un gzip your files
gzip -d *fastq.gz

chmod g+w -R /proj/kdcastil/users/mayapow/CW_RTE_22/ <- gives someone access to the folder

## Still in terminal - making a sample list based on the first phrase before the underscore in the .fastq name
ls *R1*.fastq | cut -d '_' -f 1 > samples.list

#View sample files using any of these commands:
# cat samples.list
# nl
# less (press q to exit) (this one is my fave personally)
# head = first 10 lines, tail = last 10 lines

#rename files
#this didn't work for any of the 1s (A1, B1) so had to rename those manually
#and I'm sorry I didn't figure out why xoxo hope your samples work better
for file in $(cat samples.list); do mv ${file}*R1*.fastq ${file}_R1.fastq; mv ${file}*R2*.fastq ${file}_R2.fastq; done

## Get rid of reads that still have the adapter sequence, shouldn't be there, I didn't have any
for file in $(cat samples.list); do bbduk.sh in1=${file}_R1.fastq in2=${file}_R2.fastq ref=adapters.fasta out1=${file}_R1_NoIll.fastq out2=${file}_R2_NoIll.fastq; done &>bbduk_NoIll.log
#yes removes 0 of anything nice

## Get rid of first 4 bases (degenerate primers created them)
for file in $(cat samples.list); do bbduk.sh in1=${file}_R1_NoIll.fastq in2=${file}_R2_NoIll.fastq ftl=4 out1=${file}_R1_NoIll_No4N.fastq out2=${file}_R2_NoIll_No4N.fastq; done &>bbduk_No4N.log
#should remove ~1.6 % of bases and 0% of reads
#yes did this exactly as we'd expect

## Only keep reads that start with the 16S primer
## higher k = more reads removed, but can't surpass primer length (here primers are 19 and 21)
#used k = 10 here with restrictleft=21, primers are 19 and 21
for file in $(cat samples.list); do bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq restrictleft=21 k=10 literal=GTGYCAGCMGCCGCGGTAA,GGACTACNVGGGTWTCTAAT copyundefined=t outm1=${file}_R1_NoIll_No4N_16S.fastq outu1=${file}_R1_check.fastq outm2=${file}_R2_NoIll_No4N_16S.fastq outu2=${file}_R2_check.fastq; done &>bbduk_16S.log

#ITS2 processing
#This takes all files that match to SYM_VAR_5.8S2 (FWD) and SYM_VAR_REV (REV) primers
#in your log file - all the "reads removed" are the ones saved in your outm files (ITS2_final.fastq)
#used k = 20 here, with restrictleft=25, primers are 21bp and 25bp
for file in $(cat samples.list); do bbduk.sh in1=${file}_R1_NoIll_No4N.fastq in2=${file}_R2_NoIll_No4N.fastq k=20 restrictleft=25 literal=GAATTGCAGAACTCCGTGAACC,CGGGTTCWCTTGTYTGACTTCATGC outm1=${file}_R1_ITS2_final.fastq outu1=${file}_R1_ITS2_check.fastq outm2=${file}_R2_ITS2_final.fastq outu2=${file}_R2_ITS2_check.fastq; done &>bbduk_ITS2.log

#Copy the ITS2 final files into a new folder
scp -r /proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/CW_2020_lane2/*ITS2_final.fastq /proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/CW_2020_ITS2
scp -r /proj/kdcastil/users/mayapow/CW_RTE_22/CW_RTE_June_23/*ITS2_final.fastq /proj/kdcastil/users/ander/CW_RTE_ALL/CW_RTE_June_23
#navigate into your ITS2 folder
#re-gzip your files (in your ITS2 folder)
#this compresses the files which makes it easier to download and makes it into the format that SymPortal wants
gzip *.fastq
#Download the ITS2 folder to your local computer from ondemand, or you can scp if using your computers terminal
#Submit to SymPortal!! Yay!!

#now back to 16S
#add cutadapt
module load cutadapt
cutadapt --version
#version 4.2 mp 3/18/22

## Use cutadapt to remove 16S primer
for file in $(cat samples.list); do cutadapt -g GTGYCAGCMGCCGCGGTAA -a ATTAGAWACCCBNGTAGTCC -G GGACTACHVGGGTWTCTAAT -A TTACCGCGGCKGCTGRCAC -n 2 --discard-untrimmed -o ${file}_R1_16S_final.fastq -p ${file}_R2_16S_final.fastq ${file}_R1_NoIll_No4N_16S.fastq ${file}_R2_NoIll_No4N_16S.fastq; done &>clip.log

#info about cutadapt presets:
##-g regular 5' forward primer 
##-G regular 5' reverse primer
##-o forward out
##-p reverse out
##-max-n 0 means 0 Ns allowed

#move 16S_final.fastq files into a new 16S folder for analysis
scp -r /proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/CW_2020_lane2/*16S_final.fastq /proj/kdcastil/users/mayapow/CW_2020_TUCF_ALL/CW_2020_16S
#Download 16S folder to your home computer, and put it in a folder where you will do your DADA2 analysis (along with this R script)
#you can gzip them if its taking 5years long which is what I did, and then just unzip them in terminal on your local computer

scp -r /proj/kdcastil/users/ander/ITS2_2024/*16S_final.fastq /proj/kdcastil/users/mayapow/CW_June_2022_16S

scp -r /proj/kdcastil/users/ander/ITS2_2024/*16S_final.fastq /proj/kdcastil/users/mayapow/CW_June_2022_16S

#~########################~#
##### DADA2 BEGINS #########
#~########################~#

#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.16")
library(dada2); packageVersion("dada2")
#Version 1.26.0

#install.packages("ShortRead")
library(ShortRead); packageVersion("ShortRead")
#1.56.1

library(Biostrings); packageVersion("Biostrings")
#2.66.0

path <- "~/Documents/Castillo Lab/CW_2020/CW_2020_16S" 
# CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = "_R1_16S_final.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern = "_R2_16S_final.fastq", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))

#### check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME to your reverse primer sequence

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
FWD.orients
REV.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress=TRUE)
#this step above took me a LONG TIME like 15+mins for 158 samples
#should say "creating output directory" - if it doesn't you probably already have the directory, and need to delete the filtN directory & redo
#Some input samples had no reads pass the filter.

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[4]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[4]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[4]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[4]]))
#no hits for primers yay!! if there are hits - do cutadapt part 2

#### Visualizing raw data ####

#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs.filtN[c(1:12)]) #225, all high read counts
plotQualityProfile(fnFs.filtN[c(13:24)]) #B8 175 with 5,000 reads, all else 225 high reads
plotQualityProfile(fnFs.filtN[c(25:36)]) #C2-C12 all low 2-9,000 reads, all good to 225
plotQualityProfile(fnFs.filtN[c(37:48)]) #225, D2, D5, D9 good, all else low
plotQualityProfile(fnFs.filtN[c(49,51,52,53,54,56,58)]) #none made it lol, all bad
plotQualityProfile(fnFs.filtN[c(61:72)]) #225, all high read counts
plotQualityProfile(fnFs.filtN[c(73:84)]) #225, all high read counts, G3-G5 drop a bit at 175 but not much
plotQualityProfile(fnFs.filtN[c(85:96)]) #225, all high read counts, H12 = H20 =2,000 reads
plotQualityProfile(fnFs.filtN[c(97:108)]) #225, super super high read counts, all 90,000 + except J1
plotQualityProfile(fnFs.filtN[c(109:120)]) #K2 & K10 drop at 175, K4 3,900 reads
plotQualityProfile(fnFs.filtN[c(121:132)]) #L1 703 reads, L10 7,000 reads, L3 5,900 reads, L6 drops at 175
plotQualityProfile(fnFs.filtN[c(133:144)]) #L7 9,900, L8 5,300, M1 9,700, M2 3,500 M5 1,900, M6 684
plotQualityProfile(fnFs.filtN[c(145:158)]) #lots drop at 175, M7 889 reads, N1 954 reads, N12 = H20 = 369 reads

#Then look at quality profile of R2 reads
plotQualityProfile(fnRs.filtN[c(1:12)]) #225 all good, some drop a bit at 175 but not bad
plotQualityProfile(fnRs.filtN[c(13:24)]) #B1 ONLY ONE READ WHAT THE FUCK
plotQualityProfile(fnRs.filtN[c(25:36)]) #Cs kinda low 3-9,000
plotQualityProfile(fnRs.filtN[c(37:48)]) #D's same as Cs
plotQualityProfile(fnFs.filtN[c(49,51:54,56,58)]) #none made it lol, all bad
plotQualityProfile(fnRs.filtN[c(61:72)]) #225, all great
plotQualityProfile(fnRs.filtN[c(73:84)]) #225, all great
plotQualityProfile(fnRs.filtN[c(85:96)]) #225 all great, H12=H20=2,000
plotQualityProfile(fnRs.filtN[c(97:108)]) #all great, 225 some drop at 175 a little but ok
plotQualityProfile(fnRs.filtN[c(109:120)]) #some drop at 175, K5 3,900
plotQualityProfile(fnRs.filtN[c(121:132)]) #225, L1 703 reads, L11, L4 5-7,000 reads
plotQualityProfile(fnRs.filtN[c(133:144)]) #175 for some, must cutoff here, M4 684, M6 1953, L8 5322, M10 9773, M2 3500
plotQualityProfile(fnRs.filtN[c(145:158)]) #M7 889, N1 954, N2 369, N12 112923?!?!?!?! wtf not possible based on clip log, must be some weird naming thing here

# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#maxEE = 2,2 - allows 2 expected errors
#Trunclen = 175
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(175,175), #leaves enough overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(2,2), #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)
write.csv(out,file="out.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Learn Error Rates ########
#~############################~#

setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE) #start 5:34
#107390150 total bases in 613658 reads from 11 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#107390150 total bases in 613,658 reads from 11 samples will be used for learning the error rates.
#has to be incorrect - they shouldn't be the same here...something is weird! this happened with my initial reads too

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score

plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 
#these look different though so maybe it is real

#~############################~#
##### Dereplicate reads ########
#~############################~#

#Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence. 
#Dereplication substantially reduces computation time by eliminating redundant comparisons.
#DADA2 retains a summary of the quality information associated with each unique sequence. The consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads. These quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2’s accuracy.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#~###############################~#
##### Infer Sequence Variants #####
#~###############################~#

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 

dadaFs[[1]]
dadaRs[[1]]

#~############################~#
##### Merge paired reads #######
#~############################~#

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

#~##################################~#
##### Construct sequence table #######
#~##################################~#
#a higher-resolution version of the “OTU table” produced by classical methods

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#146 samples 39280 sequences
saveRDS(seqtab, file="seqtab_CW_16S.rds")
write.csv(seqtab, file="seqtab_CW_16S.csv")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

plot(table(nchar(getSequences(seqtab)))) #real variants appear to be right in that 240-270 window

#The sequence table is a matrix with rows corresponding to (and named by) the samples, and 
#columns corresponding to (and named by) the sequence variants. 
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, and may be worth removing

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(240,270)] #again, being fairly conservative with length
saveRDS(seqtab2, file="seqtab2_CW_16S.rds")
write.csv(seqtab2, file="seqtab2_CW_16S.csv")

#~############################~#
##### Remove chimeras ##########
#~############################~#
#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) #146 30605
#Identified 8474 bimeras out of 39079 input sequences.
saveRDS(seqtab.nochim, file="seqtab.nochim_CW_16S.rds")
write.csv(seqtab.nochim, file="seqtab.nochim_CW_16S.csv")

sum(seqtab.nochim)/sum(seqtab2)
#0.9297248
#The fraction of chimeras varies based on factors including experimental procedures and sample complexity, 
#but can be substantial. 
#looks good! Keeping >92%

#~############################~#
##### Track Read Stats #########
#~############################~#

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

setwd("~/Documents/Castillo Lab/CW_2020/CW_2020_16S")
write.csv(track,file="CW_2020_16S_readstats.csv",row.names=TRUE,quote=FALSE)

#~############################~#
##### Assign Taxonomy ##########
#~############################~#

#Using package DECIPHER as an alternative to 'assignTaxonomy'
#BiocManager::install("DECIPHER")
library(DECIPHER); packageVersion("DECIPHER")
# version 2.26.0
#citation("DECIPHER")
#http://DECIPHER.codes/Downloads.html. Download the SILVA most updated files to follow along.
#Silva info and download: https://www.arb-silva.de/

#Assign Taxonomy
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/Documents/curacao_2020/bin/16S_fully_processed/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/curacao_2020/bin/16S_fully_processed/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
#this will take a long long time!!
saveRDS(taxa, file="CW_2020_16S_taxa.rds") 

unname(head(taxa))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa <- readRDS("CW_2020_16S_taxa.rds")

taxa.plus <- addSpecies(taxa, "~/Documents/curacao_2020/bin/16S_fully_processed/silva_species_assignment_v138.1.fa.gz",tryRC=TRUE,verbose=TRUE)

saveRDS(taxa.plus, file="CW_2020_16S_taxa.plus.rds") 
#365 out of 30605 were assigned to the species level.
#Of which 313 had genera consistent with the input table.

write.csv(taxa.plus, file="CW_2020_16S_taxa.plus.csv")
write.csv(taxa, file="CW_2020_16S_taxa.csv")

#### Read in previously saved datafiles ####
seqtab.nochim <- readRDS("seqtab.nochim_CW_16S.rds")
taxa <- readRDS("CW_2020_16S_taxa.rds")
taxa.plus <- readRDS("CW_2020_16S_taxa.plus.rds")

#~############################~#
##### handoff 2 phyloseq #######
#~############################~#

BiocManager::install("phyloseq")
library('phyloseq')
library('ggplot2')
library('Rmisc')
library(cowplot)
library(ShortRead)

#import dataframe holding sample information
samples.out <- rownames(seqtab.nochim)
samdf<-read.csv("sample_info_curacao_2020.csv")
head(samdf)
samdf <- samdf %>% select(-X)
rownames(samdf) <- samdf$id

#Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa.plus))

ps 
#otu_table()   OTU Table:         [ 30605 taxa and 146 samples ]
#sample_data() Sample Data:       [ 146 samples by 24 sample variables ]
#tax_table()   Taxonomy Table:    [ 30605 taxa by 7 taxonomic ranks ]

#### first look at data ####
ps_glom <- tax_glom(ps, "Family")
plot_bar(ps_glom, x="site_zone", fill="Family")+
  theme(legend.position="none")

#phyloseq object with shorter names - doing this one instead of one above
ids <- paste0("sq", seq(1, length(colnames(seqtab.nochim))))
#making output fasta file - this will be blasted for taxa in a bit!
path='~/Documents/Castillo Lab/CW_2020/CW_2020_16S/CW_2020_16S_taxa.fasta'
uniquesToFasta(seqtab.nochim, path, ids = ids, mode = "w", width = 20000)

#### blast asvs to NCBI to see if any eukaryotes got through ####
##Running blast on longleaf to make organism match files for my 16s data
##using the fasta file made above -> CW_2020_16S_taxa.fasta

## On longleaf or your university computing cluster:

#create script
nano taxa_blast.sh

#write your bash script:

#!/bin/bash

#SBATCH -p general #this will run it on the general node
#SBATCH -t 10-00:00:00 #time limit is 10 days (make sure to include)
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem=30gb
#SBATCH --output <CW_2020_16S_taxids.out>

module load blast
blastn -query CW_2020_16S_taxa.fasta -db nt -outfmt "6 std staxids sskingdoms" -evalue 1e-5 -max_target_seqs 5 -out CW_2020_16S_taxids.out

#exit out and save
#run this command below in terminal
#your formatting for this script may look quite different based on what your university or computing cluster does
sbatch blast_taxid.sh

#do this to check if your code is running
#only relevant for UNC longleaf computing cluster
squeue -u mayapow 
#-u = username

#what the prompts in this script mean for reference:
#SBATCH -p general #this will run it on a general node on the server
#SBATCH -t 10-00:00:00 #time limit is 10 days (make sure to include seconds!)
#SBATCH -N 1 (this will run it on 1 single node)
#SBATCH -n 40 (this will use 40 "CPUs" to run it)
#SBATCH --mem=30gb (specify how much memory)
#SBATCH --output <CW_2020_16S_taxids.out> (your output log file)
#other versions of this script have used the "-remote" flag to access the actual NCBI remote database
#you can add -v or -vv if you want it to be verbose and give u more information in the log file
#however, longleaf has a local version of the nt blast database (through NCBI) that stays updated, so we can use this
#the remote version will take on the order of days to run, while the local version takes ~5-10hrs

###Remove contaminants
#Now continue and pull out euks yourself here and update phyloseq objects
#also remove negative control contaminants
colnames(seqtab.nochim)<-ids
taxa2 <- cbind(taxa.plus, rownames(taxa.plus)) #retaining raw sequence info before renaming
rownames(taxa2)<-ids

#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2))
ps
#otu_table()   OTU Table:         [ 30605 taxa and 146 samples ]
#sample_data() Sample Data:       [ 146 samples by 24 sample variables ]
#tax_table()   Taxonomy Table:    [ 30605 taxa by 8 taxonomic ranks ]

#### remove mitochondria, chloroplasts, non-bacteria #### 
ps.mito <- subset_taxa(ps, (Family=="Mitochondria"))
ps.mito #1015 taxa to remove
ps.chlor <- subset_taxa(ps, (Order=="Chloroplast"))
ps.chlor #734 taxa to remove
ps.notbact <- subset_taxa(ps, (Kingdom!="Bacteria") | is.na(Kingdom))
ps.notbact #1101 taxa to remove

ps.nomito <- subset_taxa(ps, (Family!="Mitochondria") | is.na(Family))
ps.nomito #29590 taxa
ps.nochlor <- subset_taxa(ps.nomito, (Order!="Chloroplast") | is.na(Order))
ps.nochlor #28856 taxa
ps.clean <- subset_taxa(ps.nochlor, (Kingdom=="Bacteria"))
ps.clean #27755 taxa

#just archaea
ps.arch <- subset_taxa(ps.nomito, (Kingdom=="Archaea"))
ps.arch #1084 taxa

#### identifying contamination ####
BiocManager::install("decontam")
library(decontam)

df <- as.data.frame(sample_data(ps.clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=site)) + 
  geom_point()

ggplot(data=df, aes(x=Index, y=LibrarySize, color=collection_month, group = site)) + 
  geom_point(aes(shape = site))

sample_data(ps.clean)$is.neg <- sample_data(ps.clean)$full_sample_id == "blank"
contamdf.prev <- isContaminant(ps.clean, neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)

#FALSE  TRUE 
#27751     4 
#4 contaminants - remove

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps.clean, function(abund) 1*(abund>0))
#ps.pa.neg <- prune_samples(sample_data(ps.pa)$site == "neg", ps.pa) #none
ps.pa.pos <- prune_samples(sample_data(ps.pa)$full_sample_id != "blank", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove from ps.clean:
ps.clean1 <- prune_taxa(!contamdf.prev$contaminant,ps.clean)
#also remove negative controls, don't need them anymore
ps.cleaner <- subset_samples(ps.clean1,(full_sample_id!="blank"))
ps.cleaner
#27751 taxa and 144 samples (started with 146 - 2 for negative controls)
saveRDS(ps.cleaner, file ="CW_2020_16S_ps.cleaner.RDS")
saveRDS(taxa, file ="CW_2020_16S_taxa.Rdata")
saveRDS(taxa.plus, file ="CW_2020_16S_taxa.plus.Rdata")
saveRDS(taxa2, file ="CW_2020_16S_taxa2.Rdata")
ps.cleaner <- readRDS("CW_2020_16S_ps.cleaner.RDS")

###STOP HERE####
####YOU MUST WAIT FOR YOUR BLAST TO FINISH ON YOUR SERVER###
###ONLY CONTINUE ONCE YOU HAVE YOUR CW_2020_16S_taxids.out FILE####
#this took 4.5 days to run for me phew! with 30,000 ASVs

##now getting taxonomy info:

#in longleaf
#navigate to your folder where you ran the blast and your output file CW_2020_16S_taxids.out ended up
module load taxonkit

##extract taxa ids from blast output for taxonkit (in terminal):
awk -F " " '{print $13}' CW_2020_16S_taxids_copy.out > ids
taxonkit lineage ids > ids.tax
cut -f1 CW_2020_16S_taxids_copy.out > ids.seq; paste ids.seq ids.tax > ids.seq.tax
grep "Eukaryota" ids.seq.tax | cut -f1 | sort | uniq > euk.contam.asvs

#download euk.contam.asvs to your local computer, and put it in your working directory
#Moving back to R
setwd("~/Documents/Castillo Lab/CW_2020/CW_2020_16S")

#load in phyloseq objects from above
ps.cleaner <- readRDS("CW_2020_16S_ps.cleaner.RDS")
taxa2 <- readRDS("CW_2020_16S_taxa2.Rdata")

# convert euk.contam.asvs to .csv file
eukasvs <- read_table("euk.contam.asvs", col_names = FALSE)
write_csv(eukasvs, file = "euk.contam.asvs.csv", col_names = FALSE)

#remove these euk contaminants from ps.cleaner
# write # to remove here
euks <- read.csv("euk.contam.asvs.csv",header=FALSE)
euks_names <- euks$V1
alltaxa <- taxa_names(ps.cleaner) #27751
`%!in%` <- Negate("%in%")
keepers <- alltaxa[(alltaxa %!in% euks_names)] #keepers = 27653, so that means only 98 euks got through above
ps.cleanest <- prune_taxa(keepers, ps.cleaner) 
ps.cleanest
#27653 taxa and 144 samples in ps.cleanest

#save and re-read in cleanest phyloseq objects + seqtabs
seqtab.cleanest <- data.frame(otu_table(ps.cleanest))
write.csv(seqtab.cleanest,file="CW_2020_16S_seqtab.cleanest.csv")
saveRDS(ps.cleanest,file="CW_2020_16S_ps.cleanest.RDS")

seqtab.cleanest <- read.csv("CW_2020_16S_seqtab.cleanest.csv",row.names=1)
ps.cleanest <- readRDS("CW_2020_16S_ps.cleanest.RDS")
#quick reminder about what different files are:
#ps.clean = ps that has had all the SILVA database euks and other non-bacteria removed
#ps.cleaner = ps.clean that has had negative control contaminants removed
#ps.cleanest = ps.cleaner that has had all the euks from NCBI blast removed
#probably PS cleanest is what I will want to do my analyses with if I choose not to rarefy, but more on that later

#### rarefy #####
library(vegan); packageVersion("vegan")
#version 2.6.4

seqtab.cleanest <- data.frame(ps.cleanest@otu_table) #27653
samdf.cleanest <- data.frame(ps.cleanest@sam_data)
seqtab.cleaner <- data.frame(ps.cleaner@otu_table) #27751
samdf.cleaner <- data.frame(ps.cleaner@sam_data)

rarecurve(seqtab.cleanest,step=1000,label=FALSE) 
#seqtab cleaner after removing contaminants + euk blast NCBI
rarecurve(seqtab.cleaner,step=1000,label=FALSE)
#after removing contaminants based on SILVA only & negative controls

total <- rowSums(seqtab.cleanest)
#total <- rowSums(seqtab.cleaner)
#make some plots to figure out where to rarefy
total.df <- data.frame(total)
ggplot(total.df, aes(x = total)) +
  geom_histogram()

ggplot(total.df, aes(x=1, y=total))+
  geom_boxplot()

ggplot(total.df, aes(x=1, y=total))+
  geom_jitter()+
  scale_y_log10()

library(dplyr)
total <- rowSums(seqtab.cleanest)
total.df %>%
  arrange(total) %>%
  ggplot(aes(x=1:nrow(.),y=total))+
  geom_line()  #shows # of samples on x and # of seqs in each sample on y

total.df %>%
  arrange(total) %>%
  print(20)
subset(total, total <1000)
#samples range from 89 to 135,232
#YEESH
#not sure what to do about this - how can I even compare samples between?

subset(total, total <1000) #7 samples need to be removed (H11,H12,H5,H9,H10,H8,H2)

#get rid of all samples with less than 1000 seqs
#6 samples
# identified by MCMC.OTU below as being too low 

row.names.remove <- c("K4","K9","M3","M5","M6","M9","N10","N4","N8")
seqtab.less <- seqtab.cleanest[!(row.names(seqtab.cleanest) %in% row.names.remove),]
#seqtab.less <- seqtab.cleaner[!(row.names(seqtab.cleaner) %in% row.names.remove),]
samdf.less <- samdf.cleanest[!(row.names(samdf.cleanest) %in% row.names.remove), ]
#samdf.less <- samdf.cleaner[!(row.names(samdf.cleaner) %in% row.names.remove), ]
dim(samdf.less) #135 samples left
ps.less <- phyloseq(otu_table(seqtab.less, taxa_are_rows=FALSE), 
                    sample_data(samdf.less), 
                    tax_table(taxa2))
###saving
#setwd("~/oculina/data")
write.csv(seqtab.less, file="CW_2020_16S_seqtab.less.csv")
write.csv(samdf.less, file="CW_2020_16S_samdf.less.csv")
saveRDS(ps.less,file="CW_2020_16S_ps.less.RDS")

#for rarefied data, removing values up to 10,000
subset(total, total <10000)
row.names.remove <- c("B10","B8","B9","C10","C11","C12","C2","C3","C5","C6","C7","C8","C9","D1","D10","D11","D12","D2","D3","D4","D6","D7","D8","G3","H11","H9","J7","K4","K9","L10","L3","L6","L8","L9","M1","M2","M3","M4","M5","M6","M9","N10","N4","N7","N8","N9")
seqtab.r <- seqtab.cleanest[!(row.names(seqtab.cleanest) %in% row.names.remove),]
#seqtab.less <- seqtab.cleaner[!(row.names(seqtab.cleaner) %in% row.names.remove),]
samdf.r <- samdf.cleanest[!(row.names(samdf.cleanest) %in% row.names.remove), ]
#samdf.less <- samdf.cleaner[!(row.names(samdf.cleaner) %in% row.names.remove), ]
dim(samdf.r) #98 samples left
seqtab.rare <- rrarefy(seqtab.r,sample=10000) #all samples less than 40,000 not rarefied
rarecurve(seqtab.rare,step=1000,label=FALSE)
#phyloseq object but rarefied
ps.rare <- phyloseq(otu_table(seqtab.rare, taxa_are_rows=FALSE), 
                    sample_data(samdf.rare), 
                    tax_table(taxa2))
ps.rare #3433 taxa, 64 samples

#removing missing taxa - lost after rarefying
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare
#23759 taxa and 98 samples

seqtab.rare <- data.frame(otu_table(ps.rare))
samdf.rare <- data.frame(sample_data(ps.rare))

#saving
#### data files - rarefied, decontaminated ####
saveRDS(ps.rare, file="CW_2020_16S_ps.rare.RDS")
write.csv(seqtab.rare,file="CW_2020_16S_seqtab.rare.csv")
write.csv(samdf.rare,file="CW_2020_16S_samdf.rare.csv")

#### trim underrepresented otus ####
install.packages("MCMC.OTU")
library(MCMC.OTU)

#formatting the table for mcmc.otu - requires one first column that's 1 through whatever
#& has "X" as column name
nums <- 1:nrow(seqtab.less)  #originally used seqtab.cleanest but I think seqtab.less better
samples <- rownames(seqtab.less)

int <- cbind(sample = 0, seqtab.less) #27654
seq.formcmc <- cbind(X = 0, int) #27655

seq.formcmc$X <- nums
seq.formcmc$sample <- samples

#change second columns value to equal total variables
dim(seq.formcmc) #135 by 27655
seq.trim.allinfo <- purgeOutliers(seq.formcmc,count.columns=3:27655,sampleZcut=-2.5,otu.cut=0.0001,zero.cut=0.02)
#1 sample with counts below z-score
#1006 ASVs pass filters using seqtab.less

#[1] "samples with counts below z-score -2.5 :"
#[1] "B8"
#[1] "zscores:"
#B8 
#-2.65978 
#[1] "OTUs passing frequency cutoff  1e-04 : 1006"
#[1] "OTUs with counts in 0.02 of samples:"
#FALSE  TRUE 
#171   835 

#remove sample info
dim(seq.trim.allinfo) #134 by 837
seq.trim <- seq.trim.allinfo[,3:837] #columns of sequences in seq.trim.allinfo

write.csv(seq.trim,file="CW_2020_16S_seqtab.trim.csv")
seq.trim <- read.csv("CW_2020_16S_seqtab.trim.csv",row.names=1)

#remake phyloseq objects
ps.trim <- phyloseq(otu_table(seq.trim, taxa_are_rows=FALSE), 
                    sample_data(samdf.less), 
                    tax_table(taxa2))
ps.trim #835 taxa and 134 samples

saveRDS(ps.trim,file="CW_2020_16S_ps.trim.RDS")


#### rarefy - trimmed #####
library(vegan)

seqtab.trim <- data.frame(ps.trim@otu_table)
samdf.trim <- data.frame(ps.trim@sam_data)
write.csv(samdf.trim, file="CW_2020_16S_samdf.trim.csv")

rarecurve(seqtab.trim,step=100,label=FALSE) 

total <- rowSums(seqtab.trim)
subset(total, total <3498) #smallest is 1246 - try to rarefy to this and see what we get
#21 samples

row.names.remove <- c("C10",  "C11",   "C2",   "C3",   "C5",   "D4",   "D8",  "H11",  "L10",   "L3",   "L8",   "M2",   "M4",   "N9")
seqtab.t.r <- seqtab.trim[!(row.names(seqtab.trim) %in% row.names.remove),]

seqtab.trim.rare <- rrarefy(seqtab.t.r,sample=3498)
rarecurve(seqtab.trim.rare,step=100,label=FALSE)

#look at rarefaction curves for all data
rarecurve(seqtab.less,step=100,label=FALSE)
rarecurve(seqtab.rare,step=100,label=FALSE)
rarecurve(seqtab.trim,step=100,label=FALSE)
rarecurve(seqtab.trim.rare,step=100,label=FALSE)

#phyloseq object but rarefied & trimmed
ps.trim.rare <- phyloseq(otu_table(seqtab.trim.rare, taxa_are_rows=FALSE), 
                         sample_data(samdf.trim), 
                         tax_table(taxa2))
ps.trim.rare #835 taxa and 120 samples

#saving
#### data files - rarefied, decontaminated, trimmed ####
#saving
write.csv(seqtab.trim.rare, file="CW_2020_16S_seqtab.trim.rare.3.5k.csv")
saveRDS(ps.trim.rare,file="CW_2020_16S_ps.trim.rare.RDS")
#write trim rare sample data frame too
samdf.trim.rare <- data.frame(ps.trim.rare@sam_data) 
write.csv(samdf.trim.rare, file="CW_2020_16S_samdf.trim.rare.csv")

###check through phyloseq objects to make sure they have all been saved successfully
#can extract sample data, taxa data, and sequence tables

#first ps.less
ps.less <- readRDS(file="CW_2020_16S_ps.less.RDS")
ps.less #27653 taxa and 135 samples
#ps.rare = rarefied to 10,000
ps.rare <- readRDS(file="CW_2020_16S_ps.rare.RDS")
ps.rare #23759 taxa and 98 samples
#ps.trim = trimmed using MCMC OTU
ps.trim <- readRDS(file="CW_2020_16S_ps.trim.RDS")
ps.trim #835 taxa and 134 samples
#ps.trim.rare = trimmed and then rarefied to 3498
ps.trim.rare <- readRDS(file="CW_2020_16S_ps.trim.rare.RDS")
ps.trim.rare #835 taxa and 120 samples
#all files above created and saved 3/28/2023 MEP
# check read depths for all objects
#ps.less
mean(sample_sums(ps.less)) #31,752.5
sd(sample_sums(ps.less)) #26,995.04
sum(sample_sums(ps.less)) #4,286,588
#ps.rare (nice yes logical good)
mean(sample_sums(ps.rare)) #10000
sd(sample_sums(ps.rare)) #0
sum(sample_sums(ps.rare)) #980000
#ps.trim
mean(sample_sums(ps.trim)) #23,373.57
sd(sample_sums(ps.trim)) #21,066.37
sum(sample_sums(ps.trim)) #3,132,059
#ps.trim.rare (good yes)
mean(sample_sums(ps.trim.rare)) #3498
sd(sample_sums(ps.trim.rare)) #0
sum(sample_sums(ps.trim.rare)) #419,760


#### making fasta filea for picrust2 - do 4 types - less, rare, trim, trim rare ####
library(phyloseq)
library(dada2)

seqtab.less <-read.csv("CW_2020_16S_seqtab.less.csv", header = FALSE)
#ps.less
less.otu <- as.matrix(ps.less@otu_table)
less.taxa <- data.frame(ps.less@tax_table)
rownames(less.taxa)==colnames(less.otu)
colnames(less.otu) <- less.taxa$V8
ids.less <- rownames(less.taxa)
path.less="~/Documents/Castillo Lab/CW_2020/CW_2020_16S/CW_2020_16S_less.fasta"
uniquesToFasta(less.otu, path.less, ids = ids.less, mode = "w", width = 20000)
#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.less.t <- t(as.data.frame(less.otu))
write.table(seqtab.less.t,file="CW_2020_16S_seqtab.less.t.txt")
#text to column and converted to tab delimited file in Excel

#ps.rare
rare.otu <- as.matrix(ps.rare@otu_table)
rare.taxa <- data.frame(ps.rare@tax_table)
rownames(rare.taxa)==colnames(rare.otu)
colnames(rare.otu) <- rare.taxa$V8
ids.rare <- rownames(rare.taxa)
path.rare="~/Documents/Castillo Lab/CW_2020/CW_2020_16S/CW_2020_16S_rare.fasta"
uniquesToFasta(rare.otu, path.rare, ids = ids.rare, mode = "w", width = 20000)
#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.rare.t <- t(seqtab.rare)
write.table(seqtab.rare.t,file="CW_2020_16S_seqtab.rare.t.txt")
#text to column and converted to tab delimited file in Excel

#ps.trim
trim.otu <- as.matrix(ps.trim@otu_table)
trim.taxa <- data.frame(ps.trim@tax_table)
rownames(trim.taxa)==colnames(trim.otu)
colnames(trim.otu) <- trim.taxa$V8
ids.trim <- rownames(trim.taxa)
path.trim="~/Documents/Castillo Lab/CW_2020/CW_2020_16S/CW_2020_16S_trim.fasta"
uniquesToFasta(trim.otu, path.trim, ids = ids.trim, mode = "w", width = 20000)
#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.trim.t <- t(seqtab.trim)
write.table(seqtab.trim.t,file="CW_2020_16S_seqtab.trim.t.txt")
#text to column and converted to tab delimited file in Excel

#ps.trim.rare
trim.rare.otu <- as.matrix(ps.trim.rare@otu_table)
trim.rare.taxa <- data.frame(ps.trim.rare@tax_table)
rownames(trim.rare.taxa)==colnames(trim.rare.otu)
colnames(trim.rare.otu) <- trim.rare.taxa$V8
ids.trim.rare <- rownames(trim.rare.taxa)
path.trim.rare="~/Documents/Castillo Lab/CW_2020/CW_2020_16S/CW_2020_16S_trim.rare.fasta"
uniquesToFasta(trim.rare.otu, path.trim.rare, ids = ids.trim.rare, mode = "w", width = 20000)
#re-formatting seq table so picrust likes it:
#a tab-delimited table with ASV ids as the first column and sample abundances as all subsequent columns
seqtab.trim.rare.t <- t(seqtab.trim.rare)
write.table(seqtab.trim.rare.t,file="CW_2020_16S_seqtab.trim.rare.t.txt")
#text to column and converted to tab delimited file in Excel

#Now move on to analysis of community composition and diversity!!
#all info in new scripts
