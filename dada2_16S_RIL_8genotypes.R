library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

path <- "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_raw/RIL2019_raw_16S_18S/unzipped/bacteria"

# Notes after preliminary analysis:
# path <- "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria/test"
# Most sequences are between 440 and 465 basepairs. + 12 basepairs for overlap, so 477 basepairs needed
# after optimization, the ideal tuncLen setting is truncLen=c(230, 250)
# it looks like two samples were switched - a bulk and pimpi sample (samples 156 and 155)

# Viviane uses these parameters: trimLeft= c(17,21)

#####################################
############### Dada2 ###############
#####################################

list.files(path)
fnFs <- sort(list.files(path, pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001", full.names = TRUE))

##########################################################
### Identify the 'basename' from the sample file names ###
##########################################################

sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:3)
sample.names1 <- do.call(paste, as.data.frame(t(sample.names2),sep="_"))
sample.names <- sub(" ","_",sample.names1)
sample.names <- sub(" ","_",sample.names)
switched_names <- c("35_P6_16S","36_Bulk3_2")
sample.names[c(155,156)] <- switched_names[c(2,1)]

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,6), truncLen=c(270, 220), truncQ=2,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

nrr <- round(100 * (sum(out[,1])-sum(out[,2])) / sum(out[,1]))
boxplot(out, ylim=c(0,max(out[,1])+10000), main = paste("approximately",nrr,"% of reads removed"))


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# plot of the distribution of sequence lengths that have been merged
plot(sort(table(nchar(getSequences(seqtab)))), ylab = "count", xlab="length", main = "distribution of \nsequence lengths")

# remove chimeras, using method consensus as this is more sensititve in big datasets
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Calculate some statistics about how many were reads were removed at each step
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

sum(seqtab.nochim)/sum(seqtab)
average_paired_ends_kept_real <- sum(track[,6]/track[,1])/dim(track)[1]

# approximately 45% of reads are mantained

# Download the silva database and transfer to server https://zenodo.org/record/1172783#.XSm6R5MzYWo
# scp  scp /Users/benoyserman/Downloads/silva_nr_v132_train_set.fa.gz BenO@bioinf.nioo.knaw.nl:/data/ngs/ME/raaijmakers_group/RIL2019_raw/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria
# assign Taxonommy
taxa <- assignTaxonomy(seqtab.nochim, "/data/ngs/ME/raaijmakers_group/RIL2019_raw/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria/taxa/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# remove archaea, chloroplast, mitochondria, eukaryota
pos_chloroplast_mitchochondria <- c(which(taxa[,4] %in% c("Chloroplast","Mitochondria")),which(taxa[,5] %in% c("Chloroplast","Mitochondria")),which(taxa[,1] %in% c("Archaea","Eukaryota","NA")))
seqtab.nochim.bacteria <- seqtab.nochim[,-pos_chloroplast_mitchochondria]
taxa.bacteria <- taxa[-pos_chloroplast_mitchochondria,]

# Now code for the various factors of treatment and possible confounding batch effects

samples.out <- rownames(seqtab.nochim)
samples.out.trimmed <- rownames(seqtab.nochim.trimmed)
samples.bacteria.out <- rownames(seqtab.nochim.bacteria)

seqs <- getSequences(seqtab.nochim.bacteria)
names(seqs) <- rownames(seqtab.nochim.bacteria) 

write.table(t(seqtab.nochim.bacteria), "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/seqtab-nochim-bacteria_RIL_2019.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim.bacteria, fout='/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/rep-seqs-bacteria_RIL_2019.txt.fna', ids=colnames(seqtab.nochim.bacteria))



# It looks like a P and Bulk were mixed up
# mixed_up_samples <- c(which(samples.out.trimmed=="36_Bulk3_2"), which(samples.out.trimmed=="35_P6_16S"))
# samples.out.trimmed[c(156,155)] <- c("35_P6_16S","36_Bulk3_2")
# samples.out.trimmed[mixed_up_samples] <- c("35_P6_16S","36_Bulk3_2")
# rownames(seqtab.nochim) <- samples.out.trimmed



#################################################################################
##########         Create Meta Data from the names of the files        ##########
#################################################################################

Replicate <- sapply(strsplit(samples.bacteria.out[1:225], "_"), `[`, 3)
Replicate[order(RIL_treatments_long)[193:213]] <- c("Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk","Bulk")
Replicate[order(RIL_treatments_long)[214:219]] <- c("M","M","M","M","M","M")
Replicate[order(RIL_treatments_long)[220:225]] <- c("Pimpi","Pimpi","Pimpi","Pimpi","Pimpi","Pimpi")
Replicate <- toupper(Replicate)
Replicate[226:261] <- "8genotype"

RIL_treatments_long <- sapply(strsplit(samples.bacteria.out[1:225], "_"), `[`, 2)
order_of_RIL_treatments <- order(RIL_treatments_long)
RIL_treatment_short <- RIL_treatments_long
RIL_treatment_short[order(RIL_treatments_long)[1:192]] <- "RIL"
RIL_treatment_short[order(RIL_treatments_long)[193:213]] <- "BULK"
RIL_treatment_short[order(RIL_treatments_long)[214:219]] <- "M"
RIL_treatment_short[order(RIL_treatments_long)[220:225]] <- "P"

RIL_treatment_medium<-RIL_treatments_long
RIL_treatment_medium[order(RIL_treatments_long)[193:213]] <- c("Bulk1","Bulk1","Bulk1","Bulk1","Bulk2","Bulk2","Bulk2","Bulk2","Bulk3","Bulk3","Bulk3","Bulk3","Bulk4","Bulk4","Bulk4","Bulk4","Bulk5","Bulk5","Bulk5","Bulk7","Bulk7")
RIL_treatment_medium[order(RIL_treatments_long)[214:219]] <- c("M1","M2","M3","M4","M5","M6")
RIL_treatment_medium[order(RIL_treatments_long)[220:225]] <- c("P1","P2","P3","P4","P5","P6")

Eight_genotpe_treatments <- sapply(strsplit(samples.out[226:length(samples.bacteria.out)], "_"), `[`, 1)
All_treatment_types_short <- c(RIL_treatment_short,Eight_genotpe_treatments)
All_treatment_types_medium <- c(RIL_treatment_medium,Eight_genotpe_treatments)

RIL_accession_replicate <- cbind(RIL_treatments_long,Replicate[1:225])

# All_treatment_types_short[which(samples.out=="36_Bulk3_2")] <- "P5"
# All_treatment_types_short[which(samples.out=="169_P5_16S")] <- "Bulk3"
# All_treatment_types_medium[which(samples.out=="36_Bulk3_2")] <- "P5"
# All_treatment_types_medium[which(samples.out=="169_P5_16S")] <- "Bulk3"
# All_treatment_types[which(samples.out=="36_Bulk3_2")] <- "P"
# All_treatment_types[which(samples.out=="169_P5_16S")] <- "BULK"

All_Accessions <- All_treatment_types_medium
T_A_df <- data.frame(Treatment=All_treatment_types_short, Accession=All_Accessions, Replicate = Replicate)
rownames(T_A_df) <- samples.bacteria.out
# rownames(seqtab.nochim.trimmed) <- samples.out.trimmed

# must remove 52_228_P and 190_254_Y for low coverage. 148_256_Y because it is a big outlier
low_coverage_samples <- which(samples.out.trimmed %in% c("52_228_P","190_254_Y","148_256_Y"))




################################################
############ Make table for Metagenomeseq ######
################################################

seqtab.nochim.bacteria.asv <- seqtab.nochim.bacteria
seqtab.nochim.bacteria.asv <- rbind(paste0("ASV",seq(dim(seqtab.nochim.bacteria.asv)[2])),seqtab.nochim.bacteria.asv)
seqtab.nochim.bacteria.asv <- cbind(rownames(seqtab.nochim.bacteria.asv),seqtab.nochim.bacteria.asv)
colnames(seqtab.nochim.bacteria.asv)<- NULL
rownames(seqtab.nochim.bacteria.asv)<- NULL

write.table(t(seqtab.nochim.bacteria.asv), "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/seqtab-nochim-bacteria.asv_RIL_and_8genotypes_2019.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

taxa.bacteria.asv <- taxa.bacteria
taxa.bacteria.asv <- cbind(paste0("ASV",seq(dim(taxa.bacteria)[1])),taxa.bacteria)
write.table(taxa.bacteria.asv, "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/taxa.bacteria.asv_RIL_and_8genotypes_2019.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

T_A_table <- cbind(rownames(T_A_df),T_A_df) # row.names must be false for this table to fit into 
colnames(T_A_table)[1] <- "Samples"
rownames(T_A_table) <- NULL

T_A_table <- rbind(colnames(T_A_table),T_A_table)

write.table(T_A_table, "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/Treatment_Accession_replicate_table_RIL_and_8genotypes_2019.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



