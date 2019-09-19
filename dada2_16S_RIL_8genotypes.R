library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

path <- "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria"
#
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


# Optimize truncLen -------------------------------------------------------
track_heads_1 <-c()
track_heads_2 <-c()
track_heads_3 <-c()

# first set tested (decreases in forward improved numbers)
# forward_lengths <- rep(seq(from = 240, to = 280, by=10),4)
# backward_lengths <- c(rep(220,5),rep(230,5),rep(240,5),rep(250,5))
# second set tested ()
# forward_lengths <- seq(from = 220, to = 250, by=10)
# backward_lengths <- c(rep(250,4))
# third set tested ()
# forward_lengths <- c(230,230,230)
# backward_lengths <- c(250,260,270)

# all_lengths_tested_f<-c(rep(seq(from = 240, to = 280, by=10),4),seq(from = 220, to = 250, by=10),c(230,230,230))
# all_lengths_tested_b<- c(rep(220,5),rep(230,5),rep(240,5),rep(250,5),rep(250,4),c(250,260,270))

# Next the value maxEE was optimized from maxEE=c(2,i), with i equaling 2:5
# The orginal parameter of 6 was optimal

for (i in 2:5) {

# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,6), truncLen=c(forward_lengths[i], backward_lengths[i]), truncQ=2,
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,i), truncLen=c(230, 250), truncQ=2,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# previously 17% removed
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

# plot of the distribution of sequence lengths
plot(sort(table(nchar(getSequences(seqtab)))), ylab = "count", xlab="length", main = "distribution of \nsequence lengths")


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

sum(seqtab.nochim)/sum(seqtab)

track_heads_1<-rbind(track_heads_1, as.numeric(track[1,]))
track_heads_2<-rbind(track_heads_2, as.numeric(track[2,]))
track_heads_3<-rbind(track_heads_3, as.numeric(track[3,]))

}

average_paired_ends_kept <- rowSums(cbind(track_heads_1[,6]/track_heads_1[,1],track_heads_2[,6]/track_heads_2[,1],track_heads_3[,6]/track_heads_3[,1]))/3
# see track_heads_1, track_heads_2, track_heads_3 to see results for the 3 sampels


# Save the original space tested. Test some more
track_heads_1_1 <- track_heads_1
track_heads_2_1 <- track_heads_2
track_heads_3_1 <- track_heads_3



#### Test to optimize another parameters




# Real data ---------------------------------------------------------------

# remove truncLen parameter
# vivi had the trimLeft=c(12,21)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=c(2,6), truncLen=c(270, 220), truncQ=2,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# previously 17% removed
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

# now I keep approximately 45% of reads! much more than the original ~30%
# Viviane has ~50% of her reads. She uses slightly different merging statistics, perhaps I will check those out.

# Download the silva database and transfer to server https://zenodo.org/record/1172783#.XSm6R5MzYWo
# scp  scp /Users/benoyserman/Downloads/silva_nr_v132_train_set.fa.gz BenO@bioinf.nioo.knaw.nl:/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria
# assign Taxonommy
taxa <- assignTaxonomy(seqtab.nochim, "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria/taxa/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

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

write.table(t(seqtab.nochim.bacteria), "~/RIL_results/dada2_outputs/seqtab-nochim-bacteria.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim.bacteria, fout='~/RIL_results/dada2_outputs/rep-seqs-bacteria.fna', ids=colnames(seqtab.nochim.bacteria))



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

write.table(t(seqtab.nochim.bacteria.asv), "~/RIL_results/dada2_outputs/seqtab-nochim-bacteria.asv.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

taxa.bacteria.asv <- taxa.bacteria
taxa.bacteria.asv <- cbind(paste0("ASV",seq(dim(taxa.bacteria)[1])),taxa.bacteria)
write.table(taxa.bacteria.asv, "~/RIL_results/dada2_outputs/taxa.bacteria.asv.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

T_A_table <- cbind(rownames(T_A_df),T_A_df) # row.names must be false for this table to fit into 
colnames(T_A_table)[1] <- "Samples"
rownames(T_A_table) <- NULL

T_A_table <- rbind(colnames(T_A_table),T_A_table)

write.table(T_A_table, "~/RIL_results/dada2_outputs/T_A_df.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

sphingo_ASVs_with_QTL_positions <- 
uniquesToFasta(seqtab.nochim.bacteria, fout='~/RIL_results/dada2_outputs/rep-seqs-bacteria.fna', ids=colnames(seqtab.nochim.bacteria))


#######################################################################################################
##### Effective Sample Size - Metagenomeseq - filter based on this to remove low abundant reads  #####
#######################################################################################################

####################################################
# First load all the data and make an MRexperiment #
####################################################

#Loading data
tmp = loadMeta("~/RIL_results/dada2_outputs/seqtab-nochim-bacteria.asv.txt", sep = "\t")
taxa = read.delim("~/RIL_results/dada2_outputs/taxa.bacteria.asv.txt", header=FALSE, row.names = 1, stringsAsFactors = FALSE) 
# Warning: phenotypes must have the same names as the columns on the count matrix when we create the MRexperiment object for provenance purposes
colnames(tmp$count)

mapfile = loadPhenoData("~/RIL_results/dada2_outputs/T_A_df.txt", tran = TRUE)  
phenotypeData = AnnotatedDataFrame(mapfile)
OTUdata = AnnotatedDataFrame(taxa)

RIL_obj = newMRexperiment(tmp$counts, phenoData = phenotypeData, featureData = OTUdata) 

##############################################################
##### Second, do the CSS normalization - Metagenomeseq #######
##############################################################

# Calculates the percentile for which to sum counts up to and scale by. 
p = cumNormStat(RIL_obj)
# Now that we have the percentile, we then calculate the normalizion factor
RIL_obj = cumNorm(RIL_obj, p = p)

# Now it is possible to export these counts, however first we will remove the effective samples
# mat = MRcounts(lungData, norm = TRUE, log = TRUE)[1:5, 1:5]

treatmentstatus = pData(RIL_obj)$Treatment
settings = zigControl(tol=1e-10, maxit = 25, verbose = TRUE)
mod = model.matrix(~treatmentstatus)
colnames(mod) = levels(treatmentstatus)
res = fitZig(obj = RIL_obj, mod = mod, useCSSoffset = TRUE, control = settings)
names(res)
effective_samples_res <- calculateEffectiveSamples(res)
samples_to_keep <- as.numeric(which(effective_samples_res>ave(effective_samples_res)[1]))
write.table(effective_samples_res, "~/RIL_results/dada2_outputs/effsamples_RIL.txt", sep= "\t")

eff_sample_minimum <- as.numeric(ave(res1)[1])
rareFeatures = which(rowSums(MRcounts(RIL_obj)) < eff_sample_minimum)
RIL_obj_trim = RIL_obj[-rareFeatures, ]

normFactor = normFactors(RIL_obj_trim)
normFactor = log2(normFactor/median(normFactor) + 1)

#Prepare normalized data to export

mat = MRcounts(obj_trim, norm = TRUE, log = TRUE)
mat

#Export normalized data
exportMat(mat, file = file.path("~/RIL_results/dada2_outputs/CSS_final.txt"))

# Convert matrix into phyloseq compatible matrix
tmat <- t(mat)
colnames(tmat) <- colnames(ps@otu_table)[1:dim(tmat)[2]]
ps_CSS <- ps
otu_table(ps_CSS) <- otu_table(tmat, taxa_are_rows = TRUE)

CSS_trimmed_taxonomy_table <- obj_trim@featureData@data
rownames(CSS_trimmed_taxonomy_table) <- CSS_trimmed_taxonomy_table[1,]
colnames(CSS_trimmed_taxonomy_table)<- colnames(taxa.bacteria)
tmat <- t(mat)
taxa_names(tmat) <- rownames(CSS_trimmed_taxonomy_table)

####################
##### Phyloseq #####
####################

# 9477 is based on the effective sample size

ps_CSS <- phyloseq(otu_table(tmat, taxa_are_rows=FALSE), 
               sample_data(T_A_df), 
               tax_table(taxa.bacteria)[1:9477])

# write.table(taxa.bacteria[1:9477,],file="/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/results_RIL2019/effective_ASVs_taxonomy.txt", row.names=c(paste0("ASV", seq(ntaxa(ps_all)))))

ps <- phyloseq(otu_table(seqtab.nochim.bacteria, taxa_are_rows=FALSE), 
               sample_data(T_A_df), 
               tax_table(taxa.bacteria))

ps.trimmed <- phyloseq(otu_table(seqtab.nochim.trimmed, taxa_are_rows=FALSE), 
                       sample_data(T_A_df), 
                       tax_table(taxa.trimmed))

# Bulk, Money and Pimpi
non_RIL <- order(RIL_treatments_long)[193:225]
ps_noRIL <- prune_samples(sample_names(ps)[non_RIL], ps_CSS)
dna_ps_noRIL <- Biostrings::DNAStringSet(taxa_names(ps_noRIL))
names(dna_ps_noRIL) <- taxa_names(ps_noRIL)
ps_noRIL <- merge_phyloseq(ps_noRIL, dna_ps_noRIL)
taxa_names(ps_noRIL) <- paste0("ASV", seq(ntaxa(ps_noRIL)))
ps_noRIL
plot_richness(ps_noRIL, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps.prop_noRIL <- transform_sample_counts(ps_noRIL, function(otu) otu/sum(otu))
ord.nmds.bray_noRIL <- ordinate(ps.prop_noRIL, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps.prop_noRIL, ord.nmds.bray_noRIL, color="Treatment", title="NMDS",label="Accession")
ord.PCoA_noRIL <- ordinate(ps.prop_noRIL, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps.prop_noRIL, ord.PCoA_noRIL, color="Treatment", title="Bray PCoA",label="Accession")


# All samples
ps_all <- prune_samples(sample_names(ps_CSS)[c(1:225)[-low_coverage_samples]], ps_CSS)
dna_ps_all <- Biostrings::DNAStringSet(taxa_names(ps_all))
names(dna_ps_all) <- taxa_names(ps_all)
ps_all <- merge_phyloseq(ps_all, dna_ps_all)
taxa_names(ps_all) <- paste0("ASV", seq(ntaxa(ps_all)))
ps_all
plot_richness(ps_all, x="Treatment", measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"), color="Treatment")
ps.prop_all <- transform_sample_counts(ps_all, function(otu) otu/sum(otu))
ord.nmds.bray_all <- ordinate(ps.prop_all, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps.prop_all, ord.nmds.bray_all, color="Treatment", title="NMDS",label="Accession")
ord.PCoA_all <- ordinate(ps.prop_all, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps.prop_all, ord.PCoA_all, color="Treatment", title="Bray PCoA",label="Accession")

Biostrings::writeXStringSet("~/RIL_results/", ps_all)





# try it at family level
ps_all_fam_taxglom <- tax_glom(ps_all, taxrank=rank_names(ps_all)[5], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
plot_richness(ps_all_fam_taxglom, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps_all_fam_taxglom_prop <- transform_sample_counts(ps_all_fam_taxglom, function(otu) otu/sum(otu))
ord.nmds.bray_genus_taxglom <- ordinate(ps_all_fam_taxglom_prop, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps_all_fam_taxglom_prop, ord.nmds.bray_genus_taxglom, color="Treatment", title="NMDS",label="Accession")
ord.PCoA_fam_taxglom <- ordinate(ps_all_fam_taxglom_prop, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps_all_fam_taxglom_prop, ord.PCoA_fam_taxglom, color="Treatment", title="Bray PCoA",label="Accession")

# try it at genus level
ps_all_genus_taxglom <- tax_glom(ps_all, taxrank=rank_names(ps_all)[6], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
plot_richness(ps_all_genus_taxglom, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps_all_genus_taxglom_prop <- transform_sample_counts(ps_all_genus_taxglom, function(otu) otu/sum(otu))
ord.nmds.bray_genus_taxglom <- ordinate(ps_all_genus_taxglom_prop, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps_all_genus_taxglom, ord.nmds.bray_genus_taxglom, color="Treatment", title="NMDS",label="Accession")
ord.PCoA_genus_taxglom <- ordinate(ps_all_genus_taxglom_prop, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps_all_genus_taxglom_prop, ord.PCoA_genus_taxglom, color="Treatment", title="Bray PCoA",label="Accession")




# unifrac_all <- UniFrac(ps_all, TRUE)
# Just the Rhizospheres
ps_rhizosphere <- prune_samples(sample_names(ps)[c(1:225)[c(-low_coverage_samples,-order(RIL_treatments_long)[193:213])]], ps)
dna_ps_rhizosphere <- Biostrings::DNAStringSet(taxa_names(ps_rhizosphere))
names(dna_ps_rhizosphere) <- taxa_names(ps_rhizosphere)
ps_rhizosphere <- merge_phyloseq(ps_rhizosphere, dna_ps_rhizosphere)
taxa_names(ps_rhizosphere) <- paste0("ASV", seq(ntaxa(ps_rhizosphere)))
ps_rhizosphere
plot_richness(ps_rhizosphere, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps.prop_rhizosphere <- transform_sample_counts(ps_rhizosphere, function(otu) otu/sum(otu))
ord.nmds.bray_rhizosphere <- ordinate(ps.prop_rhizosphere, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps.prop_rhizosphere, ord.nmds.bray_rhizosphere, color="Replicate", title="NMDS",label="Accession")
ord.PCoA_rhizosphere <- ordinate(ps.prop_rhizosphere, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps.prop_rhizosphere, ord.PCoA_rhizosphere, color="Replicate", title="Bray PCoA Rhizospheres",label="Accession")
# try weighted and unweighted


ps.trimmed <- prune_samples(sample_names(ps.trimmed)[1:225], ps.trimmed)
dna <- Biostrings::DNAStringSet(taxa_names(ps.trimmed))
names(dna) <- taxa_names(ps.trimmed)
ps.trimmed <- merge_phyloseq(ps.trimmed, dna)
taxa_names(ps.trimmed) <- paste0("ASV", seq(ntaxa(ps.trimmed)))
ps.trimmed
plot_richness(ps.trimmed, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps.trimmed.prop <- transform_sample_counts(ps.trimmed, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.trimmed.prop, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps.trimmed.prop, ord.nmds.bray, color="Treatment", title="NMDS",label="Accession")
ord.PCoA <- ordinate(ps.trimmed.prop, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps.trimmed.prop, ord.PCoA, color="Treatment", title="Bray PCoA",label="Accession")

