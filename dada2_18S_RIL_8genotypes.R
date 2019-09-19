library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())

path_18S <- "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/eukaryota/"

# Notes after preliminary analysis:
# path <- "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/bacteria/test"
# Most sequences are between 440 and 465 basepairs. + 12 basepairs for overlap, so 477 basepairs needed
# after optimization, the ideal tuncLen setting is truncLen=c(230, 250)
# it looks like two samples were switched - a bulk and pimpi sample (samples 156 and 155)

# Viviane uses these parameters: trimLeft= c(17,21)

#####################################
############### Dada2 ###############
#####################################

list.files(path_18S)
fnFs_18S <- sort(list.files(path_18S, pattern="R1_001", full.names = TRUE))
fnRs_18S <- sort(list.files(path_18S, pattern="R2_001", full.names = TRUE))

##########################################################
### Identify the 'basename' from the sample file names ###
##########################################################

sample.names2_18S <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:3)
sample.names1_18S <- do.call(paste, as.data.frame(t(sample.names2),sep="_"))
sample.names_18S <- sub(" ","_",sample.names1)
sample.names_18S <- sub(" ","_",sample.names)
switched_names_18S <- c("35_P6_16S","36_Bulk3_2")
sample.names_18S[c(155,156)] <- switched_names_18S[c(2,1)]
  


filtFs_18S <- file.path(path_18S, "filtered", paste0(sample.names_18S, "_F_filt.fastq.gz"))
filtRs_18S <- file.path(path_18S, "filtered", paste0(sample.names_18S, "_R_filt.fastq.gz"))
names(filtFs_18S) <- sample.names_18S
names(filtRs_18S) <- sample.names_18S

plotQualityProfile(fnFs_18S[1:3])
plotQualityProfile(fnRs_18S[1:3])

# remove truncLen parameter
# vivi had the trimLeft=c(12,21)
out_18S <- filterAndTrim(fnFs_18S, filtFs_18S, fnRs_18S, filtRs_18S, maxEE=c(2,6), truncLen=c(295, 260), truncQ=2,
                     maxN=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out_18S)

# 8% removed
nrr_18S <- round(100 * (sum(out_18S[,1])-sum(out_18S[,2])) / sum(out_18S[,1]))
boxplot(out_18S, ylim=c(0,max(out_18S[,1])+10000), main = paste("approximately",nrr_18S,"% of reads removed"))


errF_18S <- learnErrors(filtFs_18S, multithread=TRUE)
errR_18S <- learnErrors(filtRs_18S, multithread=TRUE)

dadaFs_18S <- dada(filtFs_18S, err=errF_18S, multithread=TRUE)
dadaRs_18S <- dada(filtRs_18S, err=errR_18S, multithread=TRUE)

dadaFs_18S[[1]]
dadaRs_18S[[1]]

# almost no sequences were merged using 270 / 220 so it was changed to 270 and 250
# almost no sequences were merged using 270 / 250 which equals 520, the expected size is 516 ...
# almost no sequences were merged using 270 / 270 which equals 540, the expected size is 516 ...
# better using 270 / 270 which equals 540, the expected size is 516 ... so  changed to 280/270

mergers_18S <- mergePairs(dadaFs_18S, filtFs_18S, dadaRs_18S, filtRs_18S, verbose=TRUE)
head(mergers_18S)

seqtab_18S <- makeSequenceTable(mergers_18S)
dim(seqtab_18S)

# plot of the distribution of sequence lengths that have been merged
plot(sort(table(nchar(getSequences(seqtab_18S)))), ylab = "count", xlab="length", main = "distribution of \nsequence lengths")

# remove chimeras, using method consensus as this is more sensititve in big datasets
seqtab.nochim_18S <- removeBimeraDenovo(seqtab_18S, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim_18S)

# Calculate some statistics about how many were reads were removed at each step
sum(seqtab.nochim_18S)/sum(seqtab_18S)

getN <- function(x) sum(getUniques(x))
track_18S <- cbind(out_18S, sapply(dadaFs_18S, getN), sapply(dadaRs_18S, getN), sapply(mergers_18S, getN), rowSums(seqtab.nochim_18S))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_18S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_18S) <- sample.names_18S
head(track_18S)

sum(seqtab.nochim_18S)/sum(seqtab_18S)
average_paired_ends_kept_real_18S <- sum(track_18S[,6]/track_18S[,1])/dim(track_18S)[1]

# now I keep approximately 45% of reads! much more than the original ~30%
# Viviane has ~50% of her reads. She uses slightly different merging statistics, perhaps I will check those out.

# Download the silva database and transfer to server https://zenodo.org/record/1447330#.XWwYjZMzYWo
# scp /Users/benoyserman/Downloads/silva_132.18s.99_rep_set.dada2.fa.gz BenO@bioinf.nioo.knaw.nl:/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/eukaryota
# assign Taxonommy
taxa_18S <- assignTaxonomy(seqtab.nochim_18S, "/data/ngs/ME/raaijmakers_group/RIL2019/RIL2019_raw_16S_18S/raw_sequences/unzipped/eukaryota/silva_132.18s.99_rep_set.dada2.fa.gz", multithread=20)

# remove archaea, chloroplast, mitochondria, bacteria - none were found
pos_chloroplast_mitchochondria_18S <- c(which(taxa_18S[,4] %in% c("Chloroplast","Mitochondria")),which(taxa_18S[,5] %in% c("Chloroplast","Mitochondria")),which(taxa_18S[,1] %in% c("Archaea","Bacteria","NA")))
seqtab.nochim.eukaryota_18S <- seqtab.nochim_18S
taxa.bacteria_18S <- taxa_18S

# Now code for the various factors of treatment and possible confounding batch effects

samples.out_18S <- rownames(seqtab.nochim_18S)
samples.eukaryota.out_18S <- rownames(seqtab.nochim.eukaryota_18S)

seqs_18S <- getSequences(seqtab.nochim.eukaryota_18S)
names(seqs_18S) <- rownames(seqtab.nochim.eukaryota_18S) 

# It looks like a P and Bulk labels were mixed up
mixed_up_samples_18S <- c(which(samples.out_18S=="36_Bulk3_2"), which(samples.out_18S=="35_P6_16S"))
samples.out_18S[mixed_up_samples_18S] <- c("35_P6_16S","36_Bulk3_2")
rownames(seqtab.nochim.eukaryota_18S) <- samples.out_18S


write.table(t(seqtab.nochim.eukaryota_18S), "~/RIL_results/dada2_outputs/seqtab-nochim-eukaryota.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim.eukaryota_18S, fout='~/RIL_results/dada2_outputs/rep-seqs-eukaryota.fna', ids=colnames(seqtab.nochim.eukaryota_18S))



#####################
# Need to continue from here!
#####################

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




S_richness <- estimate_richness(ps.trimmed, split = TRUE, measures = "Shannon")

plot_richness(ps.trimmed.nobulk, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")

ps.prop.nobulk <- transform_sample_counts(ps.trimmed.nobulk, function(otu) otu/sum(otu))


# Treatment Type (Bulk, RIL, Money, Pimp)
Type <- c(rownames(seqtab.nochim))
# Treatment (Bulk, RIL #, Money, Pimp)
Treatment <- c(rownames(seqtab.nochim))
# Harvest Date
Harvest <- c(rownames(seqtab.nochim))
# Extraction Date
Extraction <- c(rownames(seqtab.nochim))
# Replicate (Y, P or (C)ontrol for bulk and parentals)
Replicate <- c(rownames(seqtab.nochim))



