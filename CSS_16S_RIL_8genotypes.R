#######################################################################################################
##### Effective Sample Size - Metagenomeseq - filter based on this to remove low abundant reads  #####
#######################################################################################################

####################################################
# First load all the data and make an MRexperiment #
####################################################

#Loading data
tmp = loadMeta("/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/seqtab-nochim-bacteria.asv_RIL_and_8genotypes_2019.txt", sep = "\t")
taxa = read.delim("/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/taxa.bacteria.asv_RIL_and_8genotypes_2019.txt", header=FALSE, row.names = 1, stringsAsFactors = FALSE) 
# Warning: phenotypes must have the same names as the columns on the count matrix when we create the MRexperiment object for provenance purposes
colnames(tmp$count)

mapfile = loadPhenoData("/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/Treatment_Accession_replicate_table_RIL_and_8genotypes_2019.txt", tran = TRUE)  
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
write.table(effective_samples_res, "/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/effsamples_RIL_and_8genotypes_2019.txt", sep= "\t")

eff_sample_minimum <- as.numeric(ave(res1)[1])
rareFeatures = which(rowSums(MRcounts(RIL_obj)) < eff_sample_minimum)
RIL_obj_trim = RIL_obj[-rareFeatures, ]

normFactor = normFactors(RIL_obj_trim)
normFactor = log2(normFactor/median(normFactor) + 1)

#Prepare normalized data to export

mat = MRcounts(obj_trim, norm = TRUE, log = TRUE)
mat

#Export normalized data
exportMat(mat, file = file.path("/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/RIL2019_results/16S_processing/CSS_final_RIL_and_8genotypes_2019.txt"))

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

# Bulk, Money and Pimpi
non_RIL <- order(RIL_treatments_long)[193:225]
ps_noRIL <- prune_samples(sample_names(ps_CSS)[non_RIL], ps_CSS)
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


# All samples # Except for the 8 genotypes
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

# Just the Rhizospheres
ps_rhizosphere <- prune_samples(sample_names(ps_all)[c(1:225)[c(-low_coverage_samples,-order(RIL_treatments_long)[193:213])]], ps_all)
dna_ps_rhizosphere <- Biostrings::DNAStringSet(taxa_names(ps_rhizosphere))
names(dna_ps_rhizosphere) <- taxa_names(ps_rhizosphere)
ps_rhizosphere <- merge_phyloseq(ps_rhizosphere, dna_ps_rhizosphere)
taxa_names(ps_rhizosphere) <- paste0("ASV", seq(ntaxa(ps_rhizosphere)))
ps_rhizosphere
#plot_richness(ps_rhizosphere, x="Treatment", measures=c("Shannon", "Simpson"), color="Treatment")
ps.prop_rhizosphere <- transform_sample_counts(ps_rhizosphere, function(otu) otu/sum(otu))
ord.nmds.bray_rhizosphere <- ordinate(ps.prop_rhizosphere, method="NMDS", distance="bray", trymax = 50)
plot_ordination(ps.prop_rhizosphere, ord.nmds.bray_rhizosphere, color="Replicate", title="NMDS",label="Accession")
ord.PCoA_rhizosphere <- ordinate(ps.prop_rhizosphere, method="PCoA", distance="bray", trymax = 50)
plot_ordination(ps.prop_rhizosphere, ord.PCoA_rhizosphere, color="Replicate", title="Bray PCoA Rhizospheres",label="Accession")
# try weighted and unweighted