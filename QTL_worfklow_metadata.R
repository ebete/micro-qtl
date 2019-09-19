#!/usr/bin/env Rscript
options(error=function()traceback(2))
#args = commandArgs(trailingOnly=TRUE)
#inp_folder <- args[[1]]
#prefix <- args[[2]]
inp_folder <- "/mnt/nfs/bioinfdata/home/NIOO/beno/microbiomeQTL/microbiomeQTL/data"
setwd("/mnt/nfs/bioinfdata/home/NIOO/beno/microbiomeQTL/microbiomeQTL/data")

start.time <- Sys.time()

library(devtools)
# install.packages("qtl", repos = "http://cran.us.r-project.org")
# install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
library(qtl)
library(qtl2)
library(parallel)
library(MASS)
library(base)
#install.packages('RColorBrewer')
library(RColorBrewer)
library(plyr)

##########################################################################################################################################
################################################# First put the marker table into the correct format######################################
##########################################################################################################################################

# Restructure the genotype file
marker_table <- read.table(file.path(inp_folder, "markers_new_map.csv"), sep = ",", comment.char="", header=F, check.names=F, stringsAsFactors=F)
marker_table[1,1] <- "id"
marker_table[2,1] <- ""
marker_table[3,1] <- ""
marker_table <- marker_table[-4,]
marker_table[marker_table == "BB"] <- "B"
marker_table[marker_table == "AA"] <- "A"


# Remove samples that are missing in either the genotype or the phenotype using the biomass meta data
meta_data_RIL <- read.table(file= "/mnt/nfs/bioinfdata/home/NIOO/beno/microbiomeQTL/microbiomeQTL/RIL_metadata_5_19_2019.txt", header=TRUE)

# These are the numeric names of all the RILs
RIL_id_marker_table <- marker_table$V1[-c(1,2,3)]
RIL_id_marker_table

####################################################################################################
############# Phenotypes: Section 1) Plant Dry Weight 2) ASV 3) eASV 4) Genus level ################
####################################################################################################

######################################################
############# Phenotypes: Section 1)  ################
######################################################

# The first phenotype we test is the plant dry weight. This is to be compared with two other investigations that have had PDW as a phenotype
# Plant dry weight of the entire population # 
 RIL_id <- c("id", as.character(meta_data_RIL$RIL[c(1:93,95)]))
 PDW <- c("PDW",as.numeric(meta_data_RIL$PDW[c(1:93,95)]))
 PDW2 <- c("PDW2",as.numeric(meta_data_RIL$PDW[c(97:189,191)]))

# Re-order the data into the same order as in the marker id table
 RIL_id_ordered<- c("id",RIL_id[match(RIL_id_marker_table,RIL_id)])
 PDW_ordered <- c("PDW_P",PDW[match(RIL_id_marker_table,RIL_id)])
 PDW2_ordered <- c("PDW_Y",PDW2[match(RIL_id_marker_table,RIL_id)])
 
 # build a single table with all the PDW
 phenotype_PDW <- as.matrix(cbind(RIL_id_ordered,PDW_ordered,PDW2_ordered)) 
 colnames(phenotype_PDW)<-NULL
 phenotype_PDW<-phenotype_PDW[which(!is.na(phenotype_PDW[,1])),] 
 
 # These are removed because of the date the plants were harvested. Harest was later because germination was much later.
 # The number of days they grew was variable, they were only harvested out of curiousity, and as practice samples.

 Y_removed <- c(253,226,259,275)
 P_removed <- c(217,237,249)
 
 phenotype_PDW[which(phenotype_PDW[,1] %in% P_removed),2] <- NA
 phenotype_PDW[which(phenotype_PDW[,1] %in% Y_removed),3] <- NA
 
 # Add the average as a final phenotype
 PDW_average_ordered <- (as.numeric(phenotype_PDW[2:95,2])+as.numeric(phenotype_PDW[2:95,3]))/2
 phenotype_PDW <- cbind(phenotype_PDW,NA)
 phenotype_PDW[,4] <- c("PDW_average",PDW_average_ordered)
 
 # export the final phenotype file
 phenotype_PDW_file <- file.path(inp_folder, sprintf("phenotype_PDW_file.txt"))
 write.table(phenotype_PDW, file=phenotype_PDW_file, sep=",", quote=F, col.names=F, row.names=F)
 
 ##############################################################
 ################# Phenotypes: Section 1)  ####################
 ##############################################################
 # now for the microbioas as a phenotypes
 # first with ASV

 
 P_replicate_ASV_matrix <- ps_all@otu_table[which(T_A_table$Replicate[match(rownames(ps_all@otu_table),T_A_table$Samples)]=="P"),]
 Y_replicate_ASV_matrix <- ps_all@otu_table[which(T_A_table$Replicate[match(rownames(ps_all@otu_table),T_A_table$Samples)]=="Y"),]
 
 
 RIL_P_ASV <- sapply(strsplit(rownames(P_replicate_ASV_matrix), "_"), `[`, 2)
 RIL_Y_ASV <- sapply(strsplit(rownames(Y_replicate_ASV_matrix), "_"), `[`, 2)
 
 rownames(P_replicate_ASV_matrix) <- RIL_P_ASV
 rownames(Y_replicate_ASV_matrix) <- RIL_Y_ASV
 
 phenotype_ASV_P <- cbind(rownames(P_replicate_ASV_matrix),P_replicate_ASV_matrix)
 phenotype_ASV_P <- rbind(colnames(phenotype_ASV_P),phenotype_ASV_P)
 phenotype_ASV_P[1,1]<- "id"
 
 phenotype_ASV_Y <- cbind(rownames(Y_replicate_ASV_matrix),Y_replicate_ASV_matrix)
 phenotype_ASV_Y <- rbind(colnames(phenotype_ASV_Y),phenotype_ASV_Y)
 phenotype_ASV_Y[1,1]<- "id"
 
 
 colnames(phenotype_cASV_P) <- NULL
 rownames(phenotype_cASV_P) <- NULL
 
 colnames(phenotype_cASV_Y) <- NULL
 rownames(phenotype_cASV_Y) <- NULL
 
 phenotype_ASV_P_ordered <- phenotype_ASV_P[match(c("id",RIL_id_marker_table),phenotype_ASV_P[,1]),]
 phenotype_ASV_Y_ordered <- phenotype_ASV_Y[match(c("id",RIL_id_marker_table),phenotype_ASV_Y[,1]),]
 
 # Remove outlying and replicates that grew much longer
 phenotype_ASV_P_ordered_filtered <- phenotype_ASV_P_ordered
 phenotype_ASV_Y_ordered_filtered <- phenotype_ASV_Y_ordered
 phenotype_ASV_P_ordered_filtered[which(phenotype_ASV_P_ordered[,1] %in% P_removed),2:dim(phenotype_ASV_P_ordered)[2]] <- NA
 phenotype_ASV_Y_ordered_filtered[which(phenotype_ASV_Y_ordered[,1] %in% Y_removed),2:dim(phenotype_ASV_Y_ordered)[2]] <- NA
 
 
 # Remove bulk and parentals
 
 phenotype_ASV_P_ordered_filtered <- phenotype_ASV_P_ordered_filtered[!is.na(phenotype_ASV_P_ordered_filtered[,1]),]
 phenotype_ASV_Y_ordered_filtered <- phenotype_ASV_Y_ordered_filtered[!is.na(phenotype_ASV_Y_ordered_filtered[,1]),]

 phenotype_ASV_P_file <- file.path(inp_folder, sprintf("phenotype_ASV_P_file.txt"))
 phenotype_ASV_Y_file <- file.path(inp_folder, sprintf("phenotype_ASV_Y_file.txt"))

 write.table(phenotype_ASV_P_ordered_filtered, file=phenotype_ASV_P_file, sep=",", quote=F, col.names=F, row.names=F)
 write.table(phenotype_ASV_Y_ordered_filtered, file=phenotype_ASV_Y_file, sep=",", quote=F, col.names=F, row.names=F)
 
 ##############################################################
 ################# Phenotypes: Section 2)  eASV ###############
 ##############################################################
 
 # now with the cASVs
 
 P_replicate_c_ASV_matrix <- c_ASV_matrix[which(T_A_table$Replicate[match(rownames(c_ASV_matrix),T_A_table$Samples)]=="P"),]
 Y_replicate_c_ASV_matrix <- c_ASV_matrix[which(T_A_table$Replicate[match(rownames(c_ASV_matrix),T_A_table$Samples)]=="Y"),]
 
 # modify the matrix to fit the standard format with no colnames/rownames, but instead as the first row and column. "id" must be in first position
 
 RIL_P <- sapply(strsplit(rownames(P_replicate_c_ASV_matrix), "_"), `[`, 2)
 RIL_Y <- sapply(strsplit(rownames(Y_replicate_c_ASV_matrix), "_"), `[`, 2)
 
 rownames(P_replicate_c_ASV_matrix) <- RIL_P
 rownames(Y_replicate_c_ASV_matrix) <- RIL_Y
 
 phenotype_cASV_P <- cbind(rownames(P_replicate_c_ASV_matrix),P_replicate_c_ASV_matrix)
 phenotype_cASV_P <- rbind(colnames(phenotype_cASV),phenotype_cASV)
 phenotype_cASV_P[1,1]<- "id"
 
 phenotype_cASV_Y <- cbind(rownames(Y_replicate_c_ASV_matrix),Y_replicate_c_ASV_matrix)
 phenotype_cASV_Y <- rbind(colnames(phenotype_cASV_Y),phenotype_cASV_Y)
 phenotype_cASV_Y[1,1]<- "id"
 
 
 colnames(phenotype_cASV_P) <- NULL
 rownames(phenotype_cASV_P) <- NULL
 
 colnames(phenotype_cASV_Y) <- NULL
 rownames(phenotype_cASV_Y) <- NULL
 
 phenotype_cASV_P_ordered <- phenotype_cASV_P[match(c("id",RIL_id_marker_table),phenotype_cASV_P[,1]),]
 phenotype_cASV_Y_ordered <- phenotype_cASV_Y[match(c("id",RIL_id_marker_table),phenotype_cASV_Y[,1]),]
 
 # Remove outlying and replicates that grew much longer
 phenotype_cASV_P_ordered_filtered <- phenotype_cASV_P_ordered
 phenotype_cASV_Y_ordered_filtered <- phenotype_cASV_Y_ordered
 phenotype_cASV_P_ordered_filtered[which(phenotype_cASV_P_ordered[,1] %in% P_removed),2:dim(phenotype_cASV_P_ordered)[2]] <- NA
 phenotype_cASV_Y_ordered_filtered[which(phenotype_cASV_Y_ordered[,1] %in% Y_removed),2:dim(phenotype_cASV_Y_ordered)[2]] <- NA
 
 # Take averages of these
 phenotype_cASV_P_ordered_numeric_matrix <- matrix(as.numeric(phenotype_cASV_P_ordered_filtered[2:dim(phenotype_cASV_P_ordered_filtered)[1],2:dim(phenotype_cASV_P_ordered_filtered)[2]]),100,129)
 phenotype_cASV_Y_ordered_numeric_matrix <- matrix(as.numeric(phenotype_cASV_Y_ordered_filtered[2:dim(phenotype_cASV_Y_ordered_filtered)[1],2:dim(phenotype_cASV_Y_ordered_filtered)[2]]),100,129)
 phenotype_cASV_ave_ordered_numeric_matrix <- phenotype_cASV_P_ordered_numeric_matrix+phenotype_cASV_Y_ordered_numeric_matrix  / 2
 phenotype_cASV_ave_ordered_numeric_matrix  <- cbind(phenotype_cASV_Y_ordered[2:length(phenotype_cASV_Y_ordered[,1]),1],phenotype_cASV_ave_ordered_numeric_matrix)
 phenotype_cASV_ave_ordered_numeric_matrix  <- rbind(phenotype_cASV_Y_ordered[1,],phenotype_cASV_ave_ordered_numeric_matrix)
 
 phenotype_cASV_ave_ordered_numeric_matrix_filtered <- phenotype_cASV_ave_ordered_numeric_matrix
 
 # Remove bulk and parentals
 
 phenotype_cASV_P_ordered_filtered <- phenotype_cASV_P_ordered_filtered[!is.na(phenotype_cASV_P_ordered_filtered[,1]),]
 phenotype_cASV_Y_ordered_filtered <- phenotype_cASV_Y_ordered_filtered[!is.na(phenotype_cASV_Y_ordered_filtered[,1]),]
 phenotype_cASV_ave_ordered_numeric_matrix_filtered <- phenotype_cASV_ave_ordered_numeric_matrix_filtered[!is.na(phenotype_cASV_ave_ordered_numeric_matrix_filtered[,1]),]
 
 phenotype_cASV_P_file <- file.path(inp_folder, sprintf("phenotype_cASV_P_file.txt"))
 phenotype_cASV_Y_file <- file.path(inp_folder, sprintf("phenotype_cASV_Y_file.txt"))
 phenotype_cASV_ave_file <- file.path(inp_folder, sprintf("phenotype_cASV_ave_file.txt"))
 
 write.table(phenotype_cASV_P_ordered_filtered, file=phenotype_cASV_P_file, sep=",", quote=F, col.names=F, row.names=F)
 write.table(phenotype_cASV_Y_ordered_filtered, file=phenotype_cASV_Y_file, sep=",", quote=F, col.names=F, row.names=F)
 write.table(phenotype_cASV_ave_ordered_numeric_matrix_filtered, file=phenotype_cASV_ave_file, sep=",", quote=F, col.names=F, row.names=F)
 
 
 
 
 #####################################################################
 ################# Phenotypes: Section 3)  genus level ###############
 #####################################################################
 
 # now with the genus
 
 
 P_replicate_genus_matrix <- ps_all_genus_taxglom@otu_table[which(T_A_table$Replicate[match(rownames(ps_all_genus_taxglom@otu_table),T_A_table$Samples)]=="P"),]
 Y_replicate_genus_matrix <- ps_all_genus_taxglom@otu_table[which(T_A_table$Replicate[match(rownames(ps_all_genus_taxglom@otu_table),T_A_table$Samples)]=="Y"),]
 
 
 RIL_P_genus <- sapply(strsplit(rownames(P_replicate_genus_matrix), "_"), `[`, 2)
 RIL_Y_genus <- sapply(strsplit(rownames(Y_replicate_genus_matrix), "_"), `[`, 2)
 
 rownames(P_replicate_genus_matrix) <- RIL_P_genus
 rownames(Y_replicate_genus_matrix) <- RIL_Y_genus
 
 phenotype_genus_P <- cbind(rownames(P_replicate_genus_matrix),P_replicate_genus_matrix)
 phenotype_genus_P <- rbind(colnames(phenotype_genus_P),phenotype_genus_P)
 phenotype_genus_P[1,1]<- "id"
 
 phenotype_genus_Y <- cbind(rownames(Y_replicate_genus_matrix),Y_replicate_genus_matrix)
 phenotype_genus_Y <- rbind(colnames(phenotype_genus_Y),phenotype_genus_Y)
 phenotype_genus_Y[1,1] <- "id"
 
 
 colnames(phenotype_genus_P) <- NULL
 rownames(phenotype_genus_P) <- NULL
 
 colnames(phenotype_genus_Y) <- NULL
 rownames(phenotype_genus_Y) <- NULL
 
 phenotype_genus_P_ordered <- phenotype_genus_P[match(c("id",RIL_id_marker_table),phenotype_genus_P[,1]),]
 phenotype_genus_Y_ordered <- phenotype_genus_Y[match(c("id",RIL_id_marker_table),phenotype_genus_Y[,1]),]
 
 # Remove outlying and replicates that grew much longer
 phenotype_genus_P_ordered_filtered <- phenotype_genus_P_ordered
 phenotype_genus_Y_ordered_filtered <- phenotype_genus_Y_ordered
 phenotype_genus_P_ordered_filtered[which(phenotype_genus_P_ordered[,1] %in% P_removed),2:dim(phenotype_genus_P_ordered)[2]] <- NA
 phenotype_genus_Y_ordered_filtered[which(phenotype_genus_Y_ordered[,1] %in% Y_removed),2:dim(phenotype_genus_Y_ordered)[2]] <- NA
 
 
 # Remove bulk and parentals
 
 phenotype_genus_P_ordered_filtered <- phenotype_genus_P_ordered_filtered[!is.na(phenotype_genus_P_ordered_filtered[,1]),]
 phenotype_genus_Y_ordered_filtered <- phenotype_genus_Y_ordered_filtered[!is.na(phenotype_genus_Y_ordered_filtered[,1]),]
 
 phenotype_genus_P_file <- file.path(inp_folder, sprintf("phenotype_genus_P_file.txt"))
 phenotype_genus_Y_file <- file.path(inp_folder, sprintf("phenotype_genus_Y_file.txt"))
 
 write.table(phenotype_genus_P_ordered_filtered, file=phenotype_genus_P_file, sep=",", quote=F, col.names=F, row.names=F)
 write.table(phenotype_genus_Y_ordered_filtered, file=phenotype_genus_Y_file, sep=",", quote=F, col.names=F, row.names=F)
 
 
 ########### ########## ############
 ########### Testing GGG ###########
 ########## ########## #############
 # This makes 2 extra marker files #
 ########## ########## #############
 
 # repeat with GGG subsets
 GGG_subsets <- read.table(file="GGG_subsets_2_9_2014_plant2.txt", header=TRUE)
 GGG_subset1_pos <- which(marker_table[,1] %in% GGG_subsets[,1])
 GGG_subset2_pos <- which(marker_table[,1] %in% GGG_subsets[,2])
 marker_table_GGG1 <- marker_table[c(1,2,3,GGG_subset1_pos),]
 marker_table_GGG2 <- marker_table[c(1,2,3,GGG_subset2_pos),]
 
 RIL_id_marker_table_GGG1 <- marker_table_GGG1$V1[-c(1,2,3)]
 RIL_id_marker_table_GGG2 <- marker_table_GGG2$V1[-c(1,2,3)]
 marker_file <- file.path(inp_folder, sprintf("genotype_file.txt"))
 marker_file1 <- file.path(inp_folder, sprintf("genotype_file_GGG1.txt"))
 marker_file2 <- file.path(inp_folder, sprintf("genotype_file_GGG2.txt"))
 
 write.table(marker_table, file=marker_file, sep=",", quote=F, col.names=F, row.names=F)
 write.table(marker_table_GGG1, file=marker_file1, sep=",", quote=F, col.names=F, row.names=F)
 write.table(marker_table_GGG2, file=marker_file2, sep=",", quote=F, col.names=F, row.names=F)
 
 # Read files and transform into a readable format for R/qtl2
 genfile <- file.path(inp_folder, sprintf("genotype_file.txt"))
 genfile1 <- file.path(inp_folder, sprintf("genotype_file_GGG1.txt"))
 genfile2 <- file.path(inp_folder, sprintf("genotype_file_GGG2.txt"))
 
 
 ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
 ### Marker tables (plant genotypes) and phenotype tables are now constructed. Lets find microQTLs!
 ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
 
 # the first three rows are the name of the marker, the  chromosome and the location on the  chromosome.
 # phefile <- file.path(inp_folder, sprintf("%s_phenotype_taxlevel%s.txt", prefix, taxon_level))
 
 phefile <- file.path(inp_folder, sprintf("phenotype_PDW_file.txt"))
 phefile_Pm <- file.path(inp_folder, sprintf("phenotype_cASV_P_file.txt"))
 phefile_Ym <- file.path(inp_folder, sprintf("phenotype_cASV_Y_file.txt"))
 phefile_Avem <- file.path(inp_folder, sprintf("phenotype_cASV_ave_file.txt"))
 
 phefile_Pm_ASV <- file.path(inp_folder, sprintf("phenotype_ASV_P_file.txt"))
 phefile_Ym_ASV <- file.path(inp_folder, sprintf("phenotype_ASV_Y_file.txt"))
 phefile_Pm_genus <- file.path(inp_folder, sprintf("phenotype_genus_P_file.txt"))
 phefile_Ym_genus <- file.path(inp_folder, sprintf("phenotype_genus_Y_file.txt"))
 
 # plant dry weight
 in_cross <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile, crosstype = "riself", sep=",")
 # plant dry weight GGG1
 in_cross1 <- qtl::read.cross("csvs", genfile = genfile1, phefile = phefile, crosstype = "riself", sep=",")
 # plant dry weight GGG2
 in_cross2 <- qtl::read.cross("csvs", genfile = genfile2, phefile = phefile, crosstype = "riself", sep=",")
 # cASV
 in_crossPm <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Pm, crosstype = "riself", sep=",")
 in_crossYm <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Ym, crosstype = "riself", sep=",")
 in_crossAvem <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Avem, crosstype = "riself", sep=",")
 # ASV
  in_crossPm_ASV <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Pm_ASV, crosstype = "riself", sep=",")
 in_crossYm_ASV <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Ym_ASV, crosstype = "riself", sep=",")
 # Genus
 in_crossPm_genus <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Pm_ASV, crosstype = "riself", sep=",")
 in_crossYm_genus <- qtl::read.cross("csvs", genfile = genfile, phefile = phefile_Ym_ASV, crosstype = "riself", sep=",")
 
 
 mycross <- (qtl2::convert2cross2(in_cross))
 mycross1 <- (qtl2::convert2cross2(in_cross1))
 mycross2 <- (qtl2::convert2cross2(in_cross2))
 mycrossPm <- (qtl2::convert2cross2(in_crossPm))
 mycrossYm <- (qtl2::convert2cross2(in_crossYm))
 mycrossAvem <- (qtl2::convert2cross2(in_crossAvem))
 mycrossPm_ASV <- (qtl2::convert2cross2(in_crossPm_ASV))
 mycrossYm_ASV <- (qtl2::convert2cross2(in_crossYm_ASV))
 mycrossPm_genus <- (qtl2::convert2cross2(in_crossPm_genus))
 mycrossYm_genus <- (qtl2::convert2cross2(in_crossYm_genus))
 
 ##########################################
 # Some plots that can be nice to look at 
 ##########################################
 plot(in_cross, pheno.col = NULL)
 plot(in_cross1, pheno.col = NULL)
 plot(in_cross2, pheno.col = NULL)
 
# This shows which genotypes at each marker
plotGeno(in_cross, chr=1)
plotGeno(in_cross, chr=2) 
plotGeno(in_cross, chr=3) 
plotGeno(in_cross, chr=4) 

# how much missing information
plotInfo(in_cross, chr=1)

# phenotype distribution
plotPheno(in_cross3, pheno.col=3)
plot.map(in_cross)


################################################

# Creates markers at the given interval
#map <- insert_pseudomarkers(mycross$gmap, step=1)

 # Calculates the probability of a genotype at each marker
 pr <- calc_genoprob(mycross, mycross$gmap, error_prob=0.0001, cores=3)
 pr1 <- calc_genoprob(mycross1, mycross$gmap, error_prob=0.0001, cores=3)
 pr2 <- calc_genoprob(mycross2, mycross$gmap, error_prob=0.0001, cores=3)
 pr_Pm <- calc_genoprob(mycrossPm, mycrossPm$gmap, error_prob=0.0001, cores=3)
 pr_Ym <- calc_genoprob(mycrossYm, mycrossYm$gmap, error_prob=0.0001, cores=3)
 pr_Avem <- calc_genoprob(mycrossAvem, mycrossAvem$gmap, error_prob=0.0001, cores=3)
 pr_Pm_ASV <- calc_genoprob(mycrossPm_ASV, mycrossPm_ASV$gmap, error_prob=0.0001, cores=3)
 pr_Ym_ASV <- calc_genoprob(mycrossYm_ASV, mycrossYm_ASV$gmap, error_prob=0.0001, cores=3)
 pr_Pm_genus <- calc_genoprob(mycrossPm_genus, mycrossPm_genus$gmap, error_prob=0.0001, cores=3)
 pr_Ym_genus <- calc_genoprob(mycrossYm_genus, mycrossYm_genus$gmap, error_prob=0.0001, cores=3)
 
 
 
 # Performing a genome scan
 out <- scan1(pr, mycross$pheno[,2:4])
 out1 <- scan1(pr1, mycross1$pheno[,2:4])
 out2 <- scan1(pr2, mycross2$pheno[,2:4])
 outPm <- scan1(pr_Pm, mycrossPm$pheno[,2:dim(mycrossPm$pheno)[2]])
 outYm <- scan1(pr_Ym, mycrossYm$pheno[,2:dim(mycrossYm$pheno)[2]])
 outAvem <- scan1(pr_Avem, mycrossAvem$pheno[,2:dim(mycrossAvem$pheno)[2]])
 outPm_ASV <- scan1(pr_Pm_ASV, mycrossPm_ASV$pheno[,2:dim(mycrossPm_ASV$pheno)[2]])
 outYm_ASV <- scan1(pr_Ym_ASV, mycrossYm_ASV$pheno[,2:dim(mycrossYm_ASV$pheno)[2]])
 outPm_genus <- scan1(pr_Pm_genus, mycrossPm_genus$pheno[,2:dim(mycrossPm_genus$pheno)[2]])
 outYm_genus <- scan1(pr_Ym_genus, mycrossYm_genus$pheno[,2:dim(mycrossYm_genus$pheno)[2]])
 
 
 # Contains highest LOD score for each phenotype across the entire genome
 out_permutations_Pm_ASV <- scan1perm(pr_Pm_ASV, mycrossPm_ASV$pheno[,2:dim(mycrossPm_ASV$pheno)[2]],n_perm=1000)
 out_permutations_Ym_ASV <- scan1perm(pr_Ym_ASV, mycrossYm_ASV$pheno[,2:dim(mycrossYm_ASV$pheno)[2]],n_perm=1000)
 out_permutations_Pm_genus <- scan1perm(pr_Pm, mycrossPm_genus$pheno[,2:dim(mycrossPm_genus$pheno)[2]],n_perm=1000)
 out_permutations_Ym_genus <- scan1perm(pr_Ym, mycrossYm_genus$pheno[,2:dim(mycrossYm_genus$pheno)[2]],n_perm=1000)
 
 LOD_Pm_ASV <- quantile(out_permutations_Pm_ASV, probs=c(.95,.99,.999))
 LOD_Ym_ASV <- quantile(out_permutations_Ym_ASV, probs=c(.95,.99,.999))
 LOD_Pm_genus <- quantile(out_permutations_Pm_genus, probs=c(.95,.99,.999))
 LOD_Ym_genus <- quantile(out_permutations_Ym_genus, probs=c(.95,.99,.999))
 
 # Finding LOD peaks
 # find_peaks(out, mycross$gmap, threshold=2, drop=1.5)
 # find_peaks(out1, mycross$gmap, threshold=2, drop=1.5)
 # find_peaks(out2, mycross$gmap, threshold=2, drop=1.5)
 # cASV
 micro_QTL_Pm <- find_peaks(outPm, mycrossPm$gmap, threshold= 2, drop=1.5)
 micro_QTL_Ym <-  find_peaks(outYm, mycrossYm$gmap, threshold= 2, drop=1.5)
 micro_QTL_Avem <-  find_peaks(outAvem, mycrossAvem$gmap, threshold=2, drop=1.5)
 micro_QTL_Pm_ASV <- find_peaks(outPm_ASV, mycrossPm_ASV$gmap, threshold=  LOD_Pm_ASV[1] , drop=1.5)
 micro_QTL_Ym_ASV <-  find_peaks(outYm_ASV, mycrossYm_ASV$gmap, threshold= LOD_Ym_ASV[1], drop=1.5)
 micro_QTL_Pm_genus <- find_peaks(outPm_genus, mycrossPm_genus$gmap, threshold=2, drop=1.5)
 micro_QTL_Ym_genus <-  find_peaks(outYm_genus, mycrossYm_genus$gmap, threshold=2, drop=1.5)
 

 ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  #### 
 # The number of QTLs identified is 565 in Yellow, 522 in Pink. 
 length(micro_QTL_Ym_ASV$lodindex)
 length(micro_QTL_Pm_ASV$lodindex)
 
 # Of these, we identified 256 and 278 unique locations per replicate
 # Of these, 120 were replicated between the duplicate RIL populations
 # This is more managable, and could be prioritazation (1)
 chr_pos_ASV_QTLs_Ym <- paste(micro_QTL_Ym_ASV$chr,micro_QTL_Ym_ASV$pos)
 chr_pos_ASV_QTLs_Pm <- paste(micro_QTL_Pm_ASV$chr,micro_QTL_Pm_ASV$pos)
 replicated_chr_pos_ASV <- intersect(unique(chr_pos_ASV_QTLs_Ym), unique(chr_pos_ASV_QTLs_Pm))
 
 ######### 
 
 percent_relicated_ASV <- length(replicated_chr_pos)/length(unique(c(chr_pos_ASV_QTLs_Ym,chr_pos_ASV_QTLs_Pm)))
 table_ASV_QTLs_Ym_ASV <- table(paste(micro_QTL_Ym_ASV$chr,micro_QTL_Ym_ASV$pos))
 table_ASV_QTLs_Pm_ASV <- table(paste(micro_QTL_Pm_ASV$chr,micro_QTL_Pm_ASV$pos))
 
 # but we need to confidence intervals that go with these regions because now it is just a marker
 priority1_Ym <- micro_QTL_Ym_ASV[chr_pos_ASV_QTLs_Ym %in% replicated_chr_pos_ASV,]
 priority1_Pm <- micro_QTL_Ym_ASV[chr_pos_ASV_QTLs_Pm %in% replicated_chr_pos_ASV,]
 priority1 <- rbind(priority1_Ym,priority1_Pm)
 write.table(priority1, file="prioritized_QTLs_priority1", quote=FALSE, sep="\t")
 write.table(priority1_Ym, file="prioritized_QTLs_priority1_Y", quote=FALSE, sep="\t")
 write.table(priority1_Pm, file="prioritized_QTLs_priority1_P", quote=FALSE, sep="\t")
 

 # Interestingly, we have 120 prioritized QTL marker locations - with 543 rows, suggesting approximately 4.5 ASV per marker
 # We can prioritize this further by only taking 'high-impact' positions. This can both be because it is linked to many ASV, or to few but abundant ASV
 # First we find the 'many' ASV locations
 number_ASV_per_marker <- table(c(chr_pos_ASV_QTLs_Ym,chr_pos_ASV_QTLs_Pm))
 number_ASV_per_marker_replicated <- number_ASV_per_marker[names(number_ASV_per_marker) %in% replicated_chr_pos_ASV]
 plot(as.numeric(sort(number_ASV_per_marker_replicated,decreasing=TRUE)))
 abline(h=quantile(number_ASV_per_marker_replicated)[4],lty=2)
 high_impact_manyASV <- number_ASV_per_marker_replicated[number_ASV_per_marker_replicated>5]
 high_impact_manyASV_names <- names(high_impact_manyASV)
 
 high_impact_manyASV_qtl_Ym <- micro_QTL_Ym_ASV[which(chr_pos_ASV_QTLs_Ym %in% high_impact_manyASV_names),]
 high_impact_manyASV_qtl_Pm <- micro_QTL_Pm_ASV[which(chr_pos_ASV_QTLs_Pm %in% high_impact_manyASV_names),]
 
 many_priority1 <- rbind(high_impact_manyASV_qtl_Ym,high_impact_manyASV_qtl_Pm)
 write.table(many_priority1, file="prioritized_QTLs_many_priority1", quote=FALSE, sep="\t")
 
 

 # Next we find the 'abundant' ASV locations
 abundance_of_priority1_ASVs <- colSums(ps_all@otu_table[,priority1$lodindex])
 plot(sort(abundance_of_priority1_ASVs,decreasing=TRUE))
 abline(h=quantile(abundance_of_priority1_ASVs,.9),lty=2)
 high_abundance_priority1_ASVs <- abundance_of_priority1_ASVs[abundance_of_priority1_ASVs>quantile(abundance_of_priority1_ASVs,.9)]
 high_abundance_priority1 <- priority1[priority1$lodcolumn %in% names(high_abundance_priority1_ASVs),]
 
 write.table(high_abundance_priority1, file="prioritized_QTLs_high_abundance_priority1", quote=FALSE, sep="\t")
 
 # plotting this
 plot_peaks(priority1, mycrossPm_ASV$gmap)
 plot_peaks(high_abundance_priority1, mycrossPm_ASV$gmap)
 plot_peaks(many_priority1, mycrossPm_ASV$gmap)
 
 plot_peaks(many_priority1, mycrossPm_ASV$gmap)
 
 
 ##########
 # sort(micro_QTL_Ym_ASV$lod, decreasing = TRUE)[1:10]
 # micro_QTL_Pm_ASV[rev(order(micro_QTL_Pm_ASV$lod)),][1:10,]
 # micro_QTL_Ym_ASV[rev(order(micro_QTL_Ym_ASV$lod)),][1:10,]
 #########
 
 # To further subset these 278 unique locations, we decided to look at locations linked to a particular taxonomy
 Genus_with_priority_QTLs <- unique(ps_all@tax_table[priority1$lodindex,6])
 
 # length(unique(ps_all@tax_table[priority1$lodindex,6]))/length(unique(ps_all@tax_table[,6]))
 
 all_taxonomy_with_qtls <- ps_all@tax_table[unique(c(micro_QTL_Pm_ASV$lodindex,micro_QTL_Ym_ASV$lodindex))]
 
 setwd("/mnt/nfs/bioinfdata/ngs/ME/raaijmakers_group/results_RIL2019/")
 for (i in 1:length(Genus_with_priority_QTLs)) {
 taxon_with_qtls_pos <- which(all_taxonomy_with_qtls[,6] == Genus_with_priority_QTLs[[i]])
 # taxon_with_qtls_pos <- which(all_taxonomy_with_qtls[,6] == "Sphingobium")
 
 taxon_asv_with_qtls <- rownames(all_taxonomy_with_qtls)[taxon_with_qtls_pos]
 
 all_taxon_qtls_peaks <- rbind(micro_QTL_Pm_ASV[which(micro_QTL_Pm_ASV$lodcolumn %in% taxon_asv_with_qtls),],
                               micro_QTL_Ym_ASV[which(micro_QTL_Ym_ASV$lodcolumn %in% taxon_asv_with_qtls),]) 
 
 write.table(all_taxon_qtls_peaks, file=paste(c("prioritized_QTLs_",Genus_with_priority_QTLs[[i]]), sep="",collapse=""), quote=FALSE, sep="\t")
 }
 
 
 plot_peaks(all_taxon_qtls_peaks, mycrossPm_ASV$gmap)
 
 #############
 par(mfrow=c(1,4))
 priority1_Ym[41,]
 gP_rhizobium_ASV301 <- maxmarg(pr_Ym_ASV, mycrossYm_ASV$gmap, chr=2, pos=41.27, return_char=TRUE)
 gY <- maxmarg(pr_Ym_ASV, mycrossYm_ASV$gmap, chr=5, pos=1.062, return_char=TRUE)
 gY <- maxmarg(pr_Ym_ASV, mycrossYm_ASV$gmap, chr=5, pos=1.062, return_char=TRUE)
 
 gP_sphinogbium <- maxmarg(pr_Pm_ASV, mycrossPm_ASV$gmap, chr=6, pos=102.148, return_char=TRUE)
 
 
  # Concerning the genotype information: A is money, B is pimp.
 plot_pxg(gY, mycrossYm_ASV$pheno[,"ASV15"], ylab="CSS", main = "ASV15 Y")
 plot_pxg(gY, mycrossYm_ASV$pheno[,"ASV33"], ylab="CSS", main = "ASV15 Y")
 plot_pxg(gY, mycrossYm_ASV$pheno[,"ASV34"], ylab="CSS", main = "ASV15 Y")
 plot_pxg(gY, mycrossYm_ASV$pheno[,"ASV85"], ylab="CSS", main = "ASV15 Y")

 par(mfrow=c(1,2))
 strepo_hawaiiensis<-c("ASV12","ASV15","ASV23","ASV24","ASV25","ASV26","ASV27","ASV30","ASV33","ASV34","ASV36","ASV47","ASV48","ASV50","ASV52","ASV55","ASV56","ASV96","ASV97","ASV99","ASV105","ASV85")
 strepo_eurythermus<-c("ASV57","ASV58","ASV72","ASV73","ASV76","ASV78","ASV81","ASV112","ASV114","ASV116")
 replicated_strepto<-c("ASV34","ASV105","ASV123","ASV124","ASV172","ASV188")
    
 strepto_sum_strepo_replicated_Ym <- rowSums(mycrossYm_ASV$pheno[,replicated_strepto],na.rm=TRUE)
 strepto_sum_strepo_replicated_Pm <- rowSums(mycrossPm_ASV$pheno[,replicated_strepto],na.rm=TRUE)
 
 strepto_sum_strepo_hawaiiensis_Ym <- rowSums(mycrossYm_ASV$pheno[,strepo_hawaiiensis],na.rm=TRUE)
 strepto_sum_strepo_hawaiiensis_Pm <- rowSums(mycrossPm_ASV$pheno[,strepo_hawaiiensis],na.rm=TRUE)
 
 strepto_sum_strepo_eurythermus_Ym <- rowSums(mycrossYm_ASV$pheno[,strepo_eurythermus],na.rm=TRUE)
 strepto_sum_strepo_eurythermus_Pm <- rowSums(mycrossPm_ASV$pheno[,strepo_eurythermus],na.rm=TRUE)
 
 plot_pxg(gY[accessions_NA], strepto_sum_strepo_replicated_Ym[accessions_NA], ylab="CSS", main = i)
 plot_pxg(gP[accessions_NA], strepto_sum_strepo_replicated_Pm[accessions_NA], ylab="CSS", main = i)

 
  
 par(mfrow=c(3,6), mar = c(2,2,2,2))
for (i in strepo_hawaiiensis[10:18]) {
  plot_pxg(gY, mycrossYm_ASV$pheno[,i], ylab="CSS", main = i, col = "yellow")
  plot_pxg(gP, mycrossPm_ASV$pheno[,i], ylab="CSS", main = i, col = "pink")
  
  }
 
 par(mfrow=c(4,6), mar = c(2,2,2,2))
 for (i in strepo_eurythermus[1:18]) {
   plot_pxg(gY, mycrossYm_ASV$pheno[,i], ylab="CSS", main = i, col = "yellow")
   plot_pxg(gP, mycrossPm_ASV$pheno[,i], ylab="CSS", main = i, col = "pink")
   
 }
 
 par(mfrow=c(3,6), mar = c(2,2,2,2))
 for (i in strepo_eurythermus[1:18]) {
   plot_pxg(gP, mycrossPm_ASV$pheno[,i], ylab="CSS", main = i)
 }
 
 par(mfrow=c(1,2), mar = c(2,2,2,2))
accessions_NA <- !is.na(rowSums(mycrossYm_ASV$pheno[,strepo_hawaiiensis]))
 plot_pxg(gY[accessions_NA], strepto_sum_strepo_hawaiiensis_Ym[accessions_NA], ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062")
 plot_pxg(gP[accessions_NA], strepto_sum_strepo_hawaiiensis_Pm[accessions_NA], ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062")
 
 t.test(strepto_sum_strepo_hawaiiensis_Ym[accessions_NA][which(gY[accessions_NA]=="BB")],strepto_sum_strepo_hawaiiensis_Ym[accessions_NA][which(gY[accessions_NA]=="AA")])
 t.test(strepto_sum_strepo_hawaiiensis_Pm[accessions_NA][which(gP[accessions_NA]=="BB")],strepto_sum_strepo_hawaiiensis_Pm[accessions_NA][which(gP[accessions_NA]=="AA")])
 
 
 par(mfrow=c(1,2), mar = c(2,2,2,2))
 accessions_NA <- !is.na(rowSums(mycrossYm_ASV$pheno[,strepo_eurythermus]))
 plot_pxg(gY[accessions_NA], strepto_sum_strepo_eurythermus_Ym[accessions_NA], ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062")
 plot_pxg(gP[accessions_NA], strepto_sum_strepo_eurythermus_Pm[accessions_NA], ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062")
 
 t.test(strepto_sum_strepo_eurythermus_Ym[accessions_NA][which(gY[accessions_NA]=="BB")],strepto_sum_strepo_eurythermus_Ym[accessions_NA][which(gY[accessions_NA]=="AA")])
 t.test(strepto_sum_strepo_eurythermus_Pm[accessions_NA][which(gP[accessions_NA]=="BB")],strepto_sum_strepo_eurythermus_Pm[accessions_NA][which(gP[accessions_NA]=="AA")])
 
 
 
 
 strepto_sum_Ym <- mycrossYm_ASV$pheno[,"ASV15"]+mycrossYm_ASV$pheno[,"ASV33"]+mycrossYm_ASV$pheno[,"ASV34"]+mycrossYm_ASV$pheno[,"ASV85"]+mycrossYm_ASV$pheno[,"ASV116"]+mycrossYm_ASV$pheno[,"ASV123"]+mycrossYm_ASV$pheno[,"ASV124"]+mycrossYm_ASV$pheno[,"ASV171"]+mycrossYm_ASV$pheno[,"ASV172"]+mycrossYm_ASV$pheno[,"ASV188"]
 strepto_sum_Pm <- mycrossPm_ASV$pheno[,"ASV15"]+mycrossPm_ASV$pheno[,"ASV33"]+mycrossPm_ASV$pheno[,"ASV34"]+mycrossPm_ASV$pheno[,"ASV85"]+mycrossPm_ASV$pheno[,"ASV116"]+mycrossPm_ASV$pheno[,"ASV123"]+mycrossPm_ASV$pheno[,"ASV124"]+mycrossPm_ASV$pheno[,"ASV171"]+mycrossPm_ASV$pheno[,"ASV172"]+mycrossPm_ASV$pheno[,"ASV188"]
 
 sphingo_sum_Ym <- mycrossYm_ASV$pheno[,"ASV39"]+mycrossYm_ASV$pheno[,"ASV40"]+mycrossYm_ASV$pheno[,"ASV43"]+mycrossYm_ASV$pheno[,"ASV60"]+mycrossYm_ASV$pheno[,"ASV63"]+ mycrossYm_ASV$pheno[,"ASV88"]+mycrossYm_ASV$pheno[,"ASV92"]+mycrossYm_ASV$pheno[,"ASV93"]+mycrossYm_ASV$pheno[,"ASV108"]+mycrossYm_ASV$pheno[,"ASV110"]+mycrossYm_ASV$pheno[,"ASV120"]+mycrossYm_ASV$pheno[,"ASV122"]+mycrossYm_ASV$pheno[,"ASV136"]
 sphingo_sum_Pm <- mycrossPm_ASV$pheno[,"ASV39"]+mycrossPm_ASV$pheno[,"ASV40"]+mycrossPm_ASV$pheno[,"ASV43"]+mycrossPm_ASV$pheno[,"ASV60"]+mycrossPm_ASV$pheno[,"ASV63"]+ mycrossPm_ASV$pheno[,"ASV88"]+mycrossPm_ASV$pheno[,"ASV92"]+mycrossPm_ASV$pheno[,"ASV93"]+mycrossPm_ASV$pheno[,"ASV108"]+mycrossPm_ASV$pheno[,"ASV110"]+mycrossPm_ASV$pheno[,"ASV120"]+mycrossPm_ASV$pheno[,"ASV122"]+mycrossPm_ASV$pheno[,"ASV136"]

par(mfrow=c(1,2))
plot_pxg(gY, strepto_sum_Ym, ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062", ylim=c(0,50))
# plot_pxg(gY, strepto_sum_Pm, ylab="CSS", main = "Streptomyces \nchr=5, pos=1.062", ylim=c(0,50))
# plot_pxg(gP_sphinogbium, sphingo_sum_Ym, ylab="CSS", main = "Sphingo \nchr=6, pos=102.148", ylim=c(0,50))
plot_pxg(gP_sphinogbium, sphingo_sum_Pm, ylab="CSS", main = "Sphingo \nchr=6, pos=102.148", ylim=c(0,50))

 
 pos_ASV301 <- which(colnames(ps_all@otu_table)=="ASV301")
 plot_pxg(gP_rhizobium_ASV301, mycrossYm_ASV$pheno[,"ASV301"], ylab="CSS", main = "ASV301 Rhizobium")
 ASV301_RIL <- as.numeric(ps_all@otu_table[which(ps_all@sam_data$Treatment=="RIL"),pos_ASV301])
 ASV301_M <- as.numeric(ps_all@otu_table[which(ps_all@sam_data$Treatment=="M"),pos_ASV301])
 ASV301_P <- as.numeric(ps_all@otu_table[which(ps_all@sam_data$Treatment=="P"),pos_ASV301])
 ASV301_BULK <- as.numeric(ps_all@otu_table[which(ps_all@sam_data$Treatment=="BULK"),pos_ASV301])
 
plot(ASV301_RIL,ASV301_M,ASV301_P,ASV301_BULK, names = c("RIL","Pimpi","Money","Bulk"), ylab = "Abundance", main = "ASV301 Rhizobium")
 
 
 
 
 sall_strepto_qtls_peaks[order(all_strepto_qtls_peaks$chr),]
 

  
 sort(table_ASV_QTLs_Ym_ASV, decreasing=TRUE)
 hot_spot1 <- which(chr_pos_ASV_QTLs_Ym == "5 55.277")
 micro_QTL_Ym_ASV[hot_spot1,]
 micro_QTL_Ym_ASV$lodindex[hot_spot1]
 
 
sort(table(paste(micro_QTL_Ym_ASV$chr[which(micro_QTL_Ym_ASV$lod>3)],micro_QTL_Ym_ASV$pos[which(micro_QTL_Ym_ASV$lod>3)])),decreasing=TRUE)
length(unique(paste(micro_QTL_Ym_ASV$chr[which(micro_QTL_Ym_ASV$lod>3)],micro_QTL_Ym_ASV$pos[which(micro_QTL_Ym_ASV$lod>3)])))

 
 length(micro_QTL_Ym$lodindex)
 length(micro_QTL_Pm$lodindex)
 length(micro_QTL_Avem$lodindex)
 
 unique(c(micro_QTL_Ym$lodindex,micro_QTL_Pm$lodindex,micro_QTL_Avem$lodindex))
 unique(c(micro_QTL_Ym$pos,micro_QTL_Pm$pos,micro_QTL_Avem$pos))
 intersect(micro_QTL_Ym$lodindex,micro_QTL_Pm$lodindex)
 intersect(intersect(micro_QTL_Ym$lodindex,micro_QTL_Pm$lodindex),micro_QTL_Avem$lodindex)
 
 genus_and_cluster[which(genus_and_cluster[,2]%in%micro_QTL$lodcolumn),]
 c_ASV_qtl_pos <- which(colnames(c_ASV_matrix) %in% genus_and_cluster[which(genus_and_cluster[,2]%in%micro_QTL_Pm$lodcolumn),])
 c_ASV_matrix[c_ASV_qtl_pos,]

  phenotype_cASV_P_ordered_filtered[,1]
  phenotype_cASV_P_ordered_filtered[,113]
  
  QTL_location <- names(which(mycrossPm$gmap$`6`=="102.148"))
  # QTL_location <- names(which(mycrossPm$gmap$`8`=="80.365"))
  
 allele_of_cASV_QTL <- mycrossPm$geno$`6`[,which(colnames(mycrossPm$geno$`6`) == QTL_location)]
 # allele_of_cASV_QTL <- mycrossPm$geno$`8`[,which(colnames(mycrossPm$geno$`8`) == QTL_location)]
 
 intersect_of_gen_phen_P <- match(phenotype_cASV_P_ordered_filtered[-1,1],names(allele_of_cASV_QTL))
 intersect_of_gen_phen_Y <- match(phenotype_cASV_Y_ordered_filtered[-1,1],names(allele_of_cASV_QTL))
 
 cASV_allele_Table <- cbind(allele_of_cASV_QTL,NA)
 cASV_allele_Table <- cbind(cASV_allele_Table,NA)
 cASV_allele_Table[intersect_of_gen_phen_P,2] <- as.numeric(phenotype_cASV_P_ordered_filtered[-1,113])
 cASV_allele_Table[intersect_of_gen_phen_Y,3] <- as.numeric(phenotype_cASV_Y_ordered_filtered[-1,113])
 # cASV_allele_Table[intersect_of_gen_phen,2] <- as.numeric(phenotype_cASV_P_ordered_filtered[-1,4])
 
 colnames(cASV_allele_Table)<- c("allele","rel_ab_P","rel_ab_Y")
 # cASV_allele_Table <- as.data.frame(cASV_allele_Table)
 # cASV_allele_Table <- cASV_allele_Table[!is.na(cASV_allele_Table$rel_ab_P),]
#  cASV_allele_Table$rel_ab <- as.numeric(as.character(cASV_allele_Table$rel_ab))
 boxplot(cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==1)],
         cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==2)],
         cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==1)],
         cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==2)],
         names = c("P_1","P_2","Y_1","Y_2"))
         
 
 boxplot(c(cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==1)],cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==1)]),
         c(cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==2)],cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==2)]),
         names = c("1","2"))
 
 t.test(c(cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==1)],cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==1)]),
         c(cASV_allele_Table$rel_ab_P[which(cASV_allele_Table$allele==2)],cASV_allele_Table$rel_ab_Y[which(cASV_allele_Table$allele==2)]))
 
 allele_of_cASV_QTL <- mycrossYm$geno$`6`[,which(colnames(mycrossYm$geno$`6`) == QTL_location)]
 intersect_of_gen_phen <- match(phenotype_cASV_Y_ordered_filtered[-1,1],names(allele_of_cASV_QTL))
 cASV_allele_Table <- cbind(allele_of_cASV_QTL,NA)
 cASV_allele_Table[intersect_of_gen_phen,2] <- phenotype_cASV_Y_ordered_filtered[-1,114]
 colnames(cASV_allele_Table)<- c("allele","rel_ab")
 cASV_allele_Table <- as.data.frame(cASV_allele_Table)
 cASV_allele_Table <- cASV_allele_Table[!is.na(cASV_allele_Table$rel_ab),]
 cASV_allele_Table$rel_ab <- as.numeric(cASV_allele_Table$rel_ab)
 boxplot(cASV_allele_Table$rel_ab[which(cASV_allele_Table$allele==1)],cASV_allele_Table$rel_ab[which(cASV_allele_Table$allele==2)])
 
 
 plot(out, mycross$gmap, lodcolumn=1, col="gray21", ylim=c(0, 5.0),lty=2, lwd=2, main ="Full RIL population")
 plot(out, mycross$gmap, lodcolumn=2, col="gray21", add=TRUE,lty=1, lwd=2)
 plot(out, mycross$gmap, lodcolumn=3, col="violetred", add=TRUE, lwd=3)
  
  boxplot(alle_of_cASV_QTL[intersect_of_gen_phen], phenotype_cASV_P_ordered_filtered[intersect_of_gen_phen,113])
  
  

 # cASV84
 # cASV90
 # cASV98
 # cASV112
 # cASV123
 # Permutation of genome scan using Haley-Knott regression # This takes the longest
 out_perm <- scan1perm(pr3, mycross$pheno[,2:4], n_perm=1000, cores=3)

 # Contains highest LOD score for each phenotype across the entire genome
 perm_thres <- summary(out_perm, alpha=c(0.15, 0.05))

 par(mfrow=c(3,1),mar=c(1,2,3,1))
 
 plot(out, mycross$gmap, lodcolumn=1, col="gray21", ylim=c(0, 5.0),lty=2, lwd=2, main ="Full RIL population")
 plot(out, mycross$gmap, lodcolumn=2, col="gray21", add=TRUE,lty=1, lwd=2)
 plot(out, mycross$gmap, lodcolumn=3, col="violetred", add=TRUE, lwd=3)
 abline(h=3, lty=c(5,1))
 legend("topright", lwd=c(3,3,3), lty=c(2,1,1), col=c("gray21", "gray21","violetred"), colnames(out), bg="gray90")
 
 plot(out1, mycross$gmap, lodcolumn=1, col="gray21", ylim=c(0, 5.0),lty=2, lwd=2, main ="GGG subset 1")
 plot(out1, mycross$gmap, lodcolumn=2, col="gray21", add=TRUE,lty=1, lwd=2)
 plot(out1, mycross$gmap, lodcolumn=3, col="violetred", add=TRUE, lwd=3)
 abline(h=3, lty=c(5,1))
 legend("topright", lwd=c(3,3,3), lty=c(2,1,1), col=c("gray21", "gray21","violetred"), colnames(out), bg="gray90")
 
 plot(out2, mycross$gmap, lodcolumn=1, col="gray21", ylim=c(0, 5.0),lty=2, lwd=2, main ="GGG subset 2")
 plot(out2, mycross$gmap, lodcolumn=2, col="gray21", add=TRUE,lty=1, lwd=2)
 plot(out2, mycross$gmap, lodcolumn=3, col="violetred", add=TRUE, lwd=3)
 abline(h=3, lty=c(5,1))
 legend("topright", lwd=c(3,3,3), lty=c(2,1,1), col=c("gray21", "gray21","violetred"), colnames(out), bg="gray90")
 

abline(h=perm_thres, lty=c(5,1))

############## Try to improve LOD scores

PY_differences <- as.numeric(phenotype_PDW[-1,2])-as.numeric(phenotype_PDW[-1,3])
PY_differences_nNA <- !is.na(as.numeric(phenotype_PDW[-1,2])-as.numeric(phenotype_PDW[-1,3]))
PY_differences[PY_differences_nNA]
ave_diff_PDW <- ave(PY_differences[PY_differences_nNA])[1]
sd_diff_PDW <- sd(PY_differences[PY_differences_nNA])[1]
cutoffs<-c((ave_diff_PDW+2*sd_diff_PDW),(ave_diff_PDW-2*sd_diff_PDW))

which(PY_differences > cutoffs[1])
which(PY_differences < cutoffs[2])
