# taxon level must be numeric
# list at taxon level must be at the same taxonomic level (e.g. taxon_level = 7 for genus, thus list_at_taxon_level must be Streptomyces)
# install.packages("seqinr")
library("seqinr")
grab_sig_ASV_fasta <- function(taxon_level,list_at_taxon_level) {

  for (i in 1:length(list_at_taxon_level)) {
  ASV_positions <- as.numeric(which(taxa.bacteria.asv[,taxon_level] == as.character(list_at_taxon_level[i])))
  ASV_names <- as.character(taxa.bacteria.asv[ASV_positions,1])
  taxon_priority1_Pm_Ym <- rbind(micro_QTL_Ym_ASV[which(micro_QTL_Ym_ASV[,2] %in% ASV_names),],micro_QTL_Pm_ASV[which(micro_QTL_Pm_ASV[,2] %in% ASV_names),])
  sig_genus_asv_index <- unique(taxon_priority1_Pm_Ym$lodindex)
  write.fasta(sequences=as.list(rownames(taxa.bacteria.asv)[sig_genus_asv_index]),names =ASV_names, file.out=paste(c("~/RIL_results/dada2_outputs/dir_of_sig_ASV_sequences/",list_at_taxon_level[[i]],"-sig-Y_or_P.fna"), sep="",collapse=""))
  }
}


grab_sig_ASV_fasta(7,Genus_with_priority_QTLs)
