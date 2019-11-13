#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(goseq)
  library(GO.db)
  library(dplyr)
})


args <- commandArgs(trailingOnly = T)

p_cutoff <- 1e-2
dataset_fname <- args[[1]]

message(sprintf("Reading %s ...", dataset_fname))
dataset <- read.csv(dataset_fname, colClasses = "character", sep = "\t")
goseq_geneids <- read.csv("goseq_geneids.csv", colClasses = c("factor", "numeric", "factor"), sep = "\t")

all_gi <- unique(goseq_geneids[, c(1, 2)])

significant_terms <- data.frame(
  lod = NULL,
  term = NULL,
  in_region = NULL,
  in_genome = NULL,
  pval_over = NULL,
  pval_under = NULL,
  go_name = NULL,
  go_ns = NULL
)
cat(paste(c("snp", "go_term", "in_region", "in_genome", "pval_over", "pval_under", "go_name", "go_namespace"), collapse = '\t'), "\n")
for(snp in unique(dataset$lod)) {
  lod_gi <- subset(dataset, lod == snp)$geneid
  v <- rep(0, nrow(all_gi))
  names(v) <- all_gi$geneid
  v[names(v) %in% lod_gi] <- 1

  pwf <- nullp(v, genome = "sl2.40", id = "geneid", bias.data = all_gi$gene_length, plot.fit = F)
  # pwf <- nullp(v, genome = "sl2.40", id = "geneid", bias.data = rnorm(length(v), 500, 1), plot.fit = F)

  go_overrep <- goseq(pwf = pwf, genome = "sl2.40", id = "geneid", gene2cat = goseq_geneids[, c(1,3)], repcnt = 2000, use_genes_without_cat = F)

  go_filtered <- subset(go_overrep, over_represented_pvalue <= p_cutoff) # | under_represented_pvalue <= p_cutoff)
  if(nrow(go_filtered) < 1)
    next()

  for(go in 1:nrow(go_filtered)) {
    cat(paste(c(snp, go_filtered[go, c(1, 4, 5, 2, 3, 6, 7)]), collapse = '\t'), "\n")
    significant_terms <- rbind(
      significant_terms,
      data.frame(
        lod = snp,
        term = go_filtered[go, 1],
        in_region = go_filtered[go, 4],
        in_genome = go_filtered[go, 5],
        pval_over = go_filtered[go, 2],
        pval_under = go_filtered[go, 3],
        go_name = go_filtered[go, 6],
        go_ns = go_filtered[go, 7]
      )
    )
  }
}


# genfile_reformatted <- genfile_reformatted %>%
#   arrange(chromosome, centimorgans) %>%
#   group_by(chromosome) %>%
#   mutate(
#     start = lag(locus, default = 0, order_by = chromosome),
#     end = lead(locus, default = last(locus), order_by = chromosome),
#     width = end - start,
#     start_cm = lag(centimorgans, default = 0, order_by = chromosome),
#     end_cm = lead(centimorgans, default = last(centimorgans), order_by = chromosome),
#     width_cm = end_cm - start_cm
#   )
# write.table(genfile_reformatted, file = "markers_as_table.csv", sep = "\t", row.names = F, col.names = T, quote = F)
