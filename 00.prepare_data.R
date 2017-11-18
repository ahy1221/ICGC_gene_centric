source("./functions.R")

library(dplyr)
#======= Load the data
#--- RNA alteration
gt <- readr::read_tsv("./data/master_table.tsv.gz")

#--- Clinical data
clinical <- readr::read_tsv("./data/pcawg_donor_clinical_August2016_v7.tsv")

#--- Expression 
expr <- readr::read_tsv("../pcawg_gene_fusions/data/expressions/tophat_star_fpkm_uq.v2_aliquot_gl.tsv") %>% as.data.frame()
rownames(expr) <- expr$feature
expr <- as.matrix(expr[, 2:ncol(expr)])
#--- convert To TPM
expr.tpm <- apply(expr, 2, fpkm_to_tpm)

#--- Meta
meta <- readr::read_tsv("../pcawg_gene_fusions/data/meta/pcawg.rnaseq.extended.metadata.aliquot_id.V4.tsv") %>% as.data.frame()
aliquot2donor <- meta[, c("icgc_donor_id", "aliquot_id")]
rownames(aliquot2donor) <- aliquot2donor$aliquot_id
donor2cancer <- unique(meta[, c('icgc_donor_id', 'histology_abbreviation')])

#===== Collapse multiple aliquot id
aliquots.use <- colnames(expr.tpm)[colnames(expr.tpm) %in% aliquot2donor$aliquot_id]
#--- To donor
aliquot2donor <- aliquot2donor[aliquots.use, ]
aliquot2donor <- aliquot2donor[order(aliquot2donor$icgc_donor_id, aliquot2donor$aliquot_id),]
expr.tpm <- expr.tpm[,aliquot2donor$aliquot_id]
expr.tpm.donor <- rowsum(t(expr.tpm), group = aliquot2donor$icgc_donor_id)
expr.tpm.donor <- t(sweep(expr.tpm.donor, 1, table(aliquot2donor$icgc_donor_id), "/"))


#---- subset donors
expr.tpm.use <- expr.tpm.donor[,intersect(unique(gt$ICGC_DONOR_ID), colnames(expr.tpm.donor))]
#---- subset genes
rownames(expr.tpm.use) <- sapply(strsplit(rownames(expr.tpm.use), "\\."), `[[`, 1)
#---- To symbol
gene.table <- unique(data.frame(id = gt$Ensembl_Gene_ID, symbol = gt$hgnc_symbol, stringsAsFactors = F))
expr.tpm.use <- expr.tpm.use[gene.table$id,]
rownames(expr.tpm.use) <- gene.table$symbol
donor2cancer.use <- donor2cancer[donor2cancer$icgc_donor_id %in% colnames(expr.tpm.use),]
#---- Duplicated donors
#donor2cancer.use[duplicated(donor2cancer.use$icgc_donor_id),]
donor2cancer.use <- donor2cancer.use[!(donor2cancer.use$icgc_donor_id ==  "DO45281" & donor2cancer.use$histology_abbreviation != "Liver-HCC"),]
#---- Clean
rm(list = c("expr", "expr.tpm", "expr.tpm.donor", "aliquots.use", "aliquot2donor", "donor2cancer"))

save(clinical, donor2cancer.use, expr.tpm.use, gene.table, gt, meta, file = "./data/gene.centric.rds")

