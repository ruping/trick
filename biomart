library(biomaRt)
# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id"),
                 filters = "ensembl_transcript_id", values = "ENST00000296026",
                 mart = mart)
results
#   ensembl_gene_id ensembl_transcript_id ensembl_peptide_id
# 1 ENSG00000163734       ENST00000296026    ENSP00000296026