library(biomaRt)

COVID.19.research <- read.delim("~/COVID-19 research.tsv")
COVID.19.research$Entity.Name
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")


gene_list <- unique(COVID.19.research$Entity.Name)


coords <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = gene_list,
  mart = mart)

# Keep only standard chromosomes (1-22, X, Y)
coords <- coords[coords$chromosome_name %in% c(1:22, "X", "Y"), ]

# Create BED format: chrom, start-1, end, name
bed_df <- data.frame(
  chrom = paste0("chr", coords$chromosome_name),
  start = coords$start_position - 1,
  end = coords$end_position,
  name = coords$hgnc_symbol
)

# Sort BED
bed_df <- bed_df[order(bed_df$chrom, bed_df$start), ]

head(bed_df)
