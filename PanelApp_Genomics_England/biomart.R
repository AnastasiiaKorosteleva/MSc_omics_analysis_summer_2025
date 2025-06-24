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


# > head(bed_df)
# chrom    start      end    name
# 338  chr1  1001137  1014540   ISG15
# 633  chr1  1211325  1214153 TNFRSF4
# 634  chr1  7915870  7943165 TNFRSF9
# 458  chr1  9629888  9729114  PIK3CD
# 381  chr1 11022008 11047239   MASP2
# 402  chr1 11106534 11262556    MTOR
