# Load the required libraries
library(GenomicFeatures)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.refGene)
library(org.Hs.eg.db)

# Load the RefSeq-based TxDb object
txdb <- TxDb.Hsapiens.UCSC.hg19.refGene

# Extract the TSS
tss <- promoters(txdb, upstream=0, downstream=1)

# Get the genes from the TxDb object, including genes on both strands
genes_obj <- genes(txdb, single.strand.genes.only=FALSE)

# Define a function to map the gene IDs
map_gene_ids <- function(gene_ids, column) {
  mapIds(org.Hs.eg.db, 
         keys = gene_ids, 
         column = column, 
         keytype = "ENTREZID", 
         multiVals = "first")
}

# Get the gene IDs from the genes object
gene_ids <- mcols(genes_obj)$gene_id

# Map the gene IDs to gene symbols and refseq ids
gene_symbols <- map_gene_ids(gene_ids, "SYMBOL")
gene_refseqs <- map_gene_ids(gene_ids, "REFSEQ")

# Create a data frame that maps gene ids to gene symbols and refseq ids
gene_df <- data.frame(gene_id = gene_ids, 
                      gene_symbol = gene_symbols,
                      gene_refseq = gene_refseqs,
                      stringsAsFactors = FALSE)

# Clean the data frame by removing NA gene symbols
gene_df <- gene_df[!is.na(gene_df$gene_symbol), ]

# Merge `tss` and `gene_df` on `tx_name` and `gene_refseq`
tss_with_genes <- merge(tss, gene_df, by.x="tx_name", by.y="gene_refseq")

# Remove entries not part of the main, canonical assembly
tss_with_genes <- tss_with_genes[!grepl("_", tss_with_genes$seqnames), ]

# Subset to NM_* identifiers
tss_with_genes_coding <- tss_with_genes[grepl("^NM_", tss_with_genes$tx_name), ]

# Select and reorder the columns
df <- tss_with_genes_coding[, c("seqnames", "start", "tx_name", "gene_symbol")]

# Rename the columns
names(df) <- c("chr", "tss", "refseq_id", "gene_symbol")

# Display the head of the newly formatted data frame
head(df)

# Write DataFrame to CSV file
write.csv(df, file = "assets/make/hg19_refseq.csv", row.names = FALSE)

