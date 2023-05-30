import pandas as pd

df = pd.read_csv("assets/make/hg19_refseq.csv")
df["strand"] = "."
df["gene_symbol"] = df.index.astype(str) + ":" + df["gene_symbol"]

df = df[["gene_symbol", "chr", "strand", "tss"]]
df.to_csv("assets/hg19_refseq.tsv", index=False, header=False, sep="\t")
