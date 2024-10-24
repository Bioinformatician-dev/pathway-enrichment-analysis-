# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # For human gene annotations

# Example gene list from differential expression analysis
gene_list <- c("BRCA1", "TP53", "EGFR", "VEGFA", "MYC")

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Perform pathway enrichment analysis using KEGG
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                            organism = 'hsa', 
                            pvalueCutoff = 0.05)

# View the results
head(kegg_results)
