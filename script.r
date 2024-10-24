# Load required libraries
install.packages("clusterProfiler")  # Uncomment if clusterProfiler is not installed
install.packages("org.Hs.eg.db")     # Uncomment for human gene annotations
library(clusterProfiler)
library(org.Hs.eg.db)                 # Adjust if using another organism

# Step 1: Prepare your gene list
# Replace this with your list of differentially expressed genes (DEGs)
# For example, let's assume you have a vector of gene symbols
de_genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")

# Step 2: Convert gene symbols to Entrez IDs
entrez_ids <- bitr(de_genes, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Step 3: Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                               organism = 'hsa',  # Use 'hsa' for Homo sapiens
                               pvalueCutoff = 0.05)

# Step 4: View the enrichment results
head(kegg_enrichment)

# Step 5: Visualize the enrichment results
# Dot plot
dotplot(kegg_enrichment, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

# Bar plot
barplot(kegg_enrichment, showCategory = 10) + ggtitle("KEGG Pathway Enrichment")

# Step 6: Perform Gene Ontology (GO) enrichment analysis
go_enrichment <- enrichGO(gene = entrez_ids$ENTREZID, 
                           OrgDb = org.Hs.eg.db, 
                           keyType = "ENTREZID", 
                           ont = "BP",  # Biological Process
                           pvalueCutoff = 0.05)

# Step 7: View the GO enrichment results
head(go_enrichment)

# Step 8: Visualize the GO enrichment results
# Dot plot for GO
dotplot(go_enrichment, showCategory = 10) + ggtitle("GO Biological Process Enrichment")

# Bar plot for GO
barplot(go_enrichment, showCategory = 10) + ggtitle("GO Biological Process Enrichment")

# Step 9: Save results to a file
write.csv(as.data.frame(kegg_enrichment), file = "kegg_enrichment_results.csv")
write.csv(as.data.frame(go_enrichment), file = "go_enrichment_results.csv")
