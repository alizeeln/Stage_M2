#### Script to perform cim because it didn't work in a .rmd file
## use it after charging all the object from Integrative_Analysis_RNA_chip.rmd file

pdf("heatmap_plot_3.pdf", width = 25, height = 25)  # Adjust dimensions as needed
 # Call the plot function
print(cimDiablo(final.diablo.model,
                row.names = TRUE,
                col.names = TRUE,
                color.Y = c("#D60000", "#4294D8"),
                margins = c(25, 25),
                comp = 1))

dev.off()


sim <- cimDiablo(final.diablo.model,
          row.names = TRUE,
          col.names = TRUE,
          color.Y = c("#D60000", "#4294D8"),
          margins = c(25, 25),
          comp = 1)

cluster_1_50 <- labels(sim$ddc[[1]])
cluster_2_50 <- labels(sim$ddc[[2]])

metabo_genes_all <- read.table("../Differential_Analysis/Lists_genes/all_metabo.tsv", header = TRUE, sep = "\t")
metabo_genes <- metabo_genes_all$Gene_name

## Function to extract the name of genes before __
extract_before_double_underscore <- function(x) {
  sapply(strsplit(x, "__"), `[`, 1)
}

##### Cluster 2
cluster_1_50_clean <- extract_before_double_underscore(cluster_1_50)
indices_1 <- which(cluster_1_50_clean %in% metabo_genes)
metabo_clust_1_50 <- cluster_1_50[indices_1]

##### Cluster 2
cluster_2_50_clean <- extract_before_double_underscore(cluster_2_50)
indices_2 <- which(cluster_2_50_clean %in% metabo_genes)
metabo_clust_2_50 <- cluster_2_50[indices_2]

######################################################################

