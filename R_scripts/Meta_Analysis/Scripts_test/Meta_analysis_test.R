library(mixOmics)

#######################################################################################################################
## STEPNIAK DATA
#######################################################################################################################
data_stepniak <- read.table("../../Glioma/Stepniak/RNAseq/nextflow/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE)
#--> Phenotypic matrix 
pdata_stepniak <- read.csv("pdata_DA_complete.csv", sep =";", header = TRUE)
#pdata <- pdata[-6,]
row.names(pdata_stepniak) <- pdata_stepniak$Sample

#--> Count matrix
mrna_raw_counts_stepniak <- read.table("../../Glioma/Stepniak/RNAseq/nextflow/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE)
row.names(mrna_raw_counts_stepniak) <- mrna_raw_counts_stepniak$gene_id
mrna_raw_counts_stepniak <- mrna_raw_counts_stepniak[,-c(1:2)]
mrna_raw_counts_stepniak <- mrna_raw_counts_stepniak[,c(1,4,3,6,5,7)]

#--> Order and match
mrna_raw_counts_stepniak <- mrna_raw_counts_stepniak[,match(pdata_stepniak$Sample, colnames(mrna_raw_counts_stepniak), nomatch = 0)]
pdata_stepniak <- pdata_stepniak[pdata_stepniak$Sample %in% colnames(mrna_raw_counts_stepniak),]

#--> Filter lowly expressed genes
library(HTSFilter)

mrna_filter_stepniak <- HTSFilter(round(mrna_raw_counts_stepniak),
                         pdata_stepniak$Group,
                         s.min=1, s.max=200)

mrna_filter_stepniak <- mrna_filter_stepniak$filteredData

#--> Feature annotation
mrna_gene_stepniak <- read.table("../../Glioma/Stepniak/RNAseq/nextflow/star_salmon/salmon.merged.gene_counts.tsv", header = TRUE)[,c(1:2)]

## Normalisation
library(RUVSeq)
library(reshape)
library(EDASeq)

#--> edgeR TMM normalization
library(edgeR)

tmm_set_stepniak <-  mrna_filter_stepniak

pdata_tmm_stepniak <- pdata_stepniak
design <- cbind(model.matrix(~ 0 + pdata_tmm_stepniak$Group))
colnames(design) <- gsub(".*\\$|)", "", colnames(design))

tmm_set_stepniak <- DGEList(counts = tmm_set_stepniak,
                   genes = mrna_gene_stepniak[match(
                     row.names(mrna_filter_stepniak), 
                     mrna_gene_stepniak$gene_id),],
                   group = pdata_tmm_stepniak$Group)

tmm_set_stepniak <- calcNormFactors(tmm_set_stepniak)

#--> Estimate dispersion and fit the model
tmm_set_stepniak <- estimateGLMRobustDisp(tmm_set_stepniak, design)

##### Count per million
tmm_cpm_stepniak <- cpm(tmm_set_stepniak, normalized.lib.sizes = TRUE, log = TRUE)

# --> Transpose the df to have the sample in lines
t_mrna_norm_stepniak <- as.data.frame(t(tmm_cpm_stepniak))

#######################################################################################################################
## JIGNA DATA
#######################################################################################################################
data_jigna <- mrna_raw_counts <- read.table("salmon.merged.gene_counts_Jigna.tsv", header = TRUE)
#--> Phenotypic matrix 
pdata_jigna <- read.csv("pdata_Jigna.csv", sep =";", header = TRUE)
#pdata <- pdata[-6,]
row.names(pdata_jigna) <- pdata_jigna$Sample

#--> Count matrix
mrna_raw_counts_jigna <- read.table("salmon.merged.gene_counts_Jigna.tsv", header = TRUE)
row.names(mrna_raw_counts_jigna) <- mrna_raw_counts_jigna$gene_id
mrna_raw_counts_jigna <- mrna_raw_counts_jigna[,-c(1:2)]
mrna_raw_counts_jigna <- mrna_raw_counts_jigna[,c(1,3,5,7,9,11,13,15)]

#--> Order and match
mrna_raw_counts_jigna <- mrna_raw_counts_jigna[,match(pdata_jigna$Sample, colnames(mrna_raw_counts_jigna), nomatch = 0)]
pdata_jigna <- pdata_jigna[pdata_jigna$Sample %in% colnames(mrna_raw_counts_jigna),]

#--> Filter lowly expressed genes
library(HTSFilter)

mrna_filter_jigna <- HTSFilter(round(mrna_raw_counts_jigna),
                         pdata_jigna$Group,
                         s.min=1, s.max=200)

mrna_filter_jigna <- mrna_filter_jigna$filteredData

#--> Feature annotation
mrna_gene_jigna <- read.table("salmon.merged.gene_counts_Jigna.tsv", header = TRUE)[,c(1:2)]

###### Normalization ######
library(RUVSeq)
library(reshape)
library(EDASeq)

#--> edgeR TMM normalization
library(edgeR)

tmm_set_jigna <-  mrna_filter_jigna

pdata_tmm_jigna <- pdata_jigna
design <- cbind(model.matrix(~ 0 + pdata_tmm_jigna$Group))
colnames(design) <- gsub(".*\\$|)", "", colnames(design))

tmm_set_jigna <- DGEList(counts = tmm_set_jigna,
                   genes = mrna_gene_jigna[match(
                     row.names(mrna_filter_jigna), 
                     mrna_gene_jigna$gene_id),],
                   group = pdata_tmm_jigna$Group)

tmm_set_jigna <- calcNormFactors(tmm_set_jigna)

#--> Estimate dispersion and fit the model
tmm_set_jigna <- estimateGLMRobustDisp(tmm_set_jigna, design)

##### Count per million
tmm_cpm_jigna <- cpm(tmm_set_jigna, normalized.lib.sizes = TRUE, log = TRUE)

# --> Transpose the df to have the sample in lines
t_mrna_norm_jigna <- as.data.frame(t(tmm_cpm_jigna))


#######################################################################################################################
## COMBIBNE THE 2 DATAFRAMES
common_cols <- intersect(colnames(t_mrna_norm_jigna), colnames(t_mrna_norm_stepniak))

t_mrna_norm_jigna_common <- t_mrna_norm_jigna[, common_cols]
t_mrna_norm_stepniak_common <- t_mrna_norm_stepniak[, common_cols]
  
X_combined_count_mrna <- rbind(t_mrna_norm_jigna_common,t_mrna_norm_stepniak_common)
Y_pdata <- read.csv("pdata.csv", sep =";", header = TRUE)

class <- as.factor(Y_pdata$Group)
study <- factor(Y_pdata$Study)
#######################################################################################################################
## Comparison with previous analysis
genes_Jigna_downreg <- read.csv("../Differential_Analysis/Differential_analysis_Jigna/downregulated_genes_Jigna_samples.csv", header = TRUE)
genes_Jigna_upreg <- read.csv("../Differential_Analysis/Differential_analysis_Jigna/upregulated_genes_Jigna_samples.csv", header = TRUE)
#genes_Jigna_downreg_metabo <- read.csv("../Differential_Analysis/Differential_analysis_Jigna/downregulated_metabo_genes_Jigna_samples.csv", header = FALSE)
#genes_Jigna_upreg_metabbo <- read.csv("../Differential_Analysis/Differential_analysis_Jigna/upregulated_metabo_genes_Jigna_samples.csv", header = FALSE)

genes_Stepniak_downreg <- read.csv("../Differential_Analysis/Differential_Analysis_Stepniak/RNA/downregulated_genes_DA.csv", header = TRUE)
genes_Stepniak_upreg <- read.csv("../Differential_Analysis/Differential_Analysis_Stepniak/RNA/upregulated_genes_DA.csv", header = TRUE)
#genes_Stepniak_upreg_metabo <- read.csv("../Differential_Analysis/Differential_Analysis_Stepniak/RNA/upregulated_metabo_genes_DA.csv", header = TRUE)
#genes_Stepniak_downreg_metabo <- read.csv("../Differential_Analysis/Differential_Analysis_Stepniak/RNA/downregulated_metabo_genes_DA.csv", header = TRUE)

gene_names <- mrna_gene_stepniak[match(
  common_cols, 
  mrna_gene_stepniak$gene_id),]

metabo_genes_all <- read.table("../Differential_Analysis/Lists_genes/all_metabo.tsv", header = TRUE, sep = "\t")
metabo_genes <- metabo_genes_all$Gene_name

meta_analysis_genes <- gene_names$gene_name
down_genes_Jigna <- as.character(genes_Jigna_downreg$x)
up_genes_Jigna <- as.character(genes_Jigna_upreg$x)
genes_jigna <- c(down_genes_Jigna,up_genes_Jigna)

Jigna_genes = list(Downregulated = down_genes_Jigna, meta_analysis_genes = meta_analysis_genes, Metabo = metabo_genes, Upregulated = up_genes_Jigna)

# Create the Venn diagram
venn_plot_dysreg_jigna <- ggvenn(
  Jigna_genes,
)
venn_plot_dysreg_jigna

#######
meta_analysis_genes <- gene_names$gene_name
down_genes_Stepniak <- as.character(genes_Stepniak_downreg$x)
up_genes_Stepniak <- as.character(genes_Stepniak_upreg$x)
genes_stepniak <- c(down_genes_Stepniak,up_genes_Stepniak)

Stepniak_genes = list(Downregulated = down_genes_Stepniak, meta_analysis_genes = meta_analysis_genes, Metabo = metabo_genes, Upregulated = up_genes_Stepniak)

# Create the Venn diagram
venn_plot_dysreg_stepniak <- ggvenn(
  Stepniak_genes,
)
venn_plot_dysreg_stepniak




#######################################################################################################################

## MINT sPLS-DA Model
result <- mint.splsda(X = X_combined_count_mrna, Y = class, study = study, ncomp = 2)
plotIndiv(result, comp = c(1, 2), legend = TRUE)

## TUNING THE NUMBER OF COMPONENTS
splsda_perf <- perf(result)
plot(splsda_perf)
splsda_perf$global.error

optim_ncomp <- 2

## TUNING THE NUMBER OF FEATURES
splsda_tune <- tune(X = X_combined_count_mrna, 
                    Y = class, 
                    study = study, 
                    ncomp = optim_ncomp,
                    test.keepX = seq(1, 100, 1),
                    nrepeat = 5,
                    method = "mint.splsda",
                    measure = "BER",
                    dist = "centroids.dist",
                    )

plot(splsda_tune, sd = FALSE)

optim_keepX <- splsda_tune$choice.keepX # extract optimal values

optim_keepX

## FINAL MODEL
final_splsda_model <- mint.splsda(X = X_combined_count_mrna, 
                                  Y = class, 
                                  study = study, 
                                  ncomp = optim_ncomp, 
                                  keepX = optim_keepX)

cim(final_splsda_model, comp = 2, margins=c(10,5), 
    row.sideColors = color.mixo(as.numeric(class)), 
    row.names = FALSE, title = "MINT sPLS-DA, component 2")


selected_gene_comp2 <- selectVar(final_splsda_model, comp = 2)$name

gene_names <- mrna_gene_stepniak[match(
  selected_gene_comp2, 
  mrna_gene_stepniak$gene_id),]


