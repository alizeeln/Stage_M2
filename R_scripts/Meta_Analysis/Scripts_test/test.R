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


X_combined_count_mrna <- t_mrna_norm_jigna

Y_pdata <- read.csv("pdata_test.csv", sep =";", header = TRUE)

class <- as.factor(Y_pdata$Group)
study <- factor(Y_pdata$Study)
#######################################################################################################################
## MINT sPLS-DA Model
result <- mint.splsda(X = X_combined_count_mrna, Y = class, study = study, ncomp = 3)
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
