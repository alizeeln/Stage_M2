basic_diablo_model <- block.splsda(X = data, Y = Y, ncomp = 2, design = design) 

#perf.diablo = perf(basic_diablo_model, validation = 'loo',nrepeat = 20) 

var_rna_comp1 <- selectVar(basic_diablo_model, comp = 1, block = 1)
var_rna_comp2 <- selectVar(basic_diablo_model, comp = 2, block = 1)

var_h3K4me3_comp1 <- selectVar(basic_diablo_model, comp = 1, block = 2)
var_h3K4me3_comp2 <- selectVar(basic_diablo_model, comp = 2, block = 2)

var_h3k27me3_comp1 <- selectVar(basic_diablo_model, comp = 1, block = 3)
var_h3k27me3_comp2 <- selectVar(basic_diablo_model, comp = 2, block = 3)

var_h3k27ac_comp1 <- selectVar(basic_diablo_model, comp = 1, block = 4)
var_h3k27ac_comp2 <- selectVar(basic_diablo_model, comp = 2, block = 4)

#top_var_rna_comp1 <- var_rna_comp1$RNAseq$name[order(var_rna_comp1$RNAseq$value, decreasing = TRUE)[1:100]]
#top_var_rna_comp2 <- var_rna_comp2$RNAseq$name[order(var_rna_comp1$RNAseq$value, decreasing = TRUE)[1:100]]

###### Top features RNA
top_var_rna_comp1 <- var_rna_comp1$RNAseq$value %>%
  arrange(desc(var_rna_comp1$RNAseq$value)) 

top_var_rna_comp1_100 <- rownames(top_var_rna_comp1)[1:100]

top_var_rna_comp2 <- var_rna_comp2$RNAseq$value %>%
  arrange(desc(var_rna_comp2$RNAseq$value)) 

top_var_rna_comp2_100 <- rownames(top_var_rna_comp2)[1:100]

###### Top features H3K4me3
top_var_h3k4me3_comp1 <- var_h3K4me3_comp1$H3K4me3$value %>%
  arrange(desc(var_h3K4me3_comp1$H3K4me3$value)) 

top_var_h3k4me3_comp1_100 <- rownames(top_var_h3k4me3_comp1)[1:100]

top_var_h3k4me3_comp2 <- var_h3K4me3_comp2$H3K4me3$value %>%
  arrange(desc(var_h3K4me3_comp2$H3K4me3$value)) 

top_var_h3k4me3_comp2_100 <- rownames(top_var_h3k4me3_comp2)[1:100]

###### Top features H3K27me3
top_var_h3k27me3_comp1 <- var_h3k27me3_comp1$H3K27me3$value %>%
  arrange(desc(var_h3k27me3_comp1$H3K27me3$value)) 

top_var_h3k27me3_comp1_100 <- rownames(top_var_h3k27me3_comp1)[1:100]

top_var_h3k27me3_comp2 <- var_h3k27me3_comp2$H3K27me3$value %>%
  arrange(desc(var_h3k27me3_comp2$H3K27me3$value)) 

top_var_h3k27me3_comp2_100 <- rownames(top_var_h3k27me3_comp2)[1:100]

###### Top features H3K27ac
top_var_h3k27ac_comp1 <- var_h3k27ac_comp1$H3K27ac$value %>%
  arrange(desc(var_h3k27ac_comp1$H3K27ac$value)) 

top_var_h3k27ac_comp1_100 <- rownames(top_var_h3k27ac_comp1)[1:100]

top_var_h3k27ac_comp2 <- var_h3k27ac_comp2$H3K27ac$value %>%
  arrange(desc(var_h3k27ac_comp2$H3K27ac$value)) 

top_var_h3k27ac_comp2_100 <- rownames(top_var_h3k27ac_comp2)[1:100]

########### Test to see if there was a huge drop to chose a better treshold but no)
#absolute <- abs(sorted_genes_df_desc)
#absolute_df_desc <- absolute %>%
#  arrange(desc(absolute$value.var))
#barplot(height=absolute_df_desc$value.var[1:300], ylim = c(0.018,0.0205))
###########