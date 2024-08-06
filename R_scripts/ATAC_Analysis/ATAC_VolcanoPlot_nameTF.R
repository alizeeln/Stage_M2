library(ggplot2)
library(ggrepel)

######## The raw ATACseq data were processed using nf-core v2.14.1 ATACseq v2.1.2 on the IFB 
# cluster with the following command line :
#sbatch --wrap="\
#nextflow run nf-core/atacseq -profile ifb_core -resume \
#                --input ./Glioma/Stepniak/ATACseq/samplesheet_ATACseq_2_Stepniak.csv \
#                --outdir ./Glioma/Stepniak/ATACseq/nextflow2/ \
#                --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#                --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
#                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
#                --narrow_peak \
#                --macs_gsize 2652783500"

###### TOBIAS was used as follows : 
### ATACorrect is used to correct the bias created by the Tn5 insertion, this command line
### have to be run independently for each sample

#sbatch --wrap="\
#TOBIAS ATACorrect \
#                --bam ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/DA01_ATACseq2_S26_REP1.mLb.clN.sorted.bam \
#                --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#                --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA01_ATACseq2_S26_REP1.mLb.clN_summits.bed \
#                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
#                --outdir ./Glioma/Stepniak/ATACseq/ATACorrect_DA01 \
#               --cores 8"

### Then, this tool allows to calculate a binding score of protein across the genome for each sample:

#sbatch --wrap="\
#TOBIAS FootprintScores \
#                --signal ./Glioma/Stepniak/ATACseq/ATACorrect_DA01/DA01_ATACseq2_S26_REP1.mLb.clN.sorted_corrected.bw \
#                --regions ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA01_ATACseq2_S26_REP1.mLb.clN_summits.bed \
#                --output ./Glioma/Stepniak/ATACseq/ATACScoreBW_DA01/DA01_ATAC_footprints.bw \
#                --cores 8"

### The final tool is BINDetect. It allows to associate at each footprint a possible Transcription
### Factor that binds at its place based on motifs. It is based on the output of this tool that the
### volcano plot can be generated
#sbatch --wrap="\
#TOBIAS BINDetect --motifs ./Glioma/Stepniak/ATACseq/motifs/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
#                    --signals ./Glioma/Stepniak/ATACseq/DA01_ATAC_footprints.bw ./Glioma/Stepniak/ATACseq/DA02_ATAC_footprints.bw \
#                    --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
#                    --peaks ./Glioma/Stepniak/ATACseq/consensus_peaks.mLb.clN_without_MT.bed \
#                    --outdir ./Glioma/Stepniak/ATACseq/BINDetect3 \
#                    --naming 'name' \
#                    --cores 8"

########## This script generates a volcano plot of the differentially binded Transcription Factors #####
########## It takes as an input the two files that are the output of TOBIAS BINDetect tool 

# Define your data paths
bindetect_results_path <- "../../../Glioma/Stepniak/ATACseq/Input_data_volcano_plot/bindetect_results.txt"
bindetect_results_copie_path <- "../../../Glioma/Stepniak/ATACseq/Input_data_volcano_plot/bindetect_results_copie.txt"

# Define constants
changelimitleft <- 0.075
changelimitright <- 0.2
pvallimit <- -log10(0.075)

# Read data
data <- read.table(bindetect_results_copie_path, header = TRUE, sep = " ")
name <- data[, 2]

topgenes <- read.table(bindetect_results_path, header = TRUE, sep = "")
topgenes <- cbind(name, topgenes)

# Calculate -log10(pvalue)
topgenes$log10pvalue <- -log10(topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_pvalue)

# Define the differential expression status
topgenes$diffexpressed <- "NO"
topgenes$diffexpressed[topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change > changelimitleft & topgenes$log10pvalue > pvallimit] <- "UP"
topgenes$diffexpressed[topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change < -changelimitright & topgenes$log10pvalue > pvallimit] <- "DOWN"

# Determine the thresholds for plotting
yvalues <- topgenes$log10pvalue
xvalues <- topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change

y_min <- quantile(yvalues[yvalues < -log10(1e-300)], 0.95)
x_min <- quantile(xvalues, 0.05)
x_max <- quantile(xvalues, 0.95)

# Determine which TFs to show
topgenes$show <- FALSE
topgenes$show[topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change < x_min | 
                topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change > x_max | 
                topgenes$log10pvalue > y_min] <- TRUE

# Add labels for significant genes
topgenes$labels <- NA
topgenes$labels[topgenes$show] <- topgenes$name[topgenes$show]

# Define colors
topgenes$color <- "black"
topgenes$color[topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change < x_min] <- "red"
topgenes$color[topgenes$DA01_ATAC_footprints_DA02_ATAC_footprints_change > x_max | 
                 topgenes$log10pvalue > y_min] <- "blue"

# Convert differential expression status to factor
topgenes$diffexpressed <- as.factor(topgenes$diffexpressed)

# Plot the results
p <- ggplot(data = topgenes, aes(x = DA01_ATAC_footprints_DA02_ATAC_footprints_change, y = log10pvalue, label = labels)) + 
  geom_point(aes(colour = color)) + 
  theme_classic() + 
  geom_text_repel(aes(colour = color), max.overlaps = Inf, show.legend = FALSE, size = 5) + 
  scale_colour_manual(values = c("red" = "red", "black" = "black", "blue" = "blue"), labels = c("Higher scores in DA01_ATAC_footprints", "No changes", "Higher scores in DA02_ATAC_footprints")) +
  geom_vline(xintercept = c(-changelimitleft, changelimitright), linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = pvallimit, linetype = "dashed", alpha = 0.4) + 
  theme(legend.position = c(0.8, 0.2), legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 20), axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15), axis.title = element_text(size = 18)) + 
  labs(x = "Differential binding score", y = expression("-log10(pvalue)")) + 
  expand_limits(x = c(-0.2, 0.2)) + 
  guides(colour = guide_legend(title = NULL))

print(p)
