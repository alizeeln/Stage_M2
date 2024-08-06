sbatch --wrap="\
TOBIAS BINDetect --motifs ./Glioma/Stepniak/ATACseq/motifs/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
                    --signals ./Glioma/Stepniak/ATACseq/DA01_ATAC_footprints.bw ./Glioma/Stepniak/ATACseq/DA02_ATAC_footprints.bw \
                    --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                    --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed \
                    --peak_header ./Glioma/Stepniak/ATACseq/peaks_header.txt \
                    --outdir ./Glioma/Stepniak/ATACseq/BINDetect \
                    --cond_names IDHmut IDHwt \
                    --cores 8"

sbatch --wrap="\
TOBIAS BINDetect --motifs ./Glioma/Stepniak/ATACseq/motifs/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
                    --signals ./Glioma/Stepniak/ATACseq/DA01_ATAC_footprints.bw ./Glioma/Stepniak/ATACseq/DA02_ATAC_footprints.bw \
                    --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.chr.fa \
                    --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.bed \
                    --outdir ./Glioma/Stepniak/ATACseq/BINDetect2 \
                    --naming 'name_id' \
                    --cores 8"


########
sbatch --wrap="\
TOBIAS BINDetect --motifs ./Glioma/Stepniak/ATACseq/motifs/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme \
                    --signals ./Glioma/Stepniak/ATACseq/DA01_ATAC_footprints.bw ./Glioma/Stepniak/ATACseq/DA02_ATAC_footprints.bw \
                    --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                    --peaks ./Glioma/Stepniak/ATACseq/consensus_peaks.mLb.clN_without_MT.bed \
                    --outdir ./Glioma/Stepniak/ATACseq/BINDetect3 \
                    --naming 'name' \
                    --cores 8"

