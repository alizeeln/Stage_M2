sbatch --wrap="\
TOBIAS FootprintScores \
                --signal ./Glioma/Stepniak/ATACseq/ATACorrect_DA01/DA01_ATACseq2_S26_REP1.mLb.clN.sorted_corrected.bw \
                --regions ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA01_ATACseq2_S26_REP1.mLb.clN_summits.bed \
                --output DA01_ATAC_footprints.bw \
                --cores 8"

./Glioma/Stepniak/ATACseq/ATACScoreBW_DA01/

sbatch --wrap="\
TOBIAS FootprintScores \
                --signal ./Glioma/Stepniak/ATACseq/ATACorrect_DA02/DA02_ATACseq2.mRp.clN.sorted_corrected.bw \
                --regions ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_replicate/macs2/narrow_peak/DA02_ATACseq2.mRp.clN_summits.bed \
                --output DA02_ATAC_footprints.bw \
                --cores 8"

./Glioma/Stepniak/ATACseq/ATACScoreBW_DA02/

sbatch --wrap="\
TOBIAS FootprintScores \
                --signal ./Glioma/Stepniak/ATACseq/ATACorrect_DA05/DA05_ATACseq_S26_REP1.mLb.clN.sorted_corrected.bw \
                --regions ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA05_ATACseq_S26_REP1.mLb.clN_summits.bed \
                --output DA05_ATAC_footprints.bw \
                --cores 8"

./Glioma/Stepniak/ATACseq/ATACScoreBW_DA05/