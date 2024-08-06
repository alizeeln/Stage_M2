sbatch --wrap="\
TOBIAS ATACorrect --bam ./Glioma/Stepniak/ATACseq/nextflow/bwa/merged_library/DA01_ATACseq2_S26_REP1.mLb.clN.sorted.bam --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa --peaks ./Glioma/Stepniak/ATACseq/nextflow/bwa/merged_library/macs2/narrow_peak/DA01_ATACseq2_S26_REP1.mLb.clN_summits.bed --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed --outdir ./Glioma/Stepniak/ATACseq/ATACorrect_test --cores 8"

###########################################################################################################################

sbatch --wrap="\
TOBIAS ATACorrect \
                --bam ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/DA01_ATACseq2_S26_REP1.mLb.clN.sorted.bam \
                --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA01_ATACseq2_S26_REP1.mLb.clN_summits.bed \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --outdir ./Glioma/Stepniak/ATACseq/ATACorrect_DA01 \
                --cores 8"

sbatch --wrap="\
TOBIAS ATACorrect \
                --bam ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_replicate/DA02_ATACseq2.mRp.clN.sorted.bam \
                --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_replicate/macs2/narrow_peak/DA02_ATACseq2.mRp.clN_summits.bed \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --outdir ./Glioma/Stepniak/ATACseq/ATACorrect_DA02 \
                --cores 8"

sbatch --wrap="\
TOBIAS ATACorrect \
                --bam ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/DA05_ATACseq_S26_REP1.mLb.clN.sorted.bam \
                --genome ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --peaks ./Glioma/Stepniak/ATACseq/nextflow2/bwa/merged_library/macs2/narrow_peak/DA05_ATACseq_S26_REP1.mLb.clN_summits.bed \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --outdir ./Glioma/Stepniak/ATACseq/ATACorrect_DA05 \
                --cores 8"

