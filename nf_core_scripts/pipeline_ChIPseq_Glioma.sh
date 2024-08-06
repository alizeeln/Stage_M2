sbatch --wrap="\
nextflow run nf-core/chipseq -profile ifb_core -resume -r 2.0.0\
                --input ./Glioma/Stepniak/ChiPseq/samplesheet_ChIPseq_2.csv \
                --outdir ./Glioma/Stepniak/ChiPseq/nextflow2/ \
                --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
                --bwa_index ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/bwa \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --save_unaligned \
                --narrow_peak \
                --macs_gsize 2652783500"

sbatch --wrap="\
nextflow run nf-core/chipseq -profile ifb_core -resume -r 2.0.0\
                --input ./Glioma/Stepniak/ChiPseq/samplesheet_ChIPseq_2_H3K27me3.csv \
                --outdir ./Glioma/Stepniak/ChiPseq/nextflowH3K27me3/ \
                --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
                --bwa_index ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/bwa \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --save_unaligned \
                --macs_gsize 2652783500"