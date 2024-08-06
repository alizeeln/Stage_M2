sbatch --wrap="\
nextflow run nf-core/atacseq -profile ifb_core -resume \
                --input ./Glioma/Stepniak/ATACseq/samplesheet_ATACseq_2_Stepniak.csv \
                --outdir ./Glioma/Stepniak/ATACseq/nextflow2/ \
                --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --save_reference \
                --narrow_peak \
                --macs_gsize 2652783500"

sbatch --wrap="\
nextflow run nf-core/atacseq -profile ifb_core -resume \
                --input ./Glioma/Stepniak/ATACseq/samplesheet_ATACseq_2_Stepniak.csv \
                --outdir ./Glioma/Stepniak/ATACseq/nextflow2/ \
                --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
                --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
                --blacklist ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/hg38-blacklist.v2.bed \
                --narrow_peak \
                --macs_gsize 2652783500"