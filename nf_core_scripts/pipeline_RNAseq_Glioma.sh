sbatch -A egattaca --wrap="\
nextflow run nf-core/rnaseq -r 3.11.0 -profile ifb_core -resume -c ./Glioma/Stepniak/RNAseq/custom_module.config \
       --input ./Glioma/Stepniak/RNAseq/samplesheet_RNA_glioma.csv \
       --outdir ./Glioma/Stepniak/RNAseq/nextflow/ \
       --fasta ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
       --gtf ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/Homo_sapiens.GRCh38.108.gtf \
       --star_index ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/STAR/ \
       --salmon_index ./Glioma/save/RESSOURCES/REFERENCES/GRCh38/salmon/"



############ TO test on 1 sample --> successfull ###############
sbatch -A egattaca --wrap="\
nextflow run nf-core/rnaseq -r 3.11.0 -profile ifb_core -resume -c ./Glioma/Stepniak/RNAseq/custom_module.config \
       --input ./Glioma/Stepniak/RNAseq/samplesheet_test.csv \
       --outdir ./Glioma/Stepniak/RNAseq/nextflow/ \
       --fasta ./Glioma/REFERENCES/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
       --gtf ./Glioma/REFERENCES/Homo_sapiens.GRCh38.111.gtf \
       --star_index ./Glioma/GRCh38/STAR/ \
       --salmon_index ./Glioma/GRCh38/salmon/"