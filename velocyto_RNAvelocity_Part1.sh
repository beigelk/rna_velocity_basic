
#!/bin/bash

# # RNA Velocity analysis using velocyto, Seurat, and scVelo (Part 1: velocyto)
# ## Use velocyto for getting spliced and unspliced counts for RNA velocity
# ### Kat Beigel (beigelk@chop.edu)


# Generated Loom files using the sorted BAM files
# (possorted_genome_bam.bam) output from Cell Ranger (`cellranger count`)
# as suggested by 10x Genomics Trajectory Analysis using 10x Genomics Single Cell Gene Expression Data
# (https://www.10xgenomics.com/resources/analysis-guides/trajectory-analysis-using-10x-Genomics-single-cell-gene-expression-data)

# Main velocyto docs: https://velocyto.org/velocyto.py/
# Using `run` command in velocyto: https://velocyto.org/velocyto.py/tutorial/cli.html#run-run-on-any-technique-advanced-use

# The best option for GTFFILE (genome annotation file) and the --mask file (optional) is to get the genome version that was used
# for getting counts. In this case, the CellRanger pipeline used the refdata-gex-mm10-2020-A genome,
# so we will use that here as well.

# Run velocyto
velocyto run --outputfolder /dir/where/to/write/data/ \
--samtools-threads 16 --samtools-memory 2000 \
--bcfile scRNA_data/WT/outs/filtered_feature_bc_matrix/barcodes.tsv \
--mask /home/beigelk/references/mm38/GRCm38_rmsk.gtf \
scRNA_data/WT/outs/possorted_genome_bam.bam \
/home/beigelk/references/mm38/cellranger_ref/refdata-gex-mm10-2020-A/genes/genes.gtf

# The resulting file (e.g., 'velocytyo_possorted_genome_bam_WT.loom', can be used for Part 3: scVelo RNA velocity)