# rna_velocity_basic
Basic RNA velocity analysis using velocyto (to get spliced/unspliced counts), Seurat (to get UMAP coords and other metadata), and scVelo (RNA velocity analysis and figures).

To run velocyto, you will need the .bam/.sam files for your data. See https://velocyto.org/velocyto.py/tutorial/index.html for instructions on running velocyto. See https://velocyto.org/velocyto.py/tutorial/cli.html#running-velocyto for more information on which `velocyto run` command you will need.

To run the Seurat script, you will need a fullt processed Seurat object with metadata and UMAP coordinates. Please note that this script was developed under Seurat v4 with a v4 object.

To run the scVelo script, you will need the .loom file from the velocyto output and the .h5ad file converted from the .h5Seurat file.

Note: part1 and part2 scripts both need to be run before part3, but part1 and part2 can be run in either order. 
