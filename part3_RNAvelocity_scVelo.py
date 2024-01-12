#!/usr/bin/env python
# coding: utf-8

# # RNA Velocity analysis using velocyto, Seurat, and scVelo (Part 3: RNA Velocity)
# ## After Part 1: velocyto and Part 2: Seurat extract metadata
# ### Kat Beigel (beigelk@chop.edu)

# See the following pages for guides on how to combine Seurat and velocyto data for scVelo analysis and visualization.
# https://scvelo.readthedocs.io/en/stable/VelocityBasics/
# https://smorabit.github.io/tutorials/8_velocyto/

import anndata
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import loompy as lpy
from scipy.stats import rankdata


# print(np.__version__) # there were errors using newer numpy versions. numpy v1.23.5 worked for me.
# print(sc.__version__)
# print(scv.__version__)
# print(anndata.__version__)


# Load the Seurat data object
adata_seurat = sc.read_h5ad('WT_seurat.h5ad')
print(adata_seurat)

# Load the results .loom file from velocyto
adata_loom = anndata.read_loom('velocytyo_possorted_genome_bam_WT.loom')
print(adata_loom) # to see the object
print(adata_loom.obs.index) # to see the barcodes


# Optional: Add prefixes to barcodes (barcodes between adata_seurat and adata_loom need be the same)
# In my case the adata_seurat had 'WT_' prefixing the barcodes so I added that to the adata_loom
barcodes = [('WT_'+(bc.split(':')[1]).replace('x', '')) for bc in adata_loom.obs.index.tolist()]
adata_loom.obs.index = barcodes
print(adata_loom.obs.index) # to see the barcodes


# Filter adata_loom to keep only observations in adata_seurat (cells kept after filtering in seurat)
adata_loom = adata_loom[np.isin(adata_loom.obs.index, adata_seurat.obs.index)].copy()
print(adata_loom)

# Merge the adata_loom and adata_seurat
adata_all = scv.utils.merge(adata_loom, adata_seurat, copy = True)
print(adata_all.obs)

# Check to see that the UMAP coordinates are present
print(adata_all.obsm['X_umap'])

# Need to make the metadata of clusters (or other classification data) into type 'category'
# This was required for properly plotting and running in scVelo
for col in ['orig.ident', 'seurat_clusters', 'cell_type']:
    adata_all.obs[col] = adata_all.obs[col].astype('category')

adata_all.obs['cell_type'] # check that the cluster names are present


# Need to build palettes for plotting the UMAP with specific cluster colors
# (to match the Seurat colors that I had)
# There is probably a more elegant way to do this but this was my basic solution
# Make a list of each ID type you want to plot
# Cell types
cell_type_list = [
   "celltype_1",
   "celltype_2",
   "celltype_3",
   "celltype_4",
   "celltype_5"]

# Make a palette for each labeling scheme that you want to plot
# palette for cell types
palette_cols_celltypes = {}
for key in cell_type_list:
    for value in color_list:
        palette_cols_celltypes[key] = value
        color_list.remove(value)
        break

# print(palette_cols_celltypes)


# ### RNA Velocity
# https://scvelo.readthedocs.io/en/stable/VelocityBasics.html
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')


# Plot UMAPs to see if eveyrhting looks right
plt.rcParams.update({'font.size': 5})
sc.pl.umap(adata_all, color='cell_type', frameon=False, legend_loc='on data',
           title='Wild-type', palette=palette_cols_celltypes,
           save='WT_umap.pdf')

# Proportions plot (to see proportions of spliced/unspliced)
scv.pl.proportions(adata_all, fontsize=8, figsize=(20, 5), dpi=(300), groupby='cell_type',
                   save='proportions.pdf')

# Preprocess the data
# Whether to start with filter and normalize may depend on your input data, see:
# https://github.com/theislab/scvelo/issues/1052
scv.pp.filter_and_normalize(adata_all)
scv.pp.neighbors(adata_all)
scv.pp.moments(adata_all)
scv.tl.recover_dynamics(adata_all) # if running dynamical model

# Estimate RNA velocity
scv.tl.velocity(adata_all, mode = "dynamical") # make sure to do scv.tl.recover_dynamics() before this
scv.tl.velocity_graph(adata_all)


# Project the velocities
scv.pl.velocity_embedding(adata_all, basis = 'umap', frameon=False,
                          figsize=[15,15],
                          title='Wild-type',
                          save='WT_velocity_embedding.pdf')

scv.pl.velocity_embedding_grid(adata_all, basis='umap', color='cell_type', scale=0.25,
                               figsize=[15,15], legend_loc='right margin',
                               arrow_size=1.75, size=200, alpha=1,
                               title='Wild-type', palette = palette_cols_celltypes,
                               save='WT_velocity_embedding_grid.pdf')


scv.pl.velocity_embedding_stream(adata_all, basis = 'umap', color='cell_type', legend_fontsize=10,
                                 figsize = [10, 10], arrow_size = 2, legend_loc='right margin',
                                 title='Wild-type', palette = palette_cols_celltypes,
                                 save='WT_velocity_embedding_stream.pdf')



# Dynamical modeling: https://scvelo.readthedocs.io/en/stable/DynamicalModeling/
# The driver genes (top-likelihood genes) show dynamic behavior (high likelihood in dynamic model).
top_genes = adata_all.var['fit_likelihood'].sort_values(ascending=False).index
print(top_genes)


scv.pl.scatter(adata_all, basis=top_genes[:20], ncols=1, color='cell_type', frameon=False,
               palette = palette_cols_celltypes,
               save='WT_top20_likelihood_genes_celltype.pdf')


scv.tl.latent_time(adata_all)
print(adata_all.var['velocity_genes'].sum(), adata_all.n_vars, adata_all.n_obs)


# Latent time based on transcriptional dynamics
scv.pl.scatter(adata_all, color='latent_time', color_map='gnuplot', size=80)
scv.pl.heatmap(adata_all, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100)


