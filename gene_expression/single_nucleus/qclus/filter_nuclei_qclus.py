import os
import argparse

import numpy as np
import pandas as pd
import qclus as qc
import scanpy as sc 
import matplotlib.pyplot as plt

def plot_qc_metrics(adata, outfile, show_initial_filter=True):
	"""
	Modified from qclus vignette
	"""

	# Define the desired order of qclus categories
	desired_order = ['passed', 'scrublet filter', 'outlier filter', 'clustering filter', 'initial filter']
	if not show_initial_filter:
		adata = adata[adata.obs['qclus'].isin(['passed', 'scrublet filter', 'outlier filter', 'clustering filter'])]

	# Get the unique values present in the qclus column
	unique_qclus_values = adata.obs['qclus'].unique().tolist()

	# Filter the desired order to only include categories that are present in the data
	qclus_order = [category for category in desired_order if category in unique_qclus_values]

	# Ensure the qclus column is a categorical type with the specified order
	adata.obs['qclus'] = adata.obs['qclus'].astype('category')
	adata.obs['qclus'] = adata.obs['qclus'].cat.reorder_categories(qclus_order, ordered=True)

	# Create the figure and subplots
	fig, axs = plt.subplots(2, 2, figsize=(15, 10))

	# Plot each subplot
	sc.pl.violin(adata, 'fraction_unspliced', palette='Pastel1', stripplot=False, inner='box', groupby='qclus', order=qclus_order, ax=axs[0, 0], show=False)
	axs[0, 0].set_title('Fraction Unspliced')

	sc.pl.violin(adata, 'pct_counts_MT', palette='Pastel1', stripplot=False, inner='box', groupby='qclus', order=qclus_order, ax=axs[0, 1], show=False)
	axs[0, 1].set_title('Pct Counts MT')

	sc.pl.violin(adata, 'pct_counts_nuclear', palette='Pastel1', stripplot=False, inner='box', groupby='qclus', order=qclus_order, ax=axs[1, 0], show=False)
	axs[1, 0].set_title('Pct Counts Nuclear')

	adata.obs['log10_total_counts'] = np.log10(adata.obs['total_counts']+1)

	sc.pl.violin(adata, 'log10_total_counts', palette='Pastel1', stripplot=False, inner='box', groupby='qclus', order=qclus_order, ax=axs[1, 1], show=False)
	axs[1, 1].set_title('Log10 Total Counts')

	# Adjust layout and show the plot
	plt.tight_layout()
	plt.savefig(outfile)
	plt.close('all')


def plot_scanpy_umap(adata, outfile):

	adata = adata[adata.obs['qclus'].isin(['passed', 'scrublet filter', 'outlier filter', 'clustering filter'])]

	#run standard processing for visualization purposes
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)

	#save raw dataset and filter genes
	adata.raw = adata
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	sc.pp.filter_genes(adata, min_cells=10)
	adata = adata[:, adata.var.highly_variable]

	sc.pp.regress_out(adata, ['total_counts', 'pct_counts_MT'], n_jobs = 4)
	sc.pp.scale(adata, max_value=10)

	sc.tl.pca(adata, svd_solver='randomized')
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
	sc.tl.umap(adata)
	sc.tl.leiden(adata, key_added="leiden")

	sc.pl.umap(adata, color=["fraction_unspliced", "pct_counts_MT", "pct_counts_nuclear", "qclus"], ncols=2, show=False)
	plt.tight_layout()
	plt.savefig(outfile)
	plt.close('all')


if __name__ == '__main__':

	PARSER = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	PARSER.add_argument('-i', '--cellranger_dir', required=True)

	PARSER.add_argument('-u', '--unspliced', required=True)

	PARSER.add_argument('--use_matrix', required=False, default='raw')

	PARSER.add_argument('--min_genes', required=False, default=215)

	PARSER.add_argument('-o', '--outdir', required=False, default='qclus_out')
	
	ARGS = vars(PARSER.parse_args())

	os.makedirs(ARGS['outdir'], exist_ok=True)

	assert ARGS["use_matrix"] in ['raw', 'filtered']

	adata_path = ARGS['cellranger_dir']+"/"+ARGS["use_matrix"]+"_feature_bc_matrix.h5"

	fraction_unspliced = pd.read_csv(ARGS['unspliced'], index_col=0)

	print('Running qclus')

	adata = qc.run_qclus(adata_path, fraction_unspliced, minimum_genes=ARGS['min_genes'], clustering_features=['pct_counts_nuclear', 'pct_counts_MT', 'fraction_unspliced'], clustering_k=3, clusters_to_select=["0", "1"])

	print('Plotting QC metrics')

	plot_qc_metrics(adata, ARGS['outdir']+'/qclus_metrics.pdf')

	plot_qc_metrics(adata, ARGS['outdir']+'/qclus_metrics_noinit.pdf', show_initial_filter=False)

	print('Saving qclus results to file')

	adata.obs['qclus'].to_csv(ARGS['outdir']+'/qclus_barcodes_res.tsv', sep='\t') 

	adata.obs[['fraction_unspliced', 'pct_counts_MT', 'pct_counts_nuclear', 'total_counts', 'qclus']].to_csv(ARGS['outdir']+'/qclus_metrics_all.tsv', sep='\t') 

	print('Plotting scanpy UMAP')

	plot_scanpy_umap(adata, ARGS['outdir']+'/scanpy_umap_qclus.png')


