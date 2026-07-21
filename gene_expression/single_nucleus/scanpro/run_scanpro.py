from scanpro import scanpro
import scanpy as sc
import pandas as pd

data = sc.read_10x_mtx('../expression_data/single-nucleus/final/Ni2_sparse_matrix/')

cluster_assignments = pd.read_csv('../expression_data/single-nucleus/final/Ni2_cell_clusters.tsv', sep='\t', index_col=0, usecols=[0, 1])

orig = pd.read_csv('../expression_data/single-nucleus/final/Ni2_cell_orig.tsv', sep='\t', index_col=0, usecols=[0, 1])

data.obs = data.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)

data.obs = data.obs.merge(orig, how='left', left_index=True, right_index=True)

out = scanpro(data, clusters_col='cluster', conds_col='sample', transform='arcsin')

out.results.to_csv('Ni2_scanpro_arcsin_default-5rep-100sim.tsv', sep='\t')