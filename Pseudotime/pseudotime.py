# python
## pseudotime analysis applied for the metadata analysis

# import package
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as pl
from matplotlib import rcParams

#setup working direction
os.chdir('/Users/hongyuliu/Desktop/work/meta/')

# read sample sheet
meta_data = sc.read_excel('meta_d12.xlsx', sheet ='Sheet2')

# locate the dataset
sc.tl.pca(a,svd_solver='arpack')
sc.pp.neighbors(meta_data, n_pcs=30)
sc.tl.umap(meta_data)
sc.tl.draw_graph(meta_data)

# non-liner regression
sc.tl.louvain(meta_data, resolution=1.0)
sc.tl.leiden(meta_data)
sc.pl.paga(meta_data, color=['leiden'])
sc.pl.draw_graph(meta_data, color = ['0','2','5'])
sc.pl.draw_graph(meta_data, color = ['0','2','5','10','15','20'])

# paga analysis by timepoint
sc.pl.paga(meta_data, color = ['Glutamate Metabolism'])
sc.tl.draw_graph(meta_data, init_pos='paga')
sc.pl.draw_graph(meta_data, color = ['1_2','2_2','3_2','4_2','2_2'], layout = 'fa')
sc.pl.paga_compare(
    meta_data, threshold=0.03, title='', right_margin=0.2, size=100, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True)
sc.pl.paga(meta_data, color=['louvain'], save = 'line_link_louvain_d12_1_2.png')

# setup start point
meta_data.uns['iroot'] = np.flatnonzero(meta_data.obs['leiden']  == '5')[0]
sc.pl.umap(meta_data, color = ['leiden', 'dpt_pseudotime'])
sc.pl.draw_graph(meta_data, color=['leiden', 'dpt_pseudotime'], legend_loc='on data', save = 'meta_lov_pso_d12_1_2.png')
