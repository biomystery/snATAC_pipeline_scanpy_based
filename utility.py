import os
import gzip
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import math
import seaborn as sns
import umap
import re
import statsmodels.api as sm
import sklearn.preprocessing
import scipy
import scipy.sparse
import sklearn.metrics
import sklearn.mixture
import sklearn.linear_model
import subprocess

from anndata import AnnData
sc.settings.set_figure_params(dpi=100)
# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
# %load_ext rpy2.ipython
sns.set(font_scale=1.5)
plt.style.use('seaborn-white')


def plot_QCmatrix(adata_input, clustering='leiden'):
    fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
    i = 0
    to_plot = [
        'log10_unique_usable_reads', 'frac_reads_in_peaks',
        'frac_reads_in_promoters', 'frac_promoters_used', 'frac_duplicated_reads',
        'doublet_quantile'
    ]

    for ax in axs.reshape(-1):
        sns.boxplot(x=clustering, y=to_plot[i], data=adata_input.obs, ax=ax)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title(to_plot[i])
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        i += 1
    plt.tight_layout()
    plt.show()

    fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
    i = 0

    for ax in axs.reshape(-1):
        sc.pl.umap(adata_input,
                   color=to_plot[i],
                   cmap='coolwarm',
                   size=9,
                   ax=ax,
                   show=False,
                   legend_loc='on data')
        i += 1

    plt.tight_layout()
    plt.show()


def plot_composation(adata_input, clustering='leiden'):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # 1
    sc.pl.umap(adata_input,
               color=clustering,
               size=9,
               legend_loc='on data',
               ax=axs[0, 0],
               show=False)

    # 2
    sc.pl.umap(adata_input, color=['sample_name'],
               alpha=.15, size=9, ax=axs[0, 1], show=False)

    # 3
    cell_per_cluster = adata_input.obs.groupby([
        clustering, "sample_name"
    ]).size().reset_index(name="Cells").pivot_table(index=clustering,
                                                    columns='sample_name',
                                                    values='Cells',
                                                    fill_value=0)

    pd = cell_per_cluster  # .apply(lambda x: round(x / x.sum(),3)*100)
    pd.plot(
        kind='bar',
        stacked=True,
        legend=False,
        ax=axs[1, 0],
    )
    axs[1, 0].set_ylabel('Cell count')

    # 4
    pd = cell_per_cluster.apply(
        lambda x: round(x / x.sum(), 3) * 100,
        axis=1,
    )
    pd.plot(kind='bar', stacked=True, legend=False, ax=axs[1, 1])
    axs[1, 1].set_ylabel('Percentage')
    # patches, labels = axs[1,1].get_legend_handles_labels()
    # axs[1,1].legend(patches, labels, loc='upper left', bbox_to_anchor=(1,1))

    plt.tight_layout()
    plt.show()


def plot_scatter(adata_input, clustering='leiden'):

    n_cluster = len(adata_input.obs[clustering].unique())
    n_row = math.ceil(n_cluster / 5)
    fig, axs = plt.subplots(n_row, 5, figsize=(
        10, 2 * n_row), sharex=True, sharey=True)
    j = 0

    for ax in axs.reshape(-1):
        if j >= n_cluster:
            break
        cols = ['red' if i == str(
            j) else 'grey' for i in adata_input.obs[clustering].tolist()]

        adata_input.obs.plot.scatter(x='log10_unique_usable_reads',
                                     y='frac_reads_in_promoters',
                                     s=3,
                                     c=cols,
                                     ax=ax)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_title("c{0}:{1} cells".format(str(j), str(
            cols.count('red'))), fontdict={'fontsize': 12})
        j += 1

    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)
    plt.ylabel('FRoP')
    plt.xlabel('log10_unique_usable_reads')
    plt.tight_layout()
    plt.show()


def plot_MarkerGenes(adata_input, clustering='leiden', marker_genes=None, check_all=False, use_raw_=True):
    sns.set(font_scale=1)
    plt.style.use('default')

    if not marker_genes:
        marker_genes = {
            'Tcell': ["Cd3e.1", "Cd4.1", "Cd8a", 'Tcf7'],
            'Macrophage': ['Adgre1', 'Eif4a1', 'RP23-144N15.4', 'Cx3cr1'],
            'B_cells': [
                "Cd79b",
                "Mzb1",
            ],
            'Epithelial Cell': ['Krt19.1'],
            'Fib': ['Col1a2', 'Col1a2.1', 'Col1a2.2', 'Col1a2.3', 'Col1a2.4', 'Acta2'],
            'Neutrophil':
            ["Ly6g", "Cebpe", "Csf3r", 'Lcn2', 'Ltf', 'S100a8', 'S100a9'],
            'NK': ['Klrc1'],
            # 'Itgax'-> Cd11c
            'DCs': ["Cd209a", "Cd74", "Flt3", "H2-Eb1", 'Itgax'],
            'MHC-II': ["H2-Aa", "H2-Ab1", "H2-Eb1"],  # "
        }

    if check_all and type(marker_genes) == 'dict':
        marker_genes = [x for x in adata_input.raw.var_names if re.sub(r'\.[0-9]+', '', x) in
                        set([m for v in marker_genes.values() for m in v])]
    elif check_all:
        marker_genes = [x for x in adata_input.raw.var_names if re.sub(r'\.[0-9]+', '', x) in
                        set([m for m in marker_genes])]

    ax = sc.pl.matrixplot(adata_input,
                          var_names=marker_genes,
                          cmap='Reds',
                          dendrogram=True,
                          groupby=clustering,
                          use_raw=use_raw_)

    ax = sc.pl.matrixplot(adata_input,
                          var_names=marker_genes,
                          cmap='Reds',
                          dendrogram=True,
                          groupby=clustering,
                          standard_scale='var',
                          use_raw=use_raw_)

    ax = sc.pl.dotplot(adata_input,
                       var_names=marker_genes,
                       groupby=clustering,
                       dendrogram=True,
                       use_raw=use_raw_,
                       expression_cutoff=0)


def binarize_AnnData(adata):
    '''
    binarize adata with log1p raw
    '''

    from scipy.sparse import find, csr_matrix
    adata_binary = adata.copy()
    nnz_inds = adata_binary.raw.X.nonzero()
    b = csr_matrix(
        (np.ones(len(nnz_inds[0])), (nnz_inds[0], nnz_inds[1])),
        shape=adata_binary.raw.X.shape)
    adata_binary.raw = AnnData(b, var=adata.raw.var, obs=adata.obs)
    return(adata_binary)


def run_scrublet(adata, neotic_ratio=.5):
    '''
    '''
    import scrublet as scr
    from scipy.stats import rankdata

    expected_doublet_th = adata.shape[0] / 1000 * .01 * neotic_ratio
    adata_raw = adata.raw.copy()
    adata_raw = adata_raw[:, adata_raw.var.index.isin(
        adata.var_names.tolist())]
    counts_matris_2 = adata_raw.X.expm1()
    del adata_raw
    scrub = scr.Scrublet(
        counts_matris_2, expected_doublet_rate=expected_doublet_th)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        distance_metric='cosine',
        mean_center=False,
        n_prin_comps=50,
        log_transform=True,
        min_gene_variability_pctl=0)
    scrub.plot_histogram()
    predicted_doublets = scrub.call_doublets(threshold=np.quantile(
        doublet_scores, 1 - expected_doublet_th))  # directly call by trheshold
    print('total predicted doublets:', sum(predicted_doublets))
    print('predicted doublets ratio:', sum(
        predicted_doublets) / len(predicted_doublets))
    adata.obs['doublet_score'] = doublet_scores
    adata.obs['doublet'] = predicted_doublets
    adata.obs['doublet_quantile'] = (
        rankdata(doublet_scores) / len(doublet_scores))
    return(adata)
