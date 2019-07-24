#!/usr/bin/env python

import os
import pandas as pd
from anndata import AnnData
import scanpy.api as sc
import argparse


def main():
    '''
    Given input adata and output umap and clustering results
    '''

    # parse_args
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', '--adata_file', required=True,
                        help='Andata file (h5ad) format')
    parser.add_argument('-o', '--output_file', required=True,
                        help='output file (csv) format')

    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata_file)
    adata.obs['umap_x'] = adata.obsm.X_umap[:, 0].tolist()
    adata.obs['umap_y'] = adata.obsm.X_umap[:, 1].tolist()

    adata.obs[['umap_x', 'umap_y', 'leiden']].to_csv(args.output_file)


if __name__ == '__main__':
    main()