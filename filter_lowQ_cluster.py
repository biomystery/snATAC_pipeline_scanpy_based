#!/usr/bin/env python
# coding: utf-8

# import
import os
import sys
import scanpy.api as sc
from anndata import AnnData

sample_name = sys.argv[1]
cluster_to_filter = sys.argv[2]

print('filtering clusters:{0} for {1} '.format(cluster_to_filter, sample_name))

wd = os.path.join(os.getcwd(), sample_name)
adata = sc.read_h5ad(filename=os.path.join(
    wd, '{0}.adata.h5ad'.format(sample_name)))


barcodes = adata.obs.index[adata.obs.leiden.isin([cluster_to_filter])]
filename = os.path.join(wd, '{}.filtered.txt'.format(sample_name))

print('total  {0} low quality cells will be recorded in {1}'.format(
    len(barcodes), filename))

with open(filename, 'w') as f:
    f.writelines('\n'.join(barcodes))

print('END')
