#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import papermill as pm
import sys
import subprocess
import pandas as pd

sample_name = sys.argv[1]

# Read thresholds from file
ths = pd.read_csv('{0}/{0}.qc_thresholds.txt'.format(sample_name), sep=',')
print(ths)

print("Running jupyter:" + sample_name)
pm.execute_notebook(
    'clustering_single.ipynb',
    '{0}/{0}.ipynb'.format(sample_name),
    parameters={"sample_name": sample_name,
                "min_usable_reads": str(ths['min_total'][0]),
                "min_frop": str(ths['min_FRoP'][0]),
                }
)

subprocess.call(["jupyter-nbconvert", "--to", "html",
                 '{0}/{0}.ipynb'.format(sample_name)])
# In[ ]:
