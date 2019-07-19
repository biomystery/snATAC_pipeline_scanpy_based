#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import papermill as pm
import sys
import subprocess

sample_name = sys.argv[1]
min_usable_reads = sys.argv[2]
min_frop = sys.argv[3]
print("Running jupyter:" + sample_name)
pm.execute_notebook(
    'clustering_single.ipynb',
    '{0}/{0}.ipynb'.format(sample_name),
    parameters={"sample_name": sample_name,
                "min_usable_reads": min_usable_reads,
                "min_frop": min_frop,
                }
)

subprocess.call(["jupyter-nbconvert", "--to", "html",
                 '{0}/{0}.ipynb'.format(sample_name)])
# In[ ]:
