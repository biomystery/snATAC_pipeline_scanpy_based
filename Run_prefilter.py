#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import papermill as pm
import sys
import subprocess

sample_name = sys.argv[1]
threshold = sys.argv[2]
print("Running jupyter:" + sample_name)
pm.execute_notebook(
    'prefilter.ipynb',
    '{0}/{0}_prefilter_{1}.ipynb'.format(sample_name,threshold),
    parameters={"sample": sample_name,
                "lower.th": threshold
                }
)

subprocess.call(["jupyter-nbconvert", "--to", "html",'{0}/{0}_prefilter_{1}.ipynb'.format(sample_name,threshold)])
# In[ ]:
