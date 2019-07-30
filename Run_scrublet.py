#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import papermill as pm
import sys
import subprocess

sample_name = sys.argv[1]
print("Running jupyter:" + sample_name)
pm.execute_notebook(
    'run_scrublet.ipynb',
    '{0}/{0}_scrublet.ipynb'.format(sample_name),
    parameters={"sample_name": sample_name,
                }
)

subprocess.call(["jupyter-nbconvert", "--to", "html",'{0}/{0}_scrublet.ipynb'.format(sample_name)])
# In[ ]:
