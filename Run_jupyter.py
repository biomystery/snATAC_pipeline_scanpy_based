#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import papermill as pm
import sys

sample_name = sys.argv[1]
threshold = sys.argv[2]
print("Running jupyter:" + sample_name)
pm.execute_notebook(
    'clustering_single.ipynb',
    '{0}/{0}.ipynb'.format(sample_name),
    parameters={"sample_name": sample_name,
                "min_usable_reads": threshold
                }
)


# In[ ]:
