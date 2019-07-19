snATAC pipeline 
==========

** Steps 

1. Run 10x cellranger 
2. Run preprocessing pipeline 
3. Run `Run_prefilter.py` to get thresholds 
4. Run `Run_jupyter.py` to do initial clustering
5. Run `Run_scrublet.py` to get scrublet doublet score
6. Run `Run_merge.py` to merge samples 

** Improvements 

1. Get `cut_freq_matrix` for prom x cells 
2. Group frag/reads for (cells in) clusters, then call peaks 
3. Merge frag/reads 
4. 


