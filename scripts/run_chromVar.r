############################################################
## script to run ChromVar
############################################################
## Input: 1. peaks 2. reads in bed or tagalign, 3. cell type label
## Parameters:
## Output:


##------------------------------------------------------------
## inputs
##------------------------------------------------------------
input.peaks <- getPeaks('./Islet_123.combined.merged.bed',sort_peaks = TRUE)
input.peaks <- resize(input.peaks,width=500,fix="center")
input.peaks <- as.data.table(input.peaks)[,.(seqnames,start,end)]
setkey(input.peaks,seqnames,start,end)

input.filtered.reads <-  fread('./output.filtered.reads.bed')
input.filtered.reads <- input.filtered.reads[,.(V1,V2,V3,V4)]
input.filtered.reads <- input.filtered.reads[,.(seqnames=V1,start=V2,end=V3,cell=V4)]

input.cells.final <- fread('./output.cell.celltypes.txt')



##------------------------------------------------------------
## Main
##------------------------------------------------------------
## Step 1: cut counts in - peaks x cells
# https://github.com/biomystery/islet_figs_motifs/blob/master/chromVar/count_reads_in_peaks.R



## Step 2: create SummerizedExperiment Objects for chromVar: 1. filter samples,2. filter peaks. 3. add GC 4. add cell_type assignemnt
# https://github.com/biomystery/islet_figs_motifs/blob/master/chromVar/creat_count_object_for_chromVAR.R


## Step 3:
