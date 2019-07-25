#!/usr/bin/env Rscript
suppressMessages(require(optparse))

############################################################
## input
############################################################
option_list = list(
    make_option(c("-r", "--region_file"), type="character",
                help="named region file in bed format (6 columns) "),
    make_option(c("-t", "--tagalign_file"), type="character",
                help="tn5 shifted tagAlign file (.gz) format "),
    make_option(c("-o", "--output_dir"), type="character",default='./',
                help="output dir [default %default]", metavar="character")

);

## check args
opt = parse_args( OptionParser(option_list=option_list))
region.file<-  opt$region_file
tag.file <- opt$tagalign_file
output.dir <- opt$output_dir

for (f in c(region.file,tag.file,output.dir))
    if( file.access(f) == -1) stop(sprintf("%s  does not exist", f))

## check if 6 column bed here(TODO)
if(!grepl('[.]bed$',basename(region.file))) stop(sprintf("%s  is not .bed", region.file))

############################################################
## main
############################################################
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(rtracklayer))
suppressMessages(require(Matrix))

prefix<-sub('.bed$','',basename(region.file))
output.uniq.region <- file.path(output.dir,paste0(prefix,'.uniq.bed'))
output.intersect.res <- file.path(output.dir,paste0(prefix,'.intersect.tsv'))
output.feature.sp <-file.path(output.dir,paste0(prefix,'.intersect.mtx'))

## make unqiue regions/features
regions <- fread(region.file)%>%mutate(V4=make.names(V4,unique = T))
fwrite(regions,output.uniq.region,sep = '\t',col.names = F)


## intersect features with tags/cuts: |feature_name|cell|dist_to_edge|
cmd = paste0("zcat ", tag.file, "| awk -v OFS='\t' '{if($6==\"+\") print $1,$2,$2+1,$4; else print $1,$3-1,$3,$4}'")
cmd = paste0(cmd, "|intersectBed -a ", output.uniq.region, " -b - -wa -wb |awk '{{print $4,$10,$8-$2}}'>",
             output.intersect.res)
print(cmd)
system.time(system(cmd))

res.intersect <- fread(output.intersect.res,drop=3,col.names = c('feature.name','tag.name'),stringsAsFactors = T)
res.intersect<-res.intersect[,.N,by=.(feature.name,tag.name)]


## save to parse matrix (row-cells,colun - features)
m <- with(res.intersect, sparseMatrix(i = as.numeric(tag.name), j = as.numeric(feature.name),
                                      x = N, dimnames = list( levels(tag.name),levels(feature.name))))

## DONE
