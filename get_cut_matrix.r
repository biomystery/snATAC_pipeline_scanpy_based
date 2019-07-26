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
## functions
############################################################

do.skip.step <- function(input_files, output_files) {
    isInputExists <- all(sapply(input_files, file.exists))
    isOuputExists <- all(sapply(output_files, file.exists))
    isInputNewer <- any(sapply(input_files, function(i) any(sapply(output_files, function(o) file.mtime(i) >file.mtime(o)))))

    if (!isInputExists)
        stop(paste("One or more of ", paste(input_files, collapse = ","), " does not exist"))
    if (!isOuputExists)
        return(F)
    else if(isInputNewer){ return(F)
    }else{
        return(T)
    }
}

############################################################
## main
############################################################
## input/output
prefix<-sub('.bed$','',basename(region.file))
output.uniq.region <- file.path(output.dir,paste0(prefix,'.uniq.bed'))
output.intersect.res <- file.path(output.dir,paste0(prefix,'.intersect.tsv'))
output.feature.sp <-file.path(output.dir,paste0(prefix,'.intersect.mtx'))
output.feature.sp.x <-file.path(output.dir,paste0(prefix,'.intersect.regions'))
output.feature.sp.y <-file.path(output.dir,paste0(prefix,'.intersect.cells'))
## make unqiue regions/features

if(!do.skip.step(input_files=region.file,output_files=output.uniq.region)){
    suppressMessages(require(data.table))
    suppressMessages(require(tidyverse))
    sprintf("Making features uniq")
    regions <- fread(region.file)%>%mutate(V4=make.names(V4,unique = T))
    fwrite(regions,output.uniq.region,sep = '\t',col.names = F)
}else{sprintf("Skipped making features uniq")}



## intersect features with tags/cuts: |feature_name|cell|dist_to_edge|
 if(!do.skip.step(input_files=c(output.uniq.region,tag.file),
                  output_files=output.intersect.res)){
     sprintf('bedIntersecting...')
     cmd = paste0("zcat ", tag.file, "| awk -v OFS='\t' '{if($6==\"+\") print $1,$2,$2+1,$4; else print $1,$3-1,$3,$4}'")
     cmd = paste0(cmd, "|intersectBed -a ", output.uniq.region, " -b - -wa -wb |awk '{{print $4,$10,$8-$2}}'>",
             output.intersect.res)
    print(cmd)
     system.time(system(cmd))}else{sprintf("Skipped bedIntersecting")}

if(!do.skip.step(input_files=output.intersect.res,output_files=c(output.feature.sp,
                                                                 output.feature.sp.x,
                                                                 output.feature.sp.y))){
    sprintf('converting to sparse count matrix...')
    suppressMessages(require(data.table))
    suppressMessages(require(tidyverse))
    suppressMessages(require(Matrix))
    res.intersect <- fread(output.intersect.res,drop=3,col.names = c('feature.name','tag.name'),stringsAsFactors = T)
    res.intersect<-res.intersect[,.N,by=.(feature.name,tag.name)]
    m <- with(res.intersect, sparseMatrix(i = as.numeric(feature.name),j = as.numeric(tag.name),
                                      x = N, dimnames = list( levels(feature.name),levels(tag.name))))
    t <- writeMM(m,output.feature.sp)
    write.table(data.frame(rownames(m)),file=output.feature.sp.x, col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(data.frame(colnames(m)),file=output.feature.sp.y, col.names=FALSE, row.names=FALSE, quote=FALSE)
}else{sprintf('Sparse matrix already exists. Done')}
sprintf('Finished')
