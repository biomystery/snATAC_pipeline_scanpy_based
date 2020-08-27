                                        # require(ggplot2)

require(BSgenome.Mmusculus.UCSC.mm10)
require(ChIPseeker)
require(circlize)
require(ComplexHeatmap )
require(data.table)
require(extrafont)  # fonts
require(GenomicRanges)
require(ggpubr)
require(ggrepel)
require(ggthemes)
require(Matrix)
require(pracma)#tic;toc
require(RColorBrewer)
require(SummarizedExperiment)
require(tidyverse)
require("TxDb.Mmusculus.UCSC.mm10.knownGene")
## https://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
suppressMessages(loadfonts())

theme_pubr <- function(base_size =  11,  base_family =  "") {
    ## , base_family = 'Arial'
    theme_foundation() + theme(line =  element_line(colour =  "black",  lineend =  "round",linetype =  "solid"),
                               rect =  element_rect(fill =  "white",  colour =  "black",linetype =  "solid"),
                               text =  element_text(colour =  "black",  face =  "plain",family =  base_family,
                                                    size =  base_size, vjust =  0.5, hjust =  0.5, lineheight =  0.5),
                               panel.background =  element_blank(), plot.background =  element_blank(), panel.border =  element_rect(colour =  "black",fill =  NA),
                               panel.grid =  element_blank(), strip.background =  element_rect(colour =  NA),legend.key =  element_rect(colour =  NA), title =  element_text(size =  rel(1)),
                                  plot.title =  element_text(size =  rel(1.2),  face =  "bold"), strip.text =  element_text(),axis.ticks.length =  unit(1,  "mm"))
}
"%ni%" <- Negate( "%in%")
objSize <-  function(obj,unit='Gb') object.size(obj) %>% format(unit)
format.bignumber <-  function(x) format(x,big.mark = ',')
trimRange <- function(pd,  x) {
    pd[pd > x] <- x
    pd[pd < -x] <- -x
    return(pd)
}

pal_motif <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)


ccc_fun <- function(x, y) {

    ## calculate lin's concordance correlation coefficient
    ##https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    ##x <- rand(100,1); y<- rand(100,1)
    2*cor(x, y) *sqrt(var(x))*sqrt(var(y))/(var(x)+var(y)+(mean(x)-mean(y))^2)
}



countInsertions <- function(query,  fragments,  by =  "RG") {

    ## Count By Fragments Insertions
    message("reading tags")
    tic()
    inserts <- fread(fragments,  col.names =  c("chr",  "start",  "end",  "RG",  "score",
                                                "strand")) %>%  dplyr::select(-score) %>%  mutate(start =  start + 1)
    a <-  toc()
    message("overlaping")
    inserts <- GRanges(seqnames =  inserts$chr,  ranges =  IRanges(start =  inserts$start,
                                                                   end =  inserts$end), strand =  inserts$strand, RG =  inserts$RG) %>%  resize(width =  1,
                                                                                                                                               fix =  "start")
    overlapDF <- DataFrame(findOverlaps(query,  inserts,  ignore.strand =  TRUE,  maxgap =  -1L,
                                        minoverlap =  0L, type =  "any"))
    overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
    overlapTDF <- transform(overlapDF,  id =  match(name,  unique(name)))
                                        # Calculate Overlap Stats
    inPeaks <- table(overlapDF$name)
    total <- table(mcols(inserts)[,  by])
    total <- total[names(inPeaks)]
    frip <- inPeaks/total

    ## Summarize
    sparseM <- Matrix::sparseMatrix(i =  overlapTDF[,  1],
                                    j =  overlapTDF[,  4],
                                    x =  rep(1,nrow(overlapTDF)), dims =  c(length(query),  length(unique(overlapDF$name))))
    colnames(sparseM) <- unique(overlapDF$name)
    total <- total[colnames(sparseM)]
    frip <- frip[colnames(sparseM)]
    out <- list(counts =  sparseM,  frip =  frip,  total =  total)
    toc() - a
    return(out)
}

## https://support.bioconductor.org/p/78652/
extend <- function(x,  upstream = 0,  downstream = 0)
    {

        if (any(strand(x) ==  "*"))
            warning("'*' ranges were treated as '+'")
        on <- plus <- strand(x) <<-  "+" | strand(x) <<-  "*"
        new <- start <- start(x) - ifelse(on <- plus,  upstream,  downstream)
        new <- end <- end(x) + ifelse(on <- plus,  downstream,  upstream)
        ranges(x) <- IRanges(new <- start,  new <- end)
        trim(x)
    }

grexpand <- function(inputGR,  upstream = 0,  downstream = 0) {

    inputGR <- resize(inputGR,   width =  width(inputGR)+upstream,  fix =  'end')
    inputGR <- resize(inputGR,   width =  width(inputGR)+downstream,  fix =  'start')
    return(inputGR);
    }

