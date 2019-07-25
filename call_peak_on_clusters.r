#!/usr/bin/env Rscript

suppressMessages(require(optparse))
suppressMessages(require(parallel))
############################################################
## input
############################################################

option_list = list(
    make_option(c("-c", "--cluster_file"), type="character",
                help="cluster file in csv format (col1:cell_barcodes,col2:cluster) ", metavar="character"),
    make_option(c("-t", "--tagalign_file"), type="character",
                help="tn5 shifted tagAlign file (.gz) format ", metavar="character")
);


opt = parse_args( OptionParser(option_list=option_list))

cluster.file<-  opt$cluster_file
tag.file <- opt$tagalign_file

sapply(c(cluster.file,tag.file),function(f) if( !file.access(f) == -1) stop(sprintf("Cluster file ( %s ) does not exist", f)))


############################################################
## main
############################################################
sink(tempfile())
suppressMessages(require(ggpubr))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(mixtools))
suppressMessages(require(grid))

qc <- fread(paste0(sample,'/',sample, ".qc_metrics.txt")) %>% mutate(log_uniq_usable_reads = log(unique_usable_reads +
    1,base=log.base))

set.seed(1)
wait <- qc %>% filter(unique_usable_reads > lower.th) %>% pull(log_uniq_usable_reads)
suppressMessages(mixmdl <- normalmixEM(wait, k = 2))

tot_barcodes <- nrow(qc)
tot_usable_reads <- sum(qc$unique_usable_reads)
p0 <- ggplot(qc) + geom_histogram(aes(x = log_uniq_usable_reads, ..count..), binwidth = 0.05,
    alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0, log(10^5,base=log.base)))
p1 <- p0 + annotation_custom(grobTree(textGrob(paste0("Total barcodes", "\n", tot_barcodes,
    "\nTotal unique reads:", "\n", round(tot_usable_reads/10^6), "M"), x = 0.9, y = 0.9,
    hjust = 1, vjust = 1, gp = gpar(size = 4))))+ xlab('')

tot_cells <- sum(qc$unique_usable_reads > lower.th)
tot_reads <- sum(qc %>% filter(unique_usable_reads > lower.th) %>% pull(unique_usable_reads))
p2 <- p0 %+% (qc %>% filter(unique_usable_reads > lower.th)) + geom_vline(xintercept = log(lower.th+1,base=log.base),
    color = "black", linetype = 2) + annotation_custom(grobTree(textGrob(paste0("Total barcodes(reads>",
    lower.th, "):", "\n", tot_cells, ",", round(tot_cells/tot_barcodes *
        100), "%", "\nTotal reads:", round(tot_reads/10^6), "M,", round(tot_reads/tot_usable_reads *
        100), "%"), x = 0.05, y = 0.9, hjust = 0, vjust = 1, gp = gpar(size = 4))))+ xlab('')
##
plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
}

threshold = qnorm(th.prob, mixmdl$mu[2], mixmdl$sigma[2])
final_cells <- sum(qc$log_uniq_usable_reads > threshold)
final_reads <- sum(qc %>% filter(log_uniq_usable_reads > threshold) %>% pull(unique_usable_reads))


p3 <- ggplot(qc %>% filter(unique_usable_reads > lower.th)) + geom_histogram(aes(x = log_uniq_usable_reads,
    ..density..), binwidth = 0.05, alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0,log(10^5,base=log.base)
    )) + stat_function(geom = "line", fun = plot_mix_comps, args = list(mixmdl$mu[1],
    mixmdl$sigma[1], lam = mixmdl$lambda[1]), colour = "red", lwd = 1) + stat_function(geom = "line",
    fun = plot_mix_comps, args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
    colour = "green", lwd = 1) + geom_vline(xintercept = mixmdl$mu, color = c("red",
    "green")) + geom_vline(xintercept = threshold, color = "green", linetype = 2) +
    annotation_custom(grobTree(textGrob(paste0("log likelihood:", round(mixmdl$loglik),
        "\n", "G1.median=", round(log.base^mixmdl$mu[1] - 1), "\n", "G2.median=", round(log.base^mixmdl$mu[2] -
            1), "\n", "threshold=", round(log.base^threshold - 1), ",prob=", round(th.prob,
            2), "\n", "(", final_cells, " cells,", round(final_cells/tot_barcodes *
            100), "%)", "\n(", round(final_reads/10^6), "M reads,", round(final_reads/tot_usable_reads *
            100), "%)"), x = 0.05, y = 0.8, hjust = 0, vjust = 1, gp = gpar(size = 4))))+ xlab(paste0("Usable reads (log",log.base,")"))


############################################################
## output
############################################################
pdf(paste0(sample,"_cell_threshold.pdf"))
ggarrange(p1, p2, p3, nrow = 3)
dev.off()
sink()
cat(sample,'\t',round(log.base^(threshold)-1),"\n")
