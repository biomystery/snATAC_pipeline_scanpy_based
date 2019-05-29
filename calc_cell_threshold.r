sink(tempfile())
suppressMessages(require(ggpubr))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(mixtools))
suppressMessages(require(grid))

############################################################
## input
############################################################

args <- commandArgs(trailingOnly = TRUE)

sample <-  args[1]; # "PDL1"
lower.th <- ifelse(length(args)>1,as.numeric(args[2]),2) ## default 2
th.prob <- ifelse(length(args)>2,as.numeric(args[3]),.1) ## default 0.1

qc <- fread(paste0(sample,'/',sample, ".qc_metrics.txt")) %>% mutate(log10_uniq_usable_reads = log10(unique_usable_reads +
    1))

set.seed(1)
wait <- qc %>% filter(log10_uniq_usable_reads > lower.th) %>% pull(log10_uniq_usable_reads)
suppressMessages(mixmdl <- normalmixEM(wait, k = 2))

tot_barcodes <- nrow(qc)
tot_usable_reads <- sum(qc$unique_usable_reads)
p0 <- ggplot(qc) + geom_histogram(aes(x = log10_uniq_usable_reads, ..count..), binwidth = 0.05, 
    alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0, 5))
p1 <- p0 + annotation_custom(grobTree(textGrob(paste0("Total barcodes", "\n", tot_barcodes, 
    "\nTotal unique reads:", "\n", round(tot_usable_reads/10^6), "M"), x = 0.9, y = 0.9, 
    hjust = 1, vjust = 1, gp = gpar(size = 4))))

tot_cells <- sum(qc$log10_uniq_usable_reads > lower.th)
tot_reads <- sum(qc %>% filter(log10_uniq_usable_reads > lower.th) %>% pull(unique_usable_reads))
p2 <- p0 %+% (qc %>% filter(log10_uniq_usable_reads > lower.th)) + geom_vline(xintercept = lower.th, 
    color = "black", linetype = 2) + annotation_custom(grobTree(textGrob(paste0("Total barcodes(reads>", 
    round(10^lower.th - 1), "):", "\n", tot_cells, ",", round(tot_cells/tot_barcodes * 
        100), "%", "\nTotal reads:", round(tot_reads/10^6), "M,", round(tot_reads/tot_usable_reads * 
        100), "%"), x = 0.05, y = 0.9, hjust = 0, vjust = 1, gp = gpar(size = 4))))
## 
plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
}

threshold = qnorm(th.prob, mixmdl$mu[2], mixmdl$sigma[2])
final_cells <- sum(qc$log10_uniq_usable_reads > threshold)
final_reads <- sum(qc %>% filter(log10_uniq_usable_reads > threshold) %>% pull(unique_usable_reads))


p3 <- ggplot(qc %>% filter(log10_uniq_usable_reads > lower.th)) + geom_histogram(aes(x = log10_uniq_usable_reads, 
    ..density..), binwidth = 0.05, alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0, 
    5)) + stat_function(geom = "line", fun = plot_mix_comps, args = list(mixmdl$mu[1], 
    mixmdl$sigma[1], lam = mixmdl$lambda[1]), colour = "red", lwd = 1) + stat_function(geom = "line", 
    fun = plot_mix_comps, args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]), 
    colour = "green", lwd = 1) + geom_vline(xintercept = mixmdl$mu, color = c("red", 
    "green")) + geom_vline(xintercept = threshold, color = "green", linetype = 2) + 
    annotation_custom(grobTree(textGrob(paste0("log likelihood:", round(mixmdl$loglik), 
        "\n", "G1.median=", round(10^mixmdl$mu[1] - 1), "\n", "G2.median=", round(10^mixmdl$mu[2] - 
            1), "\n", "threshold=", round(10^threshold - 1), ",prob=", round(th.prob, 
            2), "\n", "(", final_cells, " cells,", round(final_cells/tot_barcodes * 
            100), "%)", "\n(", round(final_reads/10^6), "M reads,", round(final_reads/tot_usable_reads * 
            100), "%)"), x = 0.05, y = 0.8, hjust = 0, vjust = 1, gp = gpar(size = 4))))


pdf(paste0(sample,"_cell_threshold.pdf"))
ggarrange(p1, p2, p3, nrow = 3)
dev.off()
sink()
cat(sample,'\t',round(10^(threshold)-1),"\n")
