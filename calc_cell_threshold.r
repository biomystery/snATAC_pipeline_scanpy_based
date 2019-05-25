sink(tempfile())
suppressMessages(require(ggpubr))
suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(mixtools))

############################################################
## input
############################################################

args <- commandArgs(trailingOnly = TRUE)

sample <-  args[1]; # "PDL1"
lower.th <- ifelse(length(args)>1,as.numeric(args[2]),2) ## default 2
th.prob <- ifelse(length(args)>2,as.numeric(args[3]),.05) ## default 2

qc <- fread(paste0(sample,'/',sample, ".qc_metrics.txt")) %>% mutate(log10_uniq_usable_reads = log10(unique_usable_reads +
    1))

set.seed(1)
wait <- qc %>% filter(log10_uniq_usable_reads > lower.th) %>% pull(log10_uniq_usable_reads)
suppressMessages(mixmdl <- normalmixEM(wait, k = 2))

p1 <- ggplot(qc) + geom_histogram(aes(x = log10_uniq_usable_reads, ..count..), binwidth = 0.05,
    alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0, 5))
p2 <- p1 %+% (qc %>% filter(log10_uniq_usable_reads > lower.th))
##
plot_mix_comps <- function(x, mu, sigma, lam) {
    lam * dnorm(x, mu, sigma)
}

threshold = qnorm(th.prob, mixmdl$mu[2], mixmdl$sigma[2])

tot_cells <- sum(qc$log10_uniq_usable_reads > lower.th)
final_cells <- sum(qc$log10_uniq_usable_reads > threshold)

p3 <- ggplot(qc %>% filter(log10_uniq_usable_reads > lower.th)) + geom_histogram(aes(x = log10_uniq_usable_reads,
    ..density..), binwidth = 0.05, alpha = 0.5, color = "black") + theme_bw() + coord_cartesian(xlim = c(0,
    5)) + stat_function(geom = "line", fun = plot_mix_comps, args = list(mixmdl$mu[1],
    mixmdl$sigma[1], lam = mixmdl$lambda[1]), colour = "red", lwd = 1) + stat_function(geom = "line",
    fun = plot_mix_comps, args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
    colour = "green", lwd = 1) + geom_vline(xintercept = mixmdl$mu, color = c("red",
    "green")) + geom_vline(xintercept = threshold, color = "green", linetype = 2) +
    annotate("text", 0, 0.75, size = 3, hjust = 0, label = paste0("total barcodes(>",round(10^lower.th-1),"):",
        "\n",tot_cells, "\n", "G1.median=", round(10^mixmdl$mu[1]-1), "\n", "G2.median=",
        round(10^mixmdl$mu[2]-1), "\n", "threshold=", round(10^threshold-1),",prob=",round(th.prob,2), "\n", "(",
        final_cells, " cells,", round(final_cells/tot_cells * 100), "%)"))

pdf(paste0(sample,"_cell_threshold.pdf"))
ggarrange(p1, p2, p3, nrow = 3)
dev.off()
sink()
cat(round(10^(threshold)-1),"\n")
