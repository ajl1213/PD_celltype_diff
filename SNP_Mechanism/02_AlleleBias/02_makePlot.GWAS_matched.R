
library(ggplot2)
library(RColorBrewer)


dir.create('Plots')


data1 <- read.table('AllelicReads.GWAS_matched.pval.txt', header=T)
data1$label <- factor(data1$label, levels=c('altEnrich','refEnrich','unchanged'))
col_list <- c(brewer.pal(8, 'Set2')[1],brewer.pal(8, 'Set2')[2],'grey')


# 
p1 <- ggplot(data1) +
    geom_point(aes(x=refReads, y=altReads, color=label, fill=label), shape=1, size=1) +
    geom_smooth(aes(x=refReads, y=altReads), method=lm, se=FALSE, col='red', size=1) +

    scale_color_manual(values=col_list) +    
    scale_fill_manual(values=col_list) +    
    labs(x='Reads from ref. allele',y='Reads from alt. allele') +
    xlim(min(data1$refReads, data1$altReads), max(data1$refReads, data1$altReads)) + 
    ylim(min(data1$refReads, data1$altReads), max(data1$refReads, data1$altReads)) +
    theme_classic()

pdf('Plots/allelicBias.GWAS_matched.scatter.pdf')
plot(p1)
dev.off()

print(table(data1$label))


#
pdf('Plots/allelicBias.GWAS_matched.box.pdf')
boxplot(data1$log2Enrich, las=2, notch=T, outline=T, ylab='Log2 enrichment (alt reads over expectation)', col=col_list[2], xlab='PD GWAS matched')
abline(h=0, col='black', lty=2, lwd=1.5)
dev.off()

#
p1 <- ggplot() +
    geom_violin(aes(x=1, y=data1$log2Enrich), trim=F, scale='width', width=1.1) +
    geom_hline(yintercept=0, linetype='dashed', color='black') +
    labs(x='PD GWAS matched heterozygousvariants', y='Log2 enrichment of alt. reads over expectation') +
    theme_classic()

pdf('Plots/allelicBias.GWAS_matched.violin.pdf')
plot(p1)
dev.off()

# two-sided one-sample T-test
t.test(data1$log2Enrich, mu = 0, alternative = "two.sided")




