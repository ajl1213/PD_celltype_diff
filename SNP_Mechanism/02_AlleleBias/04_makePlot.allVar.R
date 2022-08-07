
library(ggplot2)
library(RColorBrewer)


dir.create('Plots')


data1 <- read.table('AllelicReads.allVar.pval.txt', header=T)
data1$label <- factor(data1$label, levels=c('altEnrich','refEnrich','unchanged'))
col_list <- c(brewer.pal(8, 'Set2')[1],brewer.pal(8, 'Set2')[2],'grey')

print(table(data1$label))

# v1
p1 <- ggplot(data1) +
    geom_point(aes(x=log2(refReads+1), y=log2(altReads+1), color=label), shape=1) +
    geom_smooth(aes(x=log2(refReads+1), y=log2(altReads+1)), method=lm, se=FALSE, col='darkblue', size=1.5, lty=2) +

    scale_color_manual(values=col_list) +    
    labs(x='log2(reads from ref. allele)',y='log2(Reads from alt. allele)') +
    xlim(min(log2(data1$refReads+1), log2(data1$altReads+1)), max(log2(data1$refReads+1), log2(data1$altReads+1))) + 
    ylim(min(log2(data1$refReads+1), log2(data1$altReads+1)), max(log2(data1$refReads+1), log2(data1$altReads+1))) +
    theme_classic()

pdf('Plots/allelicBias.allVar.scatter.v1.pdf')
plot(p1)
dev.off()


# v2
p1 <- ggplot(data1) +
    geom_point(aes(x=log2(allReads+1), y=log2Enrich, color=label), shape=1) +
    geom_smooth(aes(x=log2(allReads+1), y=log2Enrich), se=F, col='darkblue', size=1.5, lty=2) +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    scale_color_manual(values=col_list) +    
    labs(x='log2(Read depth)',y='Log2 enrichment of alt. reads over expectaion') +
    theme_classic()

pdf('Plots/allelicBias.allVar.scatter.v2.pdf')
plot(p1)
dev.off()


#
pdf('Plots/allelicBias.allVar.box.pdf')
boxplot(data1$log2Enrich, las=2, notch=T, outline=T, ylab='Log2 enrichment of alt. reads over expectation', col=col_list[2], xlab='PD GWAS matched')
abline(h=0, col='black', lty=2, lwd=1.5)
dev.off()

# two-sided one-sample T-test
t.test(data1$log2Enrich, mu = 0, alternative = "two.sided")




