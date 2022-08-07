library(ggplot2)


args <- commandArgs(TRUE)


geneID <- args[1]
chrID <- args[2]
pt1 <- as.numeric(args[3])
pt2 <- as.numeric(args[4])

#
#data1 <- read.table(paste0(geneID,'/',chrID,'_',pt1,'_',pt2,'_',geneID,'.txt'), header=T)
data1 <- read.table(paste0(geneID,'/',geneID,'.pval.txt'), header=T)

data1$adjP <- p.adjust(data1$pval, method='BH')
data1$logP <- -log10(data1$pval)
data1$logQ <- -log10(data1$adjP)

#sigDF <- data1[which(data1$pval < 0.05),]
#nonSigDF <- data1[which(data1$pval >= 0.05),]
sigDF <- data1[which(data1$adjP < 0.05),]
nonSigDF <- data1[which(data1$adjP >= 0.05),]

##
p1 <- ggplot() +
    #geom_point(data=nonSigDF, aes(x=genomicPt, y=logP), size=1, color='black', shape=1) +
    #geom_point(data=sigDF, aes(x=genomicPt, y=logP), size=1.5, color='red', shape=1) +
    geom_point(data=nonSigDF, aes(x=genomicPt, y=logQ), size=1, color='black', shape=1) +
    geom_point(data=sigDF, aes(x=genomicPt, y=logQ), size=1.5, color='red', shape=1) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = 'black', cex=1, alpha=0.5) +
    geom_vline(xintercept=c(pt1, pt2), color = 'black', cex=1) +
    xlim(as.numeric(pt1), as.numeric(pt2)) +
    #labs(x='', y='-log10(P-value)') +
    labs(x='', y='-log10(Q-value)') +
    theme_classic() 

##
pdf(paste0(geneID,'/',chrID,'_',pt1,'_',pt2,'_',geneID,'.pdf'), width=9, height=3, pointsize=3)
plot(p1)
dev.off()



