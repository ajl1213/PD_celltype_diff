
library(ggplot2)
library(RColorBrewer)


dir.create('Plots')


data1 <- read.table('eQTL_HiC.comparison.count.txt', header=T)
celltypes <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')


# hypergeometric test
pvals <- c()
for (i in 1:nrow(data1)){
    p_hyper <- sum(dhyper(data1$common[i]:data1$hic_spec[i], data1$eqtl_spec[i], data1$allAsso[i]-data1$eqtl_spec[i], data1$hic_spec[i]))
    pvals <- c(pvals, p_hyper)   
}

df <- data.frame(data1, pvals)

print(df)

# Down Bar
downDF <- df[which(df$Datatype=='Down'),]
a <- data.frame(Celltype=downDF$Celltype, label='eqtl_spec', count=downDF$eqtl_spec)
b <- data.frame(Celltype=downDF$Celltype, label='common', count=downDF$common)
finalDF <- rbind(a, b)
finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$label <- factor(finalDF$label, levels=c('eqtl_spec','common'))

p1 <- ggplot() + 
    geom_bar(data=finalDF, aes(x=Celltype, y=count, fill=label), stat='identity') +
    scale_fill_manual(values=c('darkgrey', 'darkorange')) +
    labs(x='', y='count') +
    coord_flip() + 
    theme_classic() +
    labs(x='',y='Number of eQTLs', title='Down cRE')

pdf('Plots/Down.bar.pdf')
plot(p1)
dev.off()


# Up Bar
upDF <- df[which(df$Datatype=='Up'),]
a <- data.frame(Celltype=upDF$Celltype, label='eqtl_spec', count=upDF$eqtl_spec)
b <- data.frame(Celltype=upDF$Celltype, label='common', count=upDF$common)
finalDF <- rbind(a, b)
finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$label <- factor(finalDF$label, levels=c('eqtl_spec','common'))

p1 <- ggplot() + 
    geom_bar(data=finalDF, aes(x=Celltype, y=count, fill=label), stat='identity') +
    scale_fill_manual(values=c('darkgrey', 'darkorange')) +
    labs(x='', y='count') +
    coord_flip() + 
    theme_classic() +
    labs(x='',y='Number of eQTLs', title='Up cRE')

pdf('Plots/Up.bar.pdf')
plot(p1)
dev.off()



# MergedDys Bar
mergedDF <- df[which(df$Datatype=='MergedDysPeak'),]
a <- data.frame(Celltype=mergedDF$Celltype, label='eqtl_spec', count=mergedDF$eqtl_spec)
b <- data.frame(Celltype=mergedDF$Celltype, label='common', count=mergedDF$common)
finalDF <- rbind(a, b)
finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$label <- factor(finalDF$label, levels=c('eqtl_spec','common'))

p1 <- ggplot() + 
    geom_bar(data=finalDF, aes(x=Celltype, y=count, fill=label), stat='identity') +
    scale_fill_manual(values=c('darkgrey', 'darkorange')) +
    labs(x='', y='count') +
    coord_flip() + 
    theme_classic() +
    labs(x='',y='Number of eQTLs', title='Merged dys cRE')

pdf('Plots/MergedDysPeak.bar.pdf')
plot(p1)
dev.off()




