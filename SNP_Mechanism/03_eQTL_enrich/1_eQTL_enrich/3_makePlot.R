#!/anaconda2/scEnv


library(ggplot2)
library(RColorBrewer)



dir.create('Plots')

celltypes <- c('DopaN', 'GabaN', 'Oligo', 'OPC', 'Ast', 'Micro', 'Endo', 'Peri')


## celltype eQTL barplot
data1 <- read.table('sigEQTL.cRE_asso.txt', header=T)

data1$CellType <- factor(data1$CellType, levels=celltypes)

pdf('Plots/celltype_eQTL.bar.pdf')
barplot(rev(data1$nEQTL), names.arg=rev(data1$CellType), horiz=T, las=2, xlab='Number of eQTLs')
dev.off()


## dys cRE enrichment (forest plot)
data1 <- read.table('sigEQTL.dys_cRE_asso.txt', header=T)

downDF <- data1[which(data1$DataType=='Down'), ]
upDF <- data1[which(data1$DataType=='Up'), ]
mergedDF <- data1[which(data1$DataType=='MergedDysPeak'), ]

# Down cRE
finalDF <- NULL
for (i in 1:nrow(downDF)){
    df <- data.frame('dys_cRE' = c(downDF$nEQTL[i], downDF$expected[i]), 'genome' = c(downDF$totalEqtl[i]-downDF$nEQTL[i], downDF$totalEqtl[i]-downDF$expected[i]), row.names = c('eQTL', 'expectation'))
    fisherRes <- fisher.test(df)

    celltype <- downDF$CellType[i]
    datatype <- downDF$DataType[i]
    nEQTL <- downDF$nEQTL[i]
    expected <- downDF$expected[i]
    fisherP <- fisherRes$p.value
    OR <- fisherRes$estimate
    lowConf <- fisherRes$conf.int[1] 
    highConf <- fisherRes$conf.int[2]

    newLine <- c(celltype, datatype, nEQTL, expected, fisherP, OR, lowConf, highConf)
    finalDF <- rbind(finalDF, newLine)
}

finalDF <- data.frame(finalDF)
colnames(finalDF) <- c('Celltype','Datatype','nEQTL','expected','fisherP','OR','lowConf','highConf')
rownames(finalDF) <- NULL

finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$OR <- as.numeric(finalDF$OR)
finalDF$lowConf <- as.numeric(finalDF$lowConf)
finalDF$highConf <- as.numeric(finalDF$highConf)

print(finalDF)

# 
p1 <- ggplot(finalDF, aes(x= OR, y = Celltype, xmin = lowConf, xmax = highConf)) +
    geom_point(shape = 0, size = 3, fill='darkgrey') +  
    geom_errorbarh(height = 0) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    xlim(0, 5) +
    labs(x='Odds Ratio (95% CI)', y='', title='Down cRE /eQTL asso') +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Down_dysCre.eQTL_enrich.forestplot.pdf')
plot(p1)
dev.off()



# Up cRE
finalDF <- NULL
for (i in 1:nrow(upDF)){
    df <- data.frame('dys_cRE' = c(upDF$nEQTL[i], upDF$expected[i]), 'genome' = c(upDF$totalEqtl[i]-upDF$nEQTL[i], upDF$totalEqtl[i]-upDF$expected[i]), row.names = c('eQTL', 'expectation'))
    fisherRes <- fisher.test(df)

    celltype <- upDF$CellType[i]
    datatype <- upDF$DataType[i]
    nEQTL <- upDF$nEQTL[i]
    expected <- upDF$expected[i]
    fisherP <- fisherRes$p.value
    OR <- fisherRes$estimate
    lowConf <- fisherRes$conf.int[1] 
    highConf <- fisherRes$conf.int[2]

    newLine <- c(celltype, datatype, nEQTL, expected, fisherP, OR, lowConf, highConf)
    finalDF <- rbind(finalDF, newLine)
}

finalDF <- data.frame(finalDF)
colnames(finalDF) <- c('Celltype','Datatype','nEQTL','expected','fisherP','OR','lowConf','highConf')
rownames(finalDF) <- NULL

finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$OR <- as.numeric(finalDF$OR)
finalDF$lowConf <- as.numeric(finalDF$lowConf)
finalDF$highConf <- as.numeric(finalDF$highConf)

print(finalDF)

# 
p1 <- ggplot(finalDF, aes(x= OR, y = Celltype, xmin = lowConf, xmax = highConf)) +
    geom_point(shape = 0, size = 3, fill='darkgrey') +  
    geom_errorbarh(height = 0) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    xlim(0, 12) +
    labs(x='Odds Ratio (95% CI)', y='', title='Up cRE /eQTL asso') +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/Up_dysCre.eQTL_enrich.forestplot.pdf')
plot(p1)
dev.off()




# Merged dyspeak
finalDF <- NULL
for (i in 1:nrow(mergedDF)){
    df <- data.frame('dys_cRE' = c(mergedDF$nEQTL[i], mergedDF$expected[i]), 'genome' = c(mergedDF$totalEqtl[i]-mergedDF$nEQTL[i], mergedDF$totalEqtl[i]-mergedDF$expected[i]), row.names = c('eQTL', 'expectation'))
    fisherRes <- fisher.test(df)

    celltype <- mergedDF$CellType[i]
    datatype <- mergedDF$DataType[i]
    nEQTL <- mergedDF$nEQTL[i]
    expected <- mergedDF$expected[i]
    fisherP <- fisherRes$p.value
    OR <- fisherRes$estimate
    lowConf <- fisherRes$conf.int[1] 
    highConf <- fisherRes$conf.int[2]

    newLine <- c(celltype, datatype, nEQTL, expected, fisherP, OR, lowConf, highConf)
    finalDF <- rbind(finalDF, newLine)
}

finalDF <- data.frame(finalDF)
colnames(finalDF) <- c('Celltype','Datatype','nEQTL','expected','fisherP','OR','lowConf','highConf')
rownames(finalDF) <- NULL

finalDF$Celltype <- factor(finalDF$Celltype, levels=rev(celltypes))
finalDF$OR <- as.numeric(finalDF$OR)
finalDF$lowConf <- as.numeric(finalDF$lowConf)
finalDF$highConf <- as.numeric(finalDF$highConf)

print(finalDF)

# 
p1 <- ggplot(finalDF, aes(x= OR, y = Celltype, xmin = lowConf, xmax = highConf)) +
    geom_point(shape = 0, size = 3, fill='darkgrey') +  
    geom_errorbarh(height = 0) +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
    xlim(0, 10) +
    labs(x='Odds Ratio (95% CI)', y='', title='MergedDys_cRE /eQTL asso') +  
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

pdf('Plots/MergedDysCre.eQTL_enrich.forestplot.pdf')
plot(p1)
dev.off()



