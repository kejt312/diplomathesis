library(dplyr)
library(ALDEx2)

# Importing features and metadata dataframes
features <- read.delim("features_reduced_python.tsv", header=TRUE,  row.names=1)
metadata <- read.delim("metadata_python.tsv", header=TRUE)
metadata[which(metadata$study=="Zackular"), "Run"] <- gsub("-", "\\.", metadata[which(metadata$study=="Zackular"), "Run"]) #correcting sample names in Zackular study

# Creating conditions 
healthy <- metadata[which(metadata$diagnosis== "healthy"), "Run"]
adenoma <- metadata[which(metadata$diagnosis== "adenoma"), "Run"]
CRC <- metadata[which(metadata$diagnosis== "CRC"), "Run"]

features1 <- features[ , which(colnames(features) %in% healthy | colnames(features) %in% CRC)]
features2 <- features[ , which(colnames(features) %in% healthy | colnames(features) %in% adenoma)]

conditions1 <- colnames(features1)
conditions1[which(conditions1 %in% healthy)] <- "HTY"
conditions1[which(conditions1 %in% CRC)] <- "CRC"

conditions2 <- colnames(features2)
conditions2[which(conditions2 %in% healthy)] <- "HTY"
conditions2[which(conditions2 %in% adenoma)] <- "ADM"

# ALDEx2
aldex_CRC_HTY <- aldex(features1, conditions1, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE) 
aldex_ADM_HTY <- aldex(features2, conditions2, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE)

# Loading aldex results 
load("aldex_ADM_HTY.RData")
load("aldex_CRC_HTY.RData")

# Plots
par(mfrow=c(1,2))
aldex.plot(aldex_ADM_HTY, type="MA", test="welch", cutoff.pval = 0.05, all.col = "black")
#aldex.plot(aldex_ADM_HTY, type="MW", test="welch", cutoff.pval = 0.05, all.col="black")
title(main="Welch's test", cex.main=1)
aldex.plot(aldex_ADM_HTY, type="MA", test="wilcox", cutoff.pval = 0.05, all.col = "black")
#aldex.plot(aldex_ADM_HTY, type="MW", test="wilcox", cutoff.pval = 0.05, all.col="black")
title(main="Wilcoxon rank sum test", cex.main=1)


par(mfrow=c(1,2))
aldex.plot(aldex_CRC_HTY, type="MA", test="welch", cutoff.pval = 0.05, all.col = "black")
#aldex.plot(aldex_CRC_HTY, type="MW", test="welch", cutoff.pval = 0.05, all.col="black", all.cex = 0.25)
title(main="Welch's test", cex.main=1)
aldex.plot(aldex_CRC_HTY, type="MA", test="wilcox", cutoff.pval = 0.05, all.col = "black")
#aldex.plot(aldex_CRC_HTY, type="MW", test="wilcox", cutoff.pval = 0.05, all.col="black")
title(main="Wilcoxon rank sum test", cex.main=1)

# for CRC vs healthy
aldex <- aldex_CRC_HTY
# for adenoma vs healthy
aldex <- aldex_ADM_HTY

found.by.all <- which(aldex$we.eBH < 0.05 & aldex$wi.eBH < 0.05) #diff. abundant species found by both Welch's and Wilcoxon's test
found.by.one <- which(aldex$we.eBH < 0.05 | aldex$wi.eBH < 0.05) #diff. abundant species found by either Welch's or Wilcoxon's test

aldex$No <- seq.int(nrow(aldex))
aldex$species <- rownames(aldex)

aldex.by.all <- aldex[aldex$No %in% found.by.all, "species"]
aldex.by.one <- aldex[aldex$No %in% found.by.one, "species"]

plot_by.all <- as.data.frame(aldex.by.all)
plot_by.all$label <- rownames(plot_by.all)
names(plot_by.all)[names(plot_by.all) == "aldex.by.all"] <- "species"

aldex <- merge(aldex, plot_by.all, by="species", all.x=TRUE)

plot(aldex$diff.win, aldex$diff.btw, pch=19, cex=0.3, col=rgb(0,0,0,0.3),
     xlab="Median Log2 Dispersion", ylab="Median Log2 Difference")
points(aldex$diff.win[found.by.one], aldex$diff.btw[found.by.one], pch=19,
       cex=0.5, col=rgb(0,0,1,0.5))
points(aldex$diff.win[found.by.all], aldex$diff.btw[found.by.all], pch=19,
       cex=0.5, col=rgb(1,0,0,1))
abline(0,1,lty=2)
abline(0,-1,lty=2)
#text(diff.btw ~diff.win, data= aldex, labels=aldex$label, pos=3)
## Features identified by both tests shown in red. Features identified by only one test are shown in blue dots. 
## Non-significant features represent rare features if black and abundant features if grey dots.
## features with Difference < 0 associated with CRC / > 0 healthy