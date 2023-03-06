library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# Importing features and metadata dataframes
features <- read.delim("features_reduced_python.tsv", header=TRUE,  row.names=1)
metadata <- read.delim("metadata_python.tsv", header=TRUE)
metadata[which(metadata$study=="Zackular"), "Run"] <- gsub("-", "\\.", metadata[which(metadata$study=="Zackular"), "Run"]) #correcting sample names in Zackular study

# Separating healthy+CRC and healthy+adenoma datasets
metadata1 <- metadata[which(metadata$diagnosis == "healthy" | metadata$diagnosis == "CRC"), ]
features1 <- features[, which(colnames(features) %in% metadata1$Run)]

metadata2 <- metadata[which(metadata$diagnosis == "healthy" | metadata$diagnosis == "adenoma"), ]
features2 <- features[, which(colnames(features) %in% metadata2$Run)]

# Creating phyloseq object (with features1/2, metadata1/2)
OTU = otu_table(features1, taxa_are_rows = TRUE)
rownames(metadata1) <- metadata1$Run
META = sample_data(metadata1)
physeq = phyloseq(OTU, META)
rm(OTU, META)

# ANCOMBC
# out = ancombc(phyloseq = physeq, formula = "diagnosis+study",
#               p_adj_method = "BH", lib_cut = 1000,
#               group = "diagnosis", struc_zero = TRUE, neg_lb = FALSE,
#               tol = 1e-5, max_iter = 100, conserve = TRUE,
#               alpha = 0.05, global = TRUE)

out = ancombc2(data=physeq, fix_formula = "diagnosis+study", p_adj_method = "BH", group = "diagnosis")

#out$res$q_val #adjusted p-Values
#out$res$diff_abn #diff abundant taxa

qs <- as.data.frame(out$res$q_diagnosishealthy)
#qs <- as.data.frame(out$res$q_val)
#qs$taxa <- rownames(qs)
qs$taxa <- rownames(out$feature_table)
#significant <- qs[qs$diagnosishealthy<0.05, ]
significant <- qs[qs[, 1]<0.05, ]

#coefs <- as.data.frame(out$res$beta)
#coefs$taxa <- rownames(coefs)

#coefs <- as.data.frame(out$res$lfc)
#coefs$taxa <- rownames(coefs)
coefs <- as.data.frame(out$res$lfc_diagnosishealthy)
coefs$taxa <- rownames(out$feature_table)


# CRC vs healthy
#effect_CRC <- coefs[coefs$diagnosishealthy < (-0.3), ] #-0.3
#effect_healthy <- coefs[coefs$diagnosishealthy > 0.7, ] #0.7
effect_CRC <- coefs[coefs[, 1] < (-0.25), ] 
effect_healthy <- coefs[coefs[, 1] > 0.7, ] 

ANCOMCRC_CRC_significant <- merge(significant, effect_CRC, by="taxa", all = FALSE)
# ANCOMCRC_CRC_significant$diagnosishealthy.x <- NULL
# ANCOMCRC_CRC_significant$diagnosishealthy.y <- NULL

ANCOMCRC_healthy_significant <- merge(significant, effect_healthy, by="taxa", all = FALSE)
# ANCOMCRC_healthy_significant$diagnosishealthy.x <- NULL
# ANCOMCRC_healthy_significant$diagnosishealthy.y <- NULL

# write.csv(ANCOMCRC_CRC_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/ANCOMCRC_CRC_significant")
# write.csv(ANCOMCRC_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/ANCOMCRC_healthy_significant")
# save(ANCOMCRC_CRC_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMCRC_CRC_significant.RData")
# save(ANCOMCRC_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMCRC_healthy_significant.RData")
ANCOMBC2_CRC_HTY <- merge(significant, coefs, by="taxa", all = FALSE)
save(ANCOMBC2_CRC_HTY, file="ANCOMBC2_CRC_HTY.RData")

# Adenoma vs healthy
rm(physeq, out, features1, metadata1)

OTU = otu_table(features2, taxa_are_rows = TRUE)
rownames(metadata2) <- metadata2$Run
META = sample_data(metadata2)
physeq = phyloseq(OTU, META)
rm(OTU, META)

# out = ancombc(phyloseq = physeq, formula = "diagnosis+study",
#               p_adj_method = "BH", lib_cut = 1000,
#               group = "diagnosis", struc_zero = TRUE, neg_lb = FALSE,
#               tol = 1e-5, max_iter = 100, conserve = TRUE,
#               alpha = 0.05, global = TRUE)

out = ancombc2(data=physeq, fix_formula = "diagnosis", p_adj_method = "BH", group = "diagnosis")

qs <- as.data.frame(out$res$q_diagnosishealthy)
qs$taxa <- rownames(out$feature_table)
#qs$taxa <- rownames(qs)
#significant <- qs[qs$diagnosishealthy<0.05, ]
significant <- qs[qs[, 1]<0.05, ]

#coefs <- as.data.frame(out$res$lfc)
#coefs$taxa <- rownames(coefs)
coefs <- as.data.frame(out$res$lfc_diagnosishealthy)
coefs$taxa <- rownames(out$feature_table)

#effect_adenoma <- coefs[coefs$diagnosishealthy < (-0.2), ] #-1.2
#effect_healthy <- coefs[coefs$diagnosishealthy > 1, ] #1
effect_adenoma <- coefs[coefs[, 1] < (-1.2), ] #-1.2
effect_healthy <- coefs[coefs[, 1] > 1, ] #1

#ANCOMADM_ADM_significant <- merge(significant, coefs, by="taxa", all = FALSE)
ANCOMADM_ADM_significant <- merge(significant, effect_adenoma, by="taxa", all = FALSE)
#ANCOMADM_ADM_significant <- merge(significant, coefs, by="taxa", all = FALSE)
#ANCOMADM_ADM_significant$diagnosishealthy.x <- NULL
#ANCOMADM_ADM_significant$diagnosishealthy.y <- NULL

ANCOMADM_healthy_significant <- merge(significant, effect_healthy, by="taxa", all = FALSE)
#ANCOMADM_healthy_significant <- merge(significant, coefs, by="taxa", all = FALSE)
#ANCOMADM_healthy_significant$diagnosishealthy.x <- NULL
#ANCOMADM_healthy_significant$diagnosishealthy.y <- NULL

# write.csv(ANCOMADM_ADM_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/ANCOMADM_ADM_significant")
# write.csv(ANCOMADM_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/ANCOMADM_healthy_significant")
# save(ANCOMADM_ADM_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMADM_ADM_significant.RData")
# save(ANCOMADM_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMADM_healthy_significant.RData")
ANCOMBC2_ADM_HTY <- merge(significant, coefs, by="taxa", all = FALSE)
save(ANCOMBC2_ADM_HTY, file="ANCOMBC2_ADM_HTY.RData")

# Draw volcano plot

#dat_fig <- as.data.frame(cbind(taxa = rownames(out$res$q_val), q=out$res$q_val$diagnosishealthy, coef = out$res$beta$diagnosishealthy))
dat_fig <- as.data.frame(cbind.data.frame(taxa = qs$taxa, q=qs[, 1], coef = coefs[, 1]))
#dat_fig$q <- as.double(dat_fig$q)
dat_fig$y <- -log10(dat_fig$q)
#dat_fig$coef <- as.double(dat_fig$coef)

dat_fig <- separate(data = dat_fig, col = taxa, into = c("d", "p", "c", "o", "f", "g", "s"), sep = ";")
dat_fig$g <- gsub("g__", "", dat_fig$g)
dat_fig$s <- gsub("s__", "", dat_fig$s)
dat_fig$s <- gsub("nan", "", dat_fig$s)
dat_fig$g <- gsub("nan", "", dat_fig$g)
dat_fig$gs <- paste(dat_fig$g, dat_fig$s)

# Color for q-value and effect size 
# dat_fig$myCol <- ifelse(dat_fig$coef < -1.2 & dat_fig$y > (-log10(0.05)) | dat_fig$coef > 1 & 
#                           dat_fig$y > (-log10(0.05)), "red", "black")
# -0.3/-0.25 & 0.7 for CRC/healthy; -1.2 & 1 for adenoma/healthy 

# Color for q-value only
dat_fig$myCol <- ifelse(dat_fig$y > (-log10(0.05)) , "red", "black")

fig = ggplot(data = dat_fig) +
  geom_point(data = dat_fig, aes(x = coef, y = y, col = myCol)) + 
  scale_colour_identity() +
  labs(x = "Effect size", y = "-log10(q)") +
  #ggtitle("ANCOMBC (without bias correction) - adenoma vs. healthy") +
  ggtitle("ANCOMBC - CRC vs. healthy") +
  theme_bw()
fig
  #geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") +
  #geom_vline(xintercept= -0.3, linetype="dashed", color = "red") +
  #geom_vline(xintercept= 0.7, linetype="dashed", color = "red")  
  # geom_text(aes(label=ifelse(coef>0.7 & y > -log10(0.05),as.character(gs),'')),hjust=0,vjust=0)
  # + geom_label_repel(aes(label = ifelse(coef>0.7 & y > -log10(0.05),as.character(gs),'')),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = 20) +
  # geom_label_repel(aes(label = ifelse(coef<-0.3 & y > -log10(0.05),as.character(gs),'')),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50',
  #                  max.overlaps = 20) +
  # theme_classic()



