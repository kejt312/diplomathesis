library(rrcov)
library(dplyr)
library(ggplot2)
library(vegan)
library(compositions)
library(factoextra)
library(cluster)
library(ape)
library(ecole)

setwd("C:/Users/katar/Desktop/project/mybioma-dinghy/statistics/CRC-exploration")

colors_disease <- c("aquamarine2", "royalblue1", "deeppink3")
names(colors_disease) <- c("healthy", "adenoma", "CRC")

colors_studies <- c("seagreen3", "darkorange", "yellow2", "blueviolet")
names(colors_studies) <- c("Baxter", "Zeller", "Zackular", "Yang")

# Importing datasets from python
features <- read.delim("features_reduced_python.tsv", header=TRUE,  row.names=1)
metadata <- read.delim("metadata_python.tsv", header=TRUE)

# Removing zeros
features_raw <- features
pseudocount <- 0.000001
features <- features[rowSums(features) > 0,] # remove rows containing only 0
features <- apply(features, MARGIN = 2, function(x) {x/sum(x)}) #relative counts
features[features == 0] <- pseudocount
features <- as.data.frame(features)

# Centered log ratio transformation
features <- as.data.frame(clr(features))

# Robust PCA
features_rob_PCA <- PcaHubert(t(features), mcd = FALSE, scale = TRUE) # alpha: fraction of outliers the algorithm should resist, default: alpha=0.75, set alpha higher 
#features_rob_PCA_100 <- PcaHubert(t(features), kmax = 100, k = 100, mcd = FALSE, scale = TRUE) # alpha: fraction of outliers the algorithm should resist, default: alpha=0.75, set alpha higher 

plot(features_rob_PCA) # multivariete outliers

PCs <- as.data.frame(features_rob_PCA@scores)

screeplot(features_rob_PCA, type="lines", main="Scree plot: 10 computed PCs")
#screeplot(features_rob_PCA_100, type="lines", main="Scree plot: 100 computed PCs", k = 20)


# Variance in first 10 PCs: (cumsum(features_rob_PCA@eig0)/features_rob_PCA@totvar0)[1:10]
# % variance in first 2 PCs
pc1_var <- (cumsum(features_rob_PCA@eig0)/features_rob_PCA@totvar0)[1]*100
pc2_var <- ((cumsum(features_rob_PCA@eig0)/features_rob_PCA@totvar0)[2]-
              (cumsum(features_rob_PCA@eig0)/features_rob_PCA@totvar0)[1])*100
summary(features_rob_PCA_100)

## Visualisation

qplot(x=PC1, y=PC2, data=PCs, colour=factor(metadata$diagnosis)) +
  stat_ellipse(geom = "polygon", aes(fill = factor(metadata$diagnosis), color = factor(metadata$diagnosis)),
               alpha = 0.1) +
  scale_fill_manual(values = colors_disease, guide = "none") +
  scale_color_manual(values = colors_disease) +
  labs(x= paste("PC 1 (", round(pc1_var, 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(pc2_var, 2), "%)", sep = ""),
       color = "Group") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, colour = "black")) 

qplot(x=PC1, y=PC2, data=PCs, colour=factor(metadata$study)) +
  stat_ellipse(geom = "polygon", aes(fill = factor(metadata$study), color = factor(metadata$study)),
               alpha = 0.1) +
  scale_fill_manual(values = colors_studies, guide = "none") +
  scale_color_manual(values = colors_studies) +
  labs(x= paste("PC 1 (", round(pc1_var, 2), "%)", sep = ""), 
       y = paste("PC 2 (", round(pc2_var, 2), "%)", sep = ""),
       color = "Study") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, colour = "black"))

# PCoA

#dist.bc <- dsvdis(t(features),'bray/curtis')
#pcoa <- pco(dis.bc,k=10)
dist_aitchison <- vegdist(t(features_raw), method = "aitchison", pseudocount=0.000001)
pcoa_aitchison <- pcoa(dist_aitchison)
plot(pcoa_aitchison)

biplot.pcoa(pcoa_aitchison) +
  scale_fill_manual(values = colors_studies, guide = "none")

biplot(pcoa_aitchison, Y=NULL, plot.axes = c(1,2))

qplot(pcoa_aitchison$vectors[,"Axis.1"], pcoa_aitchison$vectors[,"Axis.2"], colour=factor(metadata$diagnosis)) +
  stat_ellipse(geom = "polygon", aes(fill = factor(metadata$diagnosis), color = factor(metadata$diagnosis)),
               alpha = 0.1) +
  scale_fill_manual(values = colors_disease, guide = "none") +
  scale_color_manual(values = colors_disease) +
  labs(x= "PCoA 1", y = "PCoA 2", color = "Group") +
  #guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, colour = "black"))

qplot(pcoa_aitchison$vectors[,"Axis.1"], pcoa_aitchison$vectors[,"Axis.2"], colour=factor(metadata$study)) +
  stat_ellipse(geom = "polygon", aes(fill = factor(metadata$study), color = factor(metadata$study)),
               alpha = 0.1) +
  scale_fill_manual(values = colors_studies, guide = "none") +
  scale_color_manual(values = colors_studies) +
  labs(x= "PCoA 1", y = "PCoA 2", color = "Study") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, colour = "black"))

permanova_diagnosis <- adonis2(dist_aitchison ~ diagnosis, data = metadata)
permanova_study <- adonis2(dist_aitchison ~ study, data = metadata)
permanova_pw_diagnosis <- permanova_pairwise(dist_aitchison, grp = metadata$diagnosis)
permanova_pw_study <- permanova_pairwise(dist_aitchison, grp = metadata$study)


# Comparison of raw and clr transformed data
histogram(features[, 1], xlab = paste("clr-transformed ASV counts in sample", names(features)[1], sep = " "), 
          type="count", col="blue", breaks=100)
histogram(features_raw[, 1], xlab = paste("ASV counts in sample", names(features_raw)[1], sep = " "), 
          type="count", col="blue", breaks=100)

## Outliers
outliers <- data.frame(No = which(features_rob_PCA@flag == FALSE))
metadata$No <- seq.int(nrow(metadata))

outliers_info <- inner_join(outliers, metadata) # outlier samples and metadata information
metadata$No <- NULL

# Study
# length(which(outliers_info$study == "Yang")) # [1] 289
# length(which(outliers_info$study == "Zeller")) # [1] 127
# length(which(outliers_info$study == "Zackular")) # [1] 36
# length(which(outliers_info$study == "Baxter")) # [1] 65

# Disease status
# length(which(outliers_info$disease_stat == "healthy")) # [1] 287
# length(which(outliers_info$disease_stat == "CRC")) # [1] 230

# PCA without outliers
selected_features <- features[, !(names(features) %in% rownames(outliers))]
selected_metadata <- metadata[!(rownames(metadata) %in% rownames(outliers)), ]

selected_features <- selected_features[rowSums(selected_features) > 0,] # remove rows containing only 0
selected_features <- apply(selected_features, MARGIN = 2, function(x) {x/sum(x)}) #relative counts
selected_features[selected_features == 0] <- pseudocount
selected_features <- as.data.frame(selected_features)

selected_features <- as.data.frame(clr(selected_features))

selected_features_rob_PCA <- PcaHubert(t(selected_features), mcd = FALSE, scale = TRUE)

PCs_2 <- as.data.frame(selected_features_rob_PCA@scores)

qplot(x=PC1, y=PC2, data=PCs_2, colour=factor(selected_metadata$study)) +
  theme(legend.position="bottom") 

length(which(selected_features_rob_PCA@flag == FALSE)) #[1] 317

screeplot(selected_features_rob_PCA, type="lines", main="Screeplot: robust PCA")

# PCA without Yang
yang <- rownames(metadata)[which(metadata$study == "Yang")]

features_no_yang <- features[, !(names(features) %in% yang)]
metadata_no_yang <- metadata[!(rownames(metadata) %in% yang), ]

features_no_yang <- features_no_yang[rowSums(features_no_yang) > 0,] # remove rows containing only 0
features_no_yang <- apply(features_no_yang, MARGIN = 2, function(x) {x/sum(x)}) #relative counts
features_no_yang[features_no_yang == 0] <- pseudocount
features_no_yang <- as.data.frame(features_no_yang)

features_no_yang <- as.data.frame(clr(features_no_yang))

features_n_y_rob_PCA <- PcaHubert(t(features_no_yang), mcd = FALSE, scale = TRUE)

PCs_3 <- as.data.frame(features_n_y_rob_PCA@scores)

screeplot(features_n_y_rob_PCA, type="lines", main="Screeplot: robust PCA")

qplot(x=PC1, y=PC2, data=PCs_3, colour=factor(metadata_no_yang$study)) +
  theme(legend.position="bottom") 

length(which(features_n_y_rob_PCA@flag == FALSE)) #[1] 244

# Clustering

## Elbow plot
wssplot <- function(data, nc=15, seed=123){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of groups",
       ylab="Sum of squares within a group")}

wssplot(t(features), nc=20)

## 2 clusters 
kmeans.2 <- kmeans(t(features), centers = 2) # Within cluster sum of squares by cluster: 8.8%
cl1 <- data.frame(No = which(kmeans.2$cluster == 1))
cl2 <- data.frame(No = which(kmeans.2$cluster == 2))
metadata$No <- seq.int(nrow(metadata))

cl1_info <- inner_join(cl1, metadata)
cl2_info <- inner_join(cl2, metadata)

# table(cl1_info$study) # Baxter     Yang Zackular   Zeller 
#                             81      957       16       13 

# table(cl2_info$study) # Baxter     Yang Zackular   Zeller 
#                            456       75       72      116

## 3 clusters
kmeans.3 <- kmeans(t(features), centers = 3)
cl1.3 <- data.frame(No = which(kmeans.3$cluster == 1))
cl2.3 <- data.frame(No = which(kmeans.3$cluster == 2))
cl3.3 <- data.frame(No = which(kmeans.3$cluster == 3))

cl1.3_info <- inner_join(cl1.3, metadata)
cl2.3_info <- inner_join(cl2.3, metadata)
cl3.3_info <- inner_join(cl3.3, metadata)

# table(cl1.3_info$study) # Baxter     Yang Zackular   Zeller 
#                              468        1       73      118 

# table(cl2.3_info$study) # Yang 
#                            568 

# table(cl3.3_info$study) # Baxter     Yang Zackular   Zeller 
#                               69      463       15       11

## Clustering validation
sil <- silhouette(kmeans.3$cluster, dist(t(features)))
summary(sil)

## Visualizing clusters in PC-plots
cl <- data.frame(cluster = kmeans.2$cluster)
cl$Run <- rownames(cl)
metadata <- inner_join(cl, metadata)

qplot(x=PC1, y=PC2, data=PCs, colour=factor(metadata$cluster)) +
  theme(legend.position="bottom")




