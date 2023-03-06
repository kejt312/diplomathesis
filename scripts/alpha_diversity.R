library(vegan)
library(ggplot2)
library(BiodiversityR)
library(FSA)
library(ggpubr)
library(rstatix)
library(dplyr)

#Importing ASV counts and metadata

featuresDf <- read.csv(
  file = "C:/Users/katar/Desktop/project/phyloseq/ASV_raw_counts_merged.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = TRUE,
  skip = 1,
  comment.char = "",
  check.names = FALSE
)

metadata <- read.delim("metadata_python.tsv", header=TRUE)

# Adjusting data in ASV data frame

## Differences between the ASV and metadata data frames
differences <- setdiff(colnames(featuresDf), unlist(metadata[,'Run']))

## Removing samples 
featuresDf <- featuresDf[, !(colnames(featuresDf) %in% differences)]
samples <- colnames(featuresDf)
#rm(featuresDf)

# Transpose
featuresDf <- as.data.frame(t(featuresDf))

## takes long to calculate so saved as RData
#alpha_div_ASV <- diversity(featuresDf, index = "shannon") # takes long to calculate so saved as RData 
#evenness_ASV <- diversityresult(featuresDf, index = "Jevenness", method = "each site")

load("alpha_div_ASV.RData")
load("evenness_ASV.RData")
load("samples_alpha_div_ASV.RData")

ASV_alpha <- as.data.frame(cbind(samples, alpha_div_ASV, disease_stat = metadata$disease_stat, diagnosis = metadata$diagnosis))
ASV_alpha$alpha_div_ASV <- as.numeric(ASV_alpha$alpha_div_ASV)

shapiro.test(ASV_alpha[which(ASV_alpha$diagnosis=="healthy"), "alpha_div_ASV"]) 
shapiro.test(ASV_alpha[which(ASV_alpha$diagnosis=="CRC"), "alpha_div_ASV"])
shapiro.test(ASV_alpha[which(ASV_alpha$diagnosis=="adenoma"), "alpha_div_ASV"])
#all not normally distributed 

kruskal.test(x=ASV_alpha$alpha_div, g=ASV_alpha$diagnosis) #p-value = 2.227e-11 < 0.05 -> significant

# p-values corrected using Benjamini-Hochberg FDR method for multiple hypothesis testing
dunnTest(x=ASV_alpha$alpha_div, g=ASV_alpha$diagnosis, method = "bh") #all significant

res.kruskal <- ASV_alpha %>% kruskal_test(alpha_div_ASV ~ diagnosis)

# Pairwise comparisons
pwc <- ASV_alpha %>% 
  dunn_test(alpha_div_ASV ~ diagnosis, p.adjust.method = "BH") 
pwc <- pwc %>% add_xy_position(x = "diagnosis")
pwc <- pwc %>% mutate(y.position = c(6.7, 7.2, 7.7))
pwc$p.adj <-  signif(pwc$p.adj,3)
pwc$xmin <- NULL
pwc$xmax <- NULL

ggboxplot(ASV_alpha, x = "diagnosis", y = "alpha_div_ASV", xlab = "Disease status", ylab = "Shannon diversity index", fill = "diagnosis", 
          palette = c("aquamarine2", "royalblue1", "deeppink3")) +
  theme(legend.position="none") + 
  stat_pvalue_manual(pwc, hide.ns = FALSE, bracket.size = 0.5, label = "p = {p.adj}") +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )  

# Evenness
ASV_evenness <- as.data.frame(cbind(samples, evenness_ASV, disease_stat = metadata$disease_stat, diagnosis = metadata$diagnosis))

shapiro.test(features_evenness[which(ASV_evenness$diagnosis=="healthy"), "Jevenness"]) 
shapiro.test(features_evenness[which(ASV_evenness$diagnosis=="CRC"), "Jevenness"]) 
shapiro.test(features_evenness[which(ASV_evenness$diagnosis=="adenoma"), "Jevenness"]) 

kruskal.test(x=ASV_evenness$Jevenness, g=ASV_evenness$diagnosis) #p-value < 2.2e-16
dunnTest(x=ASV_evenness$Jevenness, g=ASV_evenness$diagnosis, method = "bh") #all significant

res.kruskal <- ASV_evenness %>% kruskal_test(Jevenness ~ diagnosis)

# Pairwise comparisons
pwc <- ASV_evenness %>% 
  dunn_test(Jevenness ~ diagnosis, p.adjust.method = "BH") 
pwc <- pwc %>% add_xy_position(x = "diagnosis")
pwc <- pwc %>% mutate(y.position = c(1.1, 1.2, 1.3))
pwc$p.adj <-  signif(pwc$p.adj,3)
pwc$xmin <- NULL
pwc$xmax <- NULL

ggboxplot(ASV_evenness, x = "diagnosis", y = "Jevenness", xlab = "Disease status", ylab = "Shannon evenness index", fill = "diagnosis", 
          palette = c("aquamarine2", "royalblue1", "deeppink3")) +
  theme(legend.position="none") + 
  stat_pvalue_manual(pwc, hide.ns = FALSE, bracket.size = 0.5, label = "p = {p.adj}") +
  labs(
    subtitle = get_test_label(res.kruskal, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )  
