## MaAsLin2

# Differentially abundant species between CRC and healthy 

maaslin_CRC_healthy <- read.delim("maaslin_CRC_healthy.tsv", header=TRUE)
maaslin2_codes <- read.delim("maaslin2_codes.tsv")

names(maaslin_CRC_healthy)[names(maaslin_CRC_healthy) == 'feature'] <- 'code'

CRC_healthy <- merge(maaslin_CRC_healthy, maaslin2_codes, by="code", all.y = FALSE)
CRC_healthy_significant <- CRC_healthy[CRC_healthy$qval<0.05,]
CRC_healthy_significant <- CRC_healthy_significant[order(CRC_healthy_significant$coef, decreasing = TRUE), ]

# Considering exact effect size
#maaslin_CRC_CRC_significant <- CRC_healthy_significant[CRC_healthy_significant$coef<(-1.5), "species"] #associated with CRC (8)
#maaslin_CRC_healthy_significant <-  CRC_healthy_significant[CRC_healthy_significant$coef>1.5, "species"] #associated with healthy (36)

# Considering only +/- of the effect size
maaslin_CRC_CRC_significant <- CRC_healthy_significant[CRC_healthy_significant$coef<0, "species"] 
maaslin_CRC_healthy_significant <-  CRC_healthy_significant[CRC_healthy_significant$coef>0, "species"] 


# Differentially abundant species between adenoma and healthy 

maaslin_adenoma_healthy <- read.delim("maaslin_adenoma_healthy.tsv", header=TRUE)
names(maaslin_adenoma_healthy)[names(maaslin_adenoma_healthy) == 'feature'] <- 'code'

adenoma_healthy <- merge(maaslin_adenoma_healthy, maaslin2_codes, by="code", all.y = FALSE)

maaslin_ADM_adenoma_significant <- adenoma_healthy[adenoma_healthy$qval<0.05, "species"] #just one species, coef = -1.795446 -> associated with adenoma

#rm(adenoma_healthy, maaslin_adenoma_healthy, CRC_healthy, CRC_healthy_significant, maaslin_CRC_healthy, maaslin2_codes)

#write.csv(maaslin_CRC_CRC_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/maaslin_CRC_CRC_significant")
#write.csv(maaslin_CRC_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/maaslin_CRC_healthy_significant")
#write.csv(maaslin_ADM_adenoma_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/maaslin_ADM_adenoma_significant")

## ALDEx2

load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/aldex_ADM_HTY.RData")
load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/aldex_CRC_HTY.RData")

aldex_ADM_HTY$species <- rownames(aldex_ADM_HTY)
aldex_CRC_HTY$species <- rownames(aldex_CRC_HTY)

aldexADM_healthy_significant <- aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect>0), "species"]
aldexADM_adenoma_significant <- aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect<0), "species"]
# aldexADM_healthy_significant <- aldex_ADM_HTY[
#   which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect>0.3), "species"]
# aldexADM_adenoma_significant <- aldex_ADM_HTY[
#   which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect<(-0.3)), "species"]

aldexCRC_healthy_significant <- aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect>0), "species"]
aldexCRC_CRC_significant <- aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect<0), "species"]
# aldexCRC_healthy_significant <- aldex_CRC_HTY[
#   which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect>0.1), "species"]
# aldexCRC_CRC_significant <- aldex_CRC_HTY[
#   which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect<(-0.1)), "species"]

#rm(aldex_ADM_HTY, aldex_CRC_HTY)

#write.csv(aldexCRC_CRC_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/aldexCRC_CRC_significant")
#write.csv(aldexCRC_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/aldexCRC_healthy_significant")
#write.csv(aldexADM_adenoma_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/aldexADM_adenoma_significant")
#write.csv(aldexADM_healthy_significant, file = "C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/DA/aldexADM_healthy_significant")

## ANCOM

# load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMCRC_CRC_significant.RData")
# load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMCRC_healthy_significant.RData")
# load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMADM_ADM_significant.RData")
# load("C:/Users/katar/Desktop/project/phyloseq/R_objects_merged_data/ANCOM_outputs/ANCOMADM_healthy_significant.RData")

load("ANCOMBC2_CRC_HTY_sign.RData")
load("ANCOMBC2_ADM_HTY_sign.RData")

ANCOMCRC_healthy_significant <- ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                                                        ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`>0), "taxa"]
ANCOMCRC_CRC_significant <- ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                                                     ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`<0), "taxa"]

ANCOMADM_healthy_significant <- ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 &
                                                     ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`>0), "taxa"]
ANCOMADM_ADM_significant <- ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 &
                                                  ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`<0), "taxa"]


## All significant from all methods
maaslin_CRC <-  c(maaslin_CRC_CRC_significant, maaslin_CRC_healthy_significant)
maaslin_ADM <- c(maaslin_ADM_adenoma_significant)

aldex_CRC <- c(aldexCRC_CRC_significant, aldexCRC_healthy_significant)
aldex_ADM <- c(aldexADM_adenoma_significant, aldexADM_healthy_significant)

ANCOM_CRC <- ANCOMBC2_CRC_HTY$taxa
ANCOM_ADM <- ANCOMBC2_ADM_HTY$taxa

## Venn diagrams

#library(VennDiagram)
#library(dplyr)
#library(RColorBrewer)

# myCol <- brewer.pal(7, "Pastel2")
# 
# venn.diagram(
#   x = list(maaslin_CRC, aldex_CRC, ANCOM_CRC),
#   category.names = c("MaAsLin2", "ALDEx2", "ANCOM BC"),
#   filename = "DA_CRC_comparison.png",
#   #cat.default.pos = "outer",
#   lwd = 2,
#   lty = 'blank',
#   cat.pos = c(325, 35, 0),
#   cat.dist = c(0.055, 0.055, 0.055),
#   rotation = 1,
#   fill = myCol[3:5],
#   fontfamily = "sans",
#   cat.fontfamily = "sans"
# )
# 
# venn.diagram(
#   x = list(maaslin_ADM, aldex_ADM, ANCOM_ADM),
#   category.names = c("MaAsLin2", "ALDEx2", "ANCOM BC"),
#   filename = "DA_ADM_comparison.png",
#   #cat.default.pos = "outer",
#   lwd = 2,
#   lty = 'blank',
#   cat.pos = c(315, 45, 180),
#   cat.dist = c(0.015, 0.025, 0.025),
#   rotation = 1,
#   fill = myCol[3:5],
#   fontfamily = "sans",
#   cat.fontfamily = "sans"
# )

# library(ggVennDiagram)
# 
# ggVennDiagram(list(maaslin_CRC, aldex_CRC, ANCOM_CRC))
# ggVennDiagram(list(maaslin_CRC, aldex_CRC, ANCOM_CRC), color = "black", lwd = 0.8, lty = 1,
#               category.names = c("MaAsLin2", "ALDEx2", "ANCOM BC"), 
#               label = "count")
# ggVennDiagram(list(maaslin_CRC, aldex_CRC, ANCOM_CRC), color = "black", lwd = 0.8, lty = 1) + 
#   scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

library(ggplot2)
library(ggVennDiagram)

x <- list(maaslin_CRC, aldex_CRC, ANCOM_CRC)
x <- list(maaslin_ADM, aldex_ADM, ANCOM_ADM)
venn <- Venn(x)
data <- process_data(venn)

ggplot() +
  # 1. region count layer
  geom_sf(aes(fill = count), data = venn_region(data)) +
  # 2. set edge layer
  geom_sf(colour = "black", data = venn_setedge(data), size = 0.1) +
  # 3. set label layer
  geom_sf_text(aes(label = c("MaAsLin2", "ALDEx2", "ANCOM-BC")), data = venn_setlabel(data), vjust="inward",
               hjust="inward") +
  # 4. region label layer
  geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), 
                data = venn_region(data),
                size = 3) +
  scale_fill_gradient(low = "white", high = "deeppink3") +  # CRC: high = "deeppink3" / adenoma: high = "royalblue1"
  theme_void()

#intersect(intersect(aldexCRC_CRC_significant, ANCOMCRC_CRC_significant), maaslin_CRC_CRC_significant)
#intersect(intersect(aldexCRC_healthy_significant, ANCOMCRC_healthy_significant), maaslin_CRC_healthy_significant)
all_methods_CRCHTY <- intersect(intersect(maaslin_CRC, aldex_CRC), ANCOM_CRC)
write.csv(all_methods_CRCHTY, file="all_methods_CRCHTY.csv")

#intersect(aldexADM_adenoma_significant, ANCOMADM_ADM_significant)
#intersect(aldexADM_healthy_significant, ANCOMADM_healthy_significant)
all_methods_ADMHTY <- intersect(aldex_ADM, ANCOM_ADM)
write.csv(all_methods_ADMHTY, file="all_methods_ADMHTY.csv")

# Saving tables with p-values and effect-size

write.csv(CRC_healthy_significant[CRC_healthy_significant$coef>0,], file="maaslin_CRCvsHTY_HTYsignificant.csv")
write.csv(CRC_healthy_significant[CRC_healthy_significant$coef<0,], file="maaslin_CRCvsHTY_CRCsignificant.csv")
write.csv(adenoma_healthy[adenoma_healthy$qval<0.05,], file="maaslin_ADMvsHTY_ADMsignificant.csv")

write.csv(aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect>0),], 
  file="aldex_CRCvsHTY_HTYsignificant.csv")
write.csv(aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect<0),], 
  file="aldex_CRCvsHTY_ADMsignificant.csv")
write.csv(aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect>0),], 
  file="aldex_ADMvsHTY_HTYsignificant.csv")
write.csv(aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect<0),], 
  file="aldex_ADMvsHTY_ADMsignificant.csv")

write.csv(ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                                   ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`>0),], 
          file="ancombc_CRCvsHTY_HTYsignificant.csv")
write.csv(ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                                   ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`<0),], 
          file="ancombc_CRCvsHTY_CRCsignificant.csv")
write.csv(ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 &
                                   ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`>0),], 
          file="ancombc_ADMvsHTY_HTYsignificantt.csv")
write.csv(ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 &
                                   ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`<0),], 
          file="ancombc_ADMvsHTY_ADMsignificant.csv")

## For thesis, top 10 significant taxa with highest effect-size

library(knitr)

# ALDEx2

x <- aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect>0),]
x <- aldex_CRC_HTY[
  which(aldex_CRC_HTY$we.eBH < 0.05 & aldex_CRC_HTY$wi.eBH < 0.05 & aldex_CRC_HTY$effect<0),]

x <- aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect>0),]
x <- aldex_ADM_HTY[
  which(aldex_ADM_HTY$we.eBH < 0.05 & aldex_ADM_HTY$wi.eBH < 0.05 & aldex_ADM_HTY$effect<0),]


y <- x[order(abs(x$effect), decreasing = TRUE),c("we.eBH", "wi.eBH", "effect")]
y <- y[1:10,]
kable(y, "latex", digits = c(27,27,3))

# ANCOM BC

x <- ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                              ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`>0),]
x <- ANCOMBC2_CRC_HTY[which(ANCOMBC2_CRC_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                              ANCOMBC2_CRC_HTY$`out$res$lfc_diagnosishealthy`<0),]

x <- ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                              ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`>0),]
x <- ANCOMBC2_ADM_HTY[which(ANCOMBC2_ADM_HTY$`out$res$q_diagnosishealthy`<0.05 & 
                              ANCOMBC2_ADM_HTY$`out$res$lfc_diagnosishealthy`<0),]
rownames(x) <- x$taxa
x$taxa <- NULL
names(x) <- c("q", "lfc")

y <- x[order(abs(x$lfc), decreasing = TRUE),]
y <- y[1:10,]
kable(y, "latex", digits = c(37,3))

# Maaslin

x <- CRC_healthy_significant[, c("species", "qval", "coef")]

x <- adenoma_healthy[adenoma_healthy$qval<0.05, c("species", "qval", "coef")]

rownames(x) <- x$species
x$species <- NULL
x <- x[x$coef>0, ]
x <- x[x$coef<0, ]

y <- x[order(abs(x$coef), decreasing = TRUE),]
y <- y[1:10,]
kable(y, "latex", digits = c(47,3))

kable(x, "latex", digits = c(47,3))
