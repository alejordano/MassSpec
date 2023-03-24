library(limma)
library(edgeR)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(enrichplot)
library(RColorBrewer)

## DEG 1598 treated x UT ##
data <- read.csv("../DataFirst/LFQ_filteredImpNorm_AC.csv", row.names=1)

expano <- read.csv("../DataFirst/annotation_match_filtered.csv")
expano <- expano[,-1]

groups <- factor(expano$condition, levels = unique(expano$condition))
design <- model.matrix( ~ 0 + groups )
colnames(design) <- gsub("groups", "", colnames(design))
design

#Limma package:
dge <- DGEList(counts = data, lib.size = colSums(data), 
               norm.factors = rep(1,ncol(data)),
               samples = NULL, group = NULL, genes = NULL,
               remove.zeros = FALSE)
dge <- calcNormFactors(dge)

fold <- max(dge$samples$lib.size) / min(dge$samples$lib.size)
#if the ratio of the largest library size to the smallest is not more than about 3-fold, use limma-trend

#DEG
colnames(design) <- paste0("X", colnames(design))
cont <- makeContrasts(contrasts = c("X1598_8hr-X1598_UT", "X1598_12hr-X1598_UT"), levels=design)
fit <- lmFit(data, design)
fit2b <- contrasts.fit(fit, cont)
fit2b <- eBayes(fit2b, trend=TRUE)
volplt1 <- topTable(fit2b, coef=1, genelist=fit2b$genes, adjust.method="BH", p.value=1, lfc = 0, n = 3000)
volplt2 <- topTable(fit2b, coef=2, genelist=fit2b$genes, adjust.method="BH", p.value=1, lfc = 0, n = 3000)

#write.csv(volplt1, "../DataFirst/Results/DEG_1598_8hrxUT")
#write.csv(volplt2, "../DataFirst/Results/DEG_1598_12hrxUT")

#Determine thresholds. Corrected p-value <= 0.01. logFC?
thr1 <- volplt1 %>% filter(adj.P.Val <= 0.01)
thr2 <- volplt2 %>% filter(adj.P.Val <= 0.01)

p <- ggplot(thr1, aes(x = logFC)) +
  geom_histogram(bins=30,breaks=seq(-12, 20, 1), fill = "blue") +
  scale_x_continuous(breaks = seq(-15, 20, 5), lim = c(-12, 20)) +
  theme_bw() + ggtitle("1598_8hr - 1598_UT\np-value <= 0.01")

p1 <- ggplot(thr2, aes(x = logFC)) +
  geom_histogram(bins=30,breaks=seq(-12, 20, 1), fill = "blue") +
  scale_x_continuous(breaks = seq(-15, 20, 5), lim = c(-12, 20)) +
  theme_bw() + ggtitle("1598_12hr - 1598_UT\np-value <= 0.01")

panel1 <- ggarrange(p, p1)
ggsave("../DataFirst/Results/ThresholtPltsUT.png", panel1, width = 7, height = 2.5)

#volcano plots
plt8 <- ggplot(volplt1, aes(x = logFC, y = -log10(adj.P.Val), col = (adj.P.Val <= 0.01 & logFC >= 5))) + geom_point() +
  scale_color_manual(values = c('lightblue', 'blueviolet')) + ggtitle("1598_8hr-1598_UT") + theme_classic() + theme(legend.position = "none") +
  geom_label_repel(aes(label = ifelse(-log10(adj.P.Val)>30,as.character(rownames(volplt1)),'')), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 100)

plt9 <- ggplot(volplt2, aes(x = logFC, y = -log10(adj.P.Val), col = (adj.P.Val <= 0.01 & logFC >= 5))) + geom_point() +
  scale_color_manual(values = c('lightblue', 'blueviolet')) + ggtitle("1598_12hr-1598_UT") + theme_classic() + theme(legend.position = "none") + 
  geom_label_repel(aes(label = ifelse(-log10(adj.P.Val)>30,as.character(rownames(volplt1)),'')), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 100)

panel2 <- ggarrange(plt8, plt9, ncol = 2, nrow = 1)
ggsave("../DataFirst/Results/VolcanoPltsUT.png", panel2, width = 7, height = 3.5)

### Check DE method:
topg <- head(row.names(volplt1), n = 5)
exp <- dataLFQ_impNorm %>% filter(Gene %in% topg)

exp2 <- merge(exp, expano3, by = "sample_nameS")
exp2 <- exp2[,-4]
exp2$tp <- sapply(strsplit(exp2$condition, split = "_", fixed = TRUE), function(x) (x[2]))
exp2 <- exp2 %>% filter(condition %in% c("1598_8hr", "1598_UT"))

plt <- ggplot(exp2, aes(x = Gene, y = IntensityNorm, fill = tp)) + geom_boxplot()

## GOEA ##
namef <- c("1598_8hr_1598_UT", "1598_12hr_1598_UT")
inputs <- list(volplt1, volplt2)

for (i in seq_along(inputs)) {
  df <- as.data.frame(inputs[[i]]) %>% filter(logFC > 5 & adj.P.Val <= 0.01)
  list <- sort(rownames(df))
  
  ENRUp <- enrichGO(list, ont ="BP", keyType = "SYMBOL",
                    minGSSize = 3, maxGSSize = 800,
                    pvalueCutoff = 0.05, OrgDb = organism,
                    pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  ENRUp2 <- simplify(ENRUp, cutoff = 0.7, by = "p.adjust",
                     select_fun = min, measure = "Wang", semData = NULL)
  
  df2 <- ENRUp2@result
  #write.csv(df2, paste0("../DataFirst/Results/GOEA_", namef[i]))
  
  ggsave(paste0("../DataFirst/Results/GOEAPlt1_", namef[i], ".png"), dotplot(ENRUp2, showCategory=15))
  
  ENRUp3 <- pairwise_termsim(ENRUp2)
  ggsave(paste0("../DataFirst/Results/GOEAPlt2_", namef[i], ".png"), emapplot(ENRUp3), width = 11, height = 11)
  ggsave(paste0("../DataFirst/Results/GOEAPlt3_", namef[i], ".png"),treeplot(ENRUp3), width = 11, height = 10)
}