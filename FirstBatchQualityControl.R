library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(ggrepel)
library(RColorBrewer)

##Get intensity file:
data <- read.table("../DataFirst/combined_protein.csv", sep = "\t", header = TRUE, check.names=FALSE)
LFQcol <- append(4, grep("MaxLFQ", colnames(data)))
dataLFQ <- data[,LFQcol]
colnames(dataLFQ) <- gsub(" ", ".", colnames(dataLFQ))

#Filter out duplicated proteins
dup <- dataLFQ %>% group_by(Gene) %>% summarize(frequency = n()) %>% arrange(desc(frequency)) %>% filter(frequency > 1)
dataLFQ_unique <- dataLFQ[!(dataLFQ$Gene %in% dup$Gene),] # - 134 rows = 65 genes

##Get annotation file
expano <- read.csv("../DataFirst/experiment_annotation.csv", sep = "\t", row.names = 1)
expano <- expano[order(expano$sample_name),]
a <- sort(colnames(dataLFQ_unique[,-1]))
expano$sample_nameS <- a
expano2 <- expano
#write.csv(expano2, "../DataFirst/annotation_match.csv")

# Reshape data for analysis
dataLFQ_melt <- melt(dataLFQ_unique)
colnames(dataLFQ_melt)[2] <- "sample_nameS"
dataLFQ_all <- left_join(dataLFQ_melt, expano)
dataLFQ_all <- dataLFQ_all[,-4]

#Number of biological and technical replicates expressing each protein per condition
dataLFQ_all2 <- dataLFQ_all %>%
                group_by(Gene, condition) %>%
                mutate(Nbiological = n_distinct(repBio[value > 0]),
                       Ntechnical = n_distinct(replicate[value > 0]))

#Filter genes that are expressed in 2 biological replicates (4 technical replicates) in at least one sample:
filt1 <- dataLFQ_all2 %>% filter(Nbiological >= 2 & Ntechnical >= 4) %>% distinct(Gene)
filt1 <- unique(filt1$Gene)
dataLFQ_filt1 <- dataLFQ_all %>% filter(Gene %in% filt1) #need the same list of genes for all samples to run PCA.
length(unique(dataLFQ_filt1$Gene)) # Now there are 2155 proteins. lost 1204 ptns

#Number of proteins detected in each sample
dataLFQ_all3 <- dataLFQ_all2 %>%
                group_by(sample_nameS) %>%
                mutate(Nptns = n_distinct(Gene[value > 0]))
dataLFQ_all3$sample <- paste0(dataLFQ_all3$condition, "_", dataLFQ_all3$replicate)

dfplt2 <- unique(dataLFQ_all3[, c(4, 9, 10)])
plt2 <- ggplot(dfplt2, aes(y = Nptns, x = sample, fill = condition)) +
        geom_col() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
        scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
ggsave("../DataFirst/Results/NptnsAllSecBatch.png", plt2, width = 15, height = 7)

# Filter samples with low number of proteins
filt2 <- c("1598_8hr_1", "1598_8hr_2", "1597noh202_5", "1597noh202_6")
dataLFQ_filt1$sample <- paste0(dataLFQ_filt1$condition, "_", dataLFQ_filt1$replicate)
dataLFQ_filt <- dataLFQ_filt1 %>% filter(!(sample %in% filt2)) # - 4 technical replicates.

# Update anotation file:
expano2$sample <- paste0(expano2$condition, "_", expano2$replicate)
expano3 <- expano2 %>% filter(!(sample %in% filt2))
#write.csv(expano3, "../DataFirst/annotation_match_filtered.csv")

## imputate values
dataLFQ_imp <- dataLFQ_filt %>%
                group_by(Gene, condition) %>%
                mutate(imp = mean(value[value > 0])) # imputation calculated even in those with just one value among 6 replicates.

md <- unique(dataLFQ_all2[,c(1, 2, 7, 8)]) # Get NBio and Ntech with value > 0
dataLFQ_imp <- left_join(dataLFQ_imp, md, by = c("Gene", "sample_nameS"))
dataLFQ_imp$Intensity <- ifelse((dataLFQ_imp$value == 0 & dataLFQ_imp$Nbiological >= 2 & dataLFQ_imp$Ntechnical >= 4) == TRUE,
                                dataLFQ_imp$imp, dataLFQ_imp$value) # Insert imputation where appropriated.

#Check the distribution before and after imputation
plt4 <- ggplot(dataLFQ_filt, aes(x = log2(value+1), color = condition)) +
        geom_density() + ggtitle("Norm") + xlim(0,25) + xlab("Intensity")
plt5 <- ggplot(dataLFQ_imp, aes(x = log2(Intensity+1), color = condition)) +
        geom_density() + ggtitle("Imputation") + xlim(0,25) + xlab("Intensity")
panel1 <- ggarrange(plt4, plt5)

ggsave("../DataFirst/Results/ImputationFirstBatch.png", panel1, width = 8, height = 3)

#normalisation
dataLFQ_imp$IntensityNorm <- log2(dataLFQ_imp$Intensity + 1)

plt6 <- ggplot(dataLFQ_imp, aes(x = sample, y = IntensityNorm, fill = condition)) +
        geom_boxplot() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("../DataFirst/Results/NormFirstBatch.png", plt6, width = 7, height = 4)

#Caspase8 control:
plt7 <- ggplot(dataLFQ_imp %>% filter(Gene == "CASP8"), aes(x = sample, y = IntensityNorm)) +
        geom_boxplot() + theme(axis.text.x=element_text(angle = 90, vjust = 0.5))
ggsave("../DataFirst/Results/NormFirstBatch.png", plt7, width = 7, height = 2.5)

#Get protein intensity matrix back after filtering
dataLFQ_impNorm <- dataLFQ_imp[,c(1, 2, 12)]
LFQ <- dcast(dataLFQ_impNorm, Gene ~ sample_nameS)
rownames(LFQ) <- LFQ$Gene
LFQ <- LFQ[-1,-1] # deleting gene with no name in the first row and column Gene
#write.csv(LFQ, "../DataFirst/LFQ_filteredImpNorm_AC.csv")

#PCA
pca <- prcomp(t(LFQ), scale = TRUE)
pca_summary <- summary(pca)

# Make the first two PCs into a data frame for plotting with `ggplot2`
pca_df <- data.frame(pca$x[, 1:2]) %>%
  # Turn samples IDs stored as row names into a column
  tibble::rownames_to_column("sample_nameS") %>%
  # Bring only the variables that we want from the metadata into this data frame
  # here we are going to join by `refinebio_accession_code` values
  dplyr::inner_join(
    dplyr::select(expano2, sample_nameS, condition, replicate),
    by = "sample_nameS"
  )

pca_df$replicate <- as.character(pca_df$replicate)

#Determine colors of each sample in the PCA:
col <- c("#336600", "#66cc00", "#ccff99", "#99cc99", "#333399", "#6666ff", "#ccccff", "#9933ff", "#cc0033")
names(col) <- c("1597_12hr", "1597_8hr", "1597_UT", "1597noh202", "1598_12hr", "1598_8hr", "1598_UT", "1598noh202", "UT_WT")
colScale <- scale_color_manual(name = "condition", values = col)

#Determine shapes of each replicate in the PCA:
shape <- c(15, 0, 17, 24, 18, 5)
names(shape) <- seq(1, 6, 1) 
shapScale <- scale_shape_manual(name = "replicate", values = shape)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, shape = replicate)) +
  geom_point(size = 5) + theme_classic() + colScale + shapScale
ggsave("../DataFirst/Results/PCAFirstBatch.png", pca_plot)
