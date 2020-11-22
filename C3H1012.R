# remove the DEC samples
# directories
wd <- "D:/PhD/Analysis/10T/"
setwd(wd)
rawdir <- paste0(wd, "/telocal/")
outdir <- paste0(wd, "/pics/")
load("C3Honly.RData")

# libs
library(DESeq2)
library(plotly)
library(gplots)
library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(RColorBrewer)
library(ggrepel)
design <- read.table("meta", header = T)
rownames(design) <- design$file
design <- design[grepl("C3H10T1", design$batch),]
design$sample <- factor(design$sample)
design$replicate <- factor(design$replicate)
design$batch <- factor(design$batch)

meta <- data.frame(sampleName = paste0(design$file, ".cntTable"),
                   fileName = paste0(design$file, ".cntTable"),
                   sample = design$sample,
                   condition = design$batch,
                   replicate = design$replicate)

load('counts.RData')
head(gs)
head(cts)

rownames(cts) <- gs

# colnames(cts) <- sort(as.vector(design$file))

# truncated CTS
cts_t <- cts %>% dplyr::select(contains("C3H10T1"))

# Create DESeq obj
dds <- DESeqDataSetFromMatrix(countData = cts_t, 
                              colData = design, 
                              design = ~sample)
dds$sample <- relevel(dds$sample, ref = "PA")
rownames(dds) <- gs
keep <- rowSums(counts(dds)) > 0
dds <- dds[keep,]
dds <- DESeq(dds)
# transformation
vst <- varianceStabilizingTransformation(dds)

save(dds, file = "splitC3H.RData")
ddist <- dist(t(assay(vst)))
matC <- as.matrix(ddist)
rownames(mat) <- colnames(mat) <- dds$file
heatmap.2(mat, trace = "none", margin = c(15,15))
library(plotly)
library(heatmaply)
heatmapC <- heatmaply(matC, k_row = 3, k_col = 3)

features <- read.table("data/TE_split.tab", header = T)
features_trunc <- data.frame(name = features$ID, feature = features$class)
getResAndCounts <- function(dds, ... ) {
sample <- results(dds, ... )
sample <- tibble::rownames_to_column(data.frame(sample), var = "name")
sample <- sample %>% filter(padj <= 0.05) %>% left_join(., features_trunc, by = "name") 
sample <- mutate(sample, feature = replace_na(as.vector(sample$feature), "ENSMUS"))
sample_pos <- sample %>% filter(log2FoldChange > 0) %>% dplyr::select(name, feature)
sample_neg <- sample %>% filter(log2FoldChange < 0) %>% dplyr::select(name, feature)
counts <- merge(data.frame(table(sample_pos$feature)), 
                data.frame(table(sample_neg$feature)), 
                by = "Var1", all = T)
counts[is.na(counts)] <- 0
colnames(counts) <- c("Feature", "PosFoldChange", "NegFoldChange")
return(list(sample = sample, counts = counts))
}


H33K36M <- getResAndCounts(dds, contrast = c("sample", "H33K36M", "PA"))
H33K36R <- getResAndCounts(dds, contrast = c("sample", "H33K36R", "PA"))
NSD12DKO <- getResAndCounts(dds, contrast = c("sample", "NSD12DKO", "PA"))
SETD2KO <- getResAndCounts(dds, contrast = c("sample", "SETD2KO", "PA"))
TKO <- getResAndCounts(dds, contrast = c("sample", "TKO", "PA"))
updatemenusC <- list(
  list(
    type = "buttons",
    x = -0.2,
    y = 1,
    buttons = list(
    list(
      label = "H33K36M",
      method = "update",
      args = list(
        list(visible = c(T, T, F, F, F, F, F, F, F, F)),
        list(title = paste0("H33K36M DEG (total: ", nrow(H33K36M$sample),")")))),
    list(
      label = "H33K36R",
      method = "update",
      args = list(
        list(visible = c(F, F, T, T, F, F, F, F, F, F)),
        list(title = paste0("H33K36R DEG (total: ", nrow(H33K36R$sample),")")))),
    list(
      label = "NSD12DKO",
      method = "update",
      args = list(
        list(visible = c(F, F, F, F, T, T, F, F, F, F)),
        list(title = paste0("NSD12DKO DEG (total: ", nrow(NSD12DKO$sample),")")))),
    list(
      label = "SETD2KO",
      method = "update",
      args = list(
        list(visible = c(F, F, F, F, F, F, T, T, F, F)),
        list(title = paste0("SETD2KO DEG (total: ", nrow(SETD2KO$sample),")")))),
    list(
      label = "TKO",
      method = "update",
      args = list(
        list(visible = c( F, F, F, F, F, F, F, F, T, T)),
        list(title = paste0("TKO DEG (total: ", nrow(TKO$sample),")"))))
        )
      )
)
buttonplotC <- {plot_ly(H33K36M$counts, x = ~H33K36M$counts$Feature, y = ~H33K36M$counts$PosFoldChange, type = "bar", 
            name = "+FoldChange", marker = list(color = '#ff7f0e'),
            text = ~H33K36M$counts$PosFoldChange, textposition = "outside",
            textfont = list(color = '#ff7f0e'), visible = T) %>%
  add_trace(H33K36M$counts, x = ~H33K36M$counts$Feature, y = ~H33K36M$counts$NegFoldChange, type = "bar", 
            name = "-FoldChange", marker = list(color = '#1f77b4'),
            text = ~H33K36M$counts$NegFoldChange, textposition = "outside",
            textfont = list(color = '#1f77b4'), visible = T) %>%
  add_trace(H33K36R$counts, x = ~H33K36R$counts$Feature, y = ~H33K36R$counts$PosFoldChange, type = "bar", 
            name = "+FoldChange", marker = list(color = '#ff7f0e'),
            text = ~H33K36R$counts$PosFoldChange, textposition = "outside",
            textfont = list(color = '#ff7f0e'), visible = F) %>%
  add_trace(H33K36R$counts, x = ~H33K36R$counts$Feature, y = ~H33K36R$counts$NegFoldChange, type = "bar", 
            name = "-FoldChange", marker = list(color = '#1f77b4'),
            text = ~H33K36R$counts$NegFoldChange, textposition = "outside",
            textfont = list(color = '#1f77b4'), visible = F) %>%
  add_trace(NSD12DKO$counts, x = ~NSD12DKO$counts$Feature, y = ~NSD12DKO$counts$PosFoldChange, type = "bar", 
            name = "+FoldChange", marker = list(color = '#ff7f0e'),
            text = ~NSD12DKO$counts$PosFoldChange, textposition = "outside",
            textfont = list(color = '#ff7f0e'), visible = F) %>%
  add_trace(NSD12DKO$counts, x = ~NSD12DKO$counts$Feature, y = ~NSD12DKO$counts$NegFoldChange, type = "bar", 
            name = "-FoldChange", marker = list(color = '#1f77b4'),
            text = ~NSD12DKO$counts$NegFoldChange, textposition = "outside",
            textfont = list(color = '#1f77b4'), visible = F) %>%
  add_trace(SETD2KO$counts, x = ~SETD2KO$counts$Feature, y = ~SETD2KO$counts$PosFoldChange, type = "bar", 
            name = "+FoldChange", marker = list(color = '#ff7f0e'),
            text = ~SETD2KO$counts$PosFoldChange, textposition = "outside",
            textfont = list(color = '#ff7f0e'), visible = F) %>%
  add_trace(SETD2KO$counts, x = ~SETD2KO$counts$Feature, y = ~SETD2KO$counts$NegFoldChange, type = "bar", 
            name = "-FoldChange", marker = list(color = '#1f77b4'),
            text = ~SETD2KO$counts$NegFoldChange, textposition = "outside",
            textfont = list(color = '#1f77b4'), visible = F) %>%
  add_trace(TKO$counts, x = ~TKO$counts$Feature, y = ~TKO$counts$PosFoldChange, type = "bar", 
            name = "+FoldChange", marker = list(color = '#ff7f0e'),
            text = ~TKO$counts$PosFoldChange, textposition = "outside",
            textfont = list(color = '#ff7f0e'), visible = F) %>%
  add_trace(TKO$counts, x = ~TKO$counts$Feature, y = ~TKO$counts$NegFoldChange, type = "bar", 
            name = "-FoldChange", marker = list(color = '#1f77b4'),
            text = ~TKO$counts$NegFoldChange, textposition = "outside",
            textfont = list(color = '#1f77b4'), visible = F) %>%
  layout(xaxis = list(title = "Genomic Features", tickangle = 45),
         yaxis = list(title = "Differentially expressed Genes"),
         updatemenus = updatemenusC)}
buttonplotC

save(buttonplotC,  updatemenusC, H33K36M, H33K36R, TKO, NSD12DKO, SETD2KO, heatmapC, file = "C3Honly.RData")
groupInter <- intersect(intersect(intersect(H33K36M$sample$name, H33K36R$sample$name), 
                                  intersect(SETD2KO$sample$name, NSD12DKO$sample$name)), 
                        TKO$sample$name)
length(groupInter) # 491
groupUnion <- union(union(union(H33K36M$sample$name, H33K36R$sample$name), 
                                  union(SETD2KO$sample$name, NSD12DKO$sample$name)), 
                        TKO$sample$name)
length(groupUnion) # 29,353
library(ggplot2)
library(viridis)
ggplot(H33K36M$sample, aes(x = H33K36M$sample$log2FoldChange, 
                           fill = H33K36M$sample$feature)) + 
  geom_histogram(alpha = 1, binwidth = 0.5) +
  scale_fill_manual(values = viridis(n = 12)) + 
  labs(x = "Log Fold Change Distribution", y = "Number of DEGs (bin = 0.5)") +
  theme_bw() + 
  facet_wrap(~H33K36M$sample$feature, scales = "free_y")
