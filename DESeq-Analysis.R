library(DESeq2)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(openxlsx)

setwd("C:/Users/User/Downloads/")
metadata <- read.delim("metadata.tsv")
counts <- read.delim("counts.tsv")

dseq <- DESeqDataSetFromMatrix(counts, metadata, ~type)
stab <- vst(dseq)
head(assay(stab))

SampleDist <- dist(t(assay(stab)))
SampleDistMatrix <- as.matrix(SampleDist)
rownames(SampleDistMatrix) <- paste(stab$type, stab$id, sep = "-")
colnames(SampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(SampleDistMatrix,
         clustering_distance_rows = SampleDist,
         clustering_distance_cols = SampleDist,
         color = colors)
plotPCA(stab, intgroup = c("type"))

dseq$type <- relevel(dseq$type, "normal")
dseq <- DESeq(dseq)
res <- results(dseq)
T_Gene <- rownames(res)[which.min(res$padj)]
plotMA(res, ylim = c(-5.5,5.5))
with(res[T_Gene,], {
  points(baseMean, log2FoldChange, col = "dodgerblue", cex = 2, lwd = 2)
  text(baseMean, log2FoldChange, pos = 2, col = "pink")
})

EnhancedVolcano(res,
                lab = rownames(counts),
                x = "log2FoldChange",
                y = "pvalue",
                FCcutoff = 1.333,
                pCutoff = 10e-5,
                xlim = c(-5.7,5.7),
                ylim = c(0, -log10(10.2e-12)),
                labSize = 3,
                title = "The Results",
                subtitle = "Differential Expression",
                caption = "FCcutoff = 1.333; Pcutoff = 10e-5",
                legendPosition = "right",
                legendLabSize = 14,
                col = c("pink", "dodgerblue", "green", "yellow"),
                colAlpha = 0.7,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                hline = c(10e-8))

resord <- as.data.frame(res)
finaltable1 <- resord[order(resord$padj), ]
write.xlsx(finaltable1, "FINAL FILE.xlsx", rowNames = TRUE)

