#Title: Analisis Ekspresi Gen pada Jantung Tikus Hipertensi
#Dataset: GSE2116 (SHR vs WKY)
#Platform: Microarray (Affymetrix Rat Genome U34 - GPL85)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)


############################################
#PART A. PERSIAPAN LINGKUNGAN KERJA
############################################

# 1. Install BiocManager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 2. Install paket Bioconductor
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 
#rgu34a.db: database anotasi khusus untuk chip affymetrix rat (GSE2116)
BiocManager::install("rgu34a.db", ask = FALSE, update = FALSE)

# 3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 
install.packages(c("pheatmap", "ggplot2", "dplyr"))
#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

# 4. Memanggil library
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(rgu34a.db)
library(AnnotationDbi)
library(umap)


############################################
#PART B. PENGAMBILAN DATA DARI GEO
############################################

gset <- getGEO("GSE2116", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]


############################################
#PART C. PRE-PROCESSING DATA EKSPRESI
############################################

#log2 transformasi
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}


############################################
#PART D. DEFINISI KELOMPOK SAMPEL
############################################

group_info <- pData(gset)[["title"]]
groups <- ifelse(grepl("^SHR", group_info),
                 "SHR", "WKY")
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)


############################################
#PART E. DESIGN MATRIX (KERANGKA STATISTIK)
############################################

# 1. Design matrix
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# 2. Fit model
fit <- lmFit(ex, design)

# 3. Buat contrast
contrast_matrix <- makeContrasts(
  SHRvsWKY = SHR - WKY,
  levels = design
)

# 4. Fit contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 5. topTable
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)

head(topTableResults)


############################################
#PART F. ANOTASI NAMA GEN
############################################

probe_ids <- rownames(topTableResults)
gene_annotation <- AnnotationDbi::select(
  rgu34a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])


############################################
#PART G. VISUALISASI
############################################

# 1. Boxplot distribusi nilai ekspresi
#Boxplot digunakan untuk mengecek distribusi nilai ekspresi antar sampel.

group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

# 2. Density plot
#Density plot menunjukkan sebaran global nilai ekspresi gen

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen: SHR vs WKY (Jantung Tikus)",
    x = "Expression Value (log2)",
    y = "Density"
  )

# 3. UMAP
#UMAP digunakan untuk melihat pemisahan sampel secara global

umap_input <- t(ex)
umap_input[!is.finite(umap_input)] <- 0  

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: SHR vs WKY",
    x = "UMAP 1",
    y = "UMAP 2"
  )

# 4. Volcano plot
#Volcano plot digunakan untuk memvisualisasikan gen yang berbeda ekspresinya

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG SHR vs WKY")

# 5. Heatmap
#Heatmap digunakan untuk melihat pola ekspresi gen antar sampel

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

# 6. Menyimpan hasil analisis DEG
write.csv(topTableResults, "Hasil_GSE2116_DEG.csv")
message("Analisis selesai. File hasil telah disimpan.")


############################################
#PART H. ENRICHMENT ANALYSIS
############################################

BiocManager::install(c("clusterProfiler","org.Rn.eg.db","enrichplot","DOSE"), ask=FALSE, update=FALSE)
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(dplyr)

#Ambil DEG 
deg_genes <- topTableResults %>%
  filter(adj.P.Val < 0.05 & abs(logFC) >= 1)
gene_symbols <- deg_genes$SYMBOL[!is.na(deg_genes$SYMBOL)]
#Mapping SYMBOL -> ENTREZID
gene_df <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)
entrez_ids <- gene_df$ENTREZID

# 1. GO Enrichment
ego_bp <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Rn.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego_bp, showCategory = 15) + ggtitle("GO Biological Process Enrichment")
barplot(ego_bp, showCategory = 15, title = "GO Enrichment Barplot")

# 2. KEGG Enrichment
ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "rno",
  pvalueCutoff = 0.05
)

dotplot(ekegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")

# 3. Simpan hasil
write.csv(as.data.frame(ego_bp), "GO_BP_Enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "KEGG_Enrichment.csv", row.names = FALSE)
message("DONE.")
