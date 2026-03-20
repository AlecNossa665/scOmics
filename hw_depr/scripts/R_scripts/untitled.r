# Parameter sweep: Leiden vs Louvain, UMAP, resolution, ARI, Clustree, sc-SHC
# Uses log-normalized data (RNA assay) - NOT SCT

library(Seurat)
library(ggplot2)
library(clustree)

# Load checkpoint (seur_obj_log)
seur <- readRDS("checkpoint_after_clustree.rds")
DefaultAssay(seur) <- "RNA"

# 1. Leiden vs Louvain + resolution sweep
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
seur <- FindNeighbors(seur, dims = 1:20, verbose = FALSE)

for (res in resolutions) {
  seur <- FindClusters(seur, algorithm = 4, resolution = res, verbose = FALSE)
  seur@meta.data[[paste0("leiden_res_", res)]] <- seur@meta.data[[paste0("RNA_snn_res.", res)]]
}
for (res in resolutions) {
  seur <- FindClusters(seur, algorithm = 1, resolution = res, verbose = FALSE)
  seur@meta.data[[paste0("louvain_res_", res)]] <- seur@meta.data[[paste0("RNA_snn_res.", res)]]
}

# 2. ARI between resolutions (Leiden)
library(aricode)

n_res <- length(resolutions)
ari_leiden <- matrix(NA, n_res, n_res)
rownames(ari_leiden) <- colnames(ari_leiden) <- paste0("L", resolutions)

for (i in seq_along(resolutions)) {
  for (j in seq_along(resolutions)) {
    c1 <- seur@meta.data[[paste0("leiden_res_", resolutions[i])]]
    c2 <- seur@meta.data[[paste0("leiden_res_", resolutions[j])]]
    ari_leiden[i, j] <- ARI(c1, c2)
  }
}
print(round(ari_leiden, 3))

if (file.exists("annotation.csv")) {
  ann <- read.csv("annotation.csv", stringsAsFactors = FALSE)
  bc_col <- intersect(c("barcode", "cell", "Cell"), names(ann))[1]
  ct_col <- intersect(c("cell_type", "celltype", "CellType"), names(ann))[1]
  if (!is.na(bc_col) && !is.na(ct_col)) {
    m <- match(colnames(seur), ann[[bc_col]])
    ref <- ann[[ct_col]][m]
    ref[is.na(m)] <- "NA"
    ari_vs_ref <- sapply(resolutions, function(r) {
      ARI(seur@meta.data[[paste0("leiden_res_", r)]], ref)
    })
    names(ari_vs_ref) <- paste0("res_", resolutions)
    print("ARI vs annotation:"); print(round(ari_vs_ref, 3))
  }
}

# 3. UMAP parameter sweep (n.neighbors: 5,20,50; min.dist: 0.01, 0.1, 0.5)
n_neighbors <- c(5, 20, 50)
min_dists   <- c(0.01, 0.1, 0.5)

for (nn in n_neighbors) {
  for (md in min_dists) {
    seur <- RunUMAP(seur, dims = 1:20, n.neighbors = nn, min.dist = md,
                    reduction.name = paste0("umap_nn", nn, "_md", md),
                    reduction.key = paste0("UMAPnn", nn, "md", md, "_"),
                    verbose = FALSE)
  }
}

dir.create("umap_sweep", showWarnings = FALSE)
for (nn in n_neighbors) {
  for (md in min_dists) {
    rname <- paste0("umap_nn", nn, "_md", md)
    p <- DimPlot(seur, reduction = rname, group.by = "leiden_res_0.5") +
      ggtitle(paste0("n.neighbors=", nn, ", min.dist=", md))
    ggsave(file.path("umap_sweep", paste0(rname, ".png")), p, width = 6, height = 5, dpi = 150)
  }
}

# 4. UMAP by cluster + tissue × cluster table
seur$seurat_clusters <- seur$leiden_res_0.5

p_cluster <- DimPlot(seur, group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Clusters (Leiden res=0.5)")
ggsave("umap_clusters.png", p_cluster, width = 7, height = 5, dpi = 200)

tissue_col <- if ("Sample" %in% names(seur@meta.data)) "Sample" else "orig.ident"
tab <- table(seur@meta.data[[tissue_col]], seur@meta.data$seurat_clusters)
tab_pct <- round(100 * prop.table(tab, margin = 1), 1)
write.csv(tab, "tissue_cluster_counts.csv")
write.csv(tab_pct, "tissue_cluster_pct.csv")

# 5. Clustree
prefix <- "leiden_res_"
md <- seur@meta.data[, grep(paste0("^", prefix), names(seur@meta.data), value = TRUE)]
md <- md[, order(as.numeric(gsub(".*_([0-9.]+)$", "\\1", names(md))))]

p_tree <- clustree(md, prefix = prefix, node_colour = "# of cells",
                   node_colour_aggr = "n", node_size = 15)
ggsave("clustree_leiden.png", p_tree, width = 8, height = 6, dpi = 200)

# 6. sc-SHC significance evaluation
library(scSHC)

mat <- GetAssayData(seur, assay = "RNA", layer = "data")
cl <- as.character(seur$leiden_res_0.5)

sc_res <- testClusters(mat, cl, alpha = 0.05, num_features = 2000,
                       num_PCs = 20, parallel = TRUE, cores = 4)
seur$scSHC_clusters <- sc_res

p_shc <- DimPlot(seur, group.by = "scSHC_clusters", label = TRUE) +
  ggtitle("sc-SHC significant clusters")
ggsave("umap_scSHC.png", p_shc, width = 7, height = 5, dpi = 200)

# Save object
saveRDS(seur, "seur_sweep_results.rds")
