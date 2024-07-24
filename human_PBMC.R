library("Seurat")
filtered.counts <- Read10X_h5("/home/harith/scRNASeq/PBMC_scRNASeq/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
pbmc <- CreateSeuratObject(counts = filtered.counts, project = "Human PBMC", min.cells = 3, min.features = 200)

rm(filtered.counts)

pbmc[["percent.MT"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.MT")

pbmc <- subset(x = pbmc, subset = (nFeature_RNA > 700) & (nFeature_RNA < 5500) & (percent.MT < 10))
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
v_features_plt <- VariableFeaturePlot(pbmc)
v_features_plt
LabelPoints(plot = v_features_plt, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(pbmc)

pbmc <- RunUMAP(pbmc, dims = 1:13)
DimPlot(pbmc, reduction = "umap")

pbmc <- FindNeighbors(pbmc, dims = 1:13)

pbmc <- FindClusters(pbmc, resolution = 1.5)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.1.5")

pbmc <- FindClusters(pbmc, resolution = 1)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.1")

pbmc <- FindClusters(pbmc, resolution = 0.5)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.5")

pbmc <- FindClusters(pbmc, resolution = 0.3)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.3")

pbmc <- FindClusters(pbmc, resolution = 0.1)
DimPlot(pbmc, reduction = "umap", group.by = "RNA_snn_res.0.1")

FeaturePlot(pbmc, features = "percent.MT")
FeaturePlot(pbmc, features = "nCount_RNA")
FeaturePlot(pbmc, features = "nFeature_RNA")

markers <- read.csv("/home/harith/scRNASeq/PBMC_scRNASeq/markers_PBMC.csv")
unique_cells <- unique(markers$Cell.name)
marker_list <- list()#
for (cell in unique_cells) {
  subset_markers <- markers[markers$Cell.name == cell, ]
  cell_marker <- subset_markers$Cell.marker
  marker_list[[cell]] <- unique(cell_marker)
}

rm(markers, subset_markers, cell, cell_marker, unique_cells)

FeaturePlot(pbmc, features = marker_list$`Plasmacytoid dendritic cell`)
VlnPlot(pbmc, features = marker_list$`Plasmacytoid dendritic cell`, 
        group.by = "RNA_snn_res.0.5")
DotPlot(pbmc, features = marker_list$`Plasmacytoid dendritic cell`, 
        group.by = "RNA_snn_res.0.5")
DoHeatmap(pbmc, features = marker_list$`Plasmacytoid dendritic cell`, 
          group.by = "RNA_snn_res.0.5")

Idents(pbmc)
cluster.markers <- FindMarkers(pbmc, ident.1 = "1", ident.2 = "2")
cluster1.markers <- FindMarkers(pbmc, ident.1 = "1")
cluster.markers_roc <- FindMarkers(pbmc, ident.1 = "1", test.use = "roc")

FeaturePlot(pbmc, features = marker_list$`Memory CD4+ T cell`)
pbmc <- AddModuleScore(object = pbmc, features = list(marker_list$`Memory CD4+ T cell`), name = "tirosh_MemCD4_T")
FeaturePlot(pbmc, features = "tirosh_MemCD4_T1") 

auc_res <- AUCell::AUCell_run(exprMat = pbmc@assays$RNA$counts, geneSets = marker_list)
aucs <- AUCell::getAUC(auc_res)
pbmc <- AddMetaData(pbmc, metadata = t(aucs), 
                    col.name = paste0(rownames(aucs), "_score"))
pbmc[["newAssay"]] <- CreateAssayObject(counts = aucs)

FeaturePlot(pbmc, features = rownames(pbmc@assays$newAssay))

FeaturePlot(pbmc, features = c("B cell", 
                               "Natural killer cell", 
                               "Macrophage", 
                               "CD14+CD16+ monocyte", 
                               "T cell", 
                               "Plasmacytoid dendritic cell", 
                               "Endothelial cell", 
                               "Dendritic cell", 
                               "Fibroblast", 
                               "Cancer cell", 
                               "Pro-Natural killer cell (pro-NK cell)", 
                               "Megakaryocyte progenitor cell"))
FeaturePlot(pbmc, features = c("Naive CD4+ T cell", 
                               "Naive CD8+ T cell", 
                               "Effector CD8+ memory T (Tem) cell", 
                               "Regulatory T (Treg) cell", 
                               "Memory CD4+ T cell", 
                               "Naive T cell", 
                               "Macrophage", 
                               "Monocyte derived dendritic cell", 
                               "CD4+ T cell", 
                               "CD8+ T cell", 
                               "Natural killer T (NKT) cell"))

pbmc$curated_celltypes <- "unknown"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '0'] <- "Naive CD4+/CD8+ cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '1'] <- "PLasmacytoid Dendritic cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '2'] <- "CD14+CD16+ monocyte"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '3'] <- "Naive CD4+/CD8+ cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '4'] <- "CD8+ memory T (Tem) cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '5'] <- "Naive CD4+/CD8+ cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '6'] <- "B cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '7'] <- "Naive CD4+/CD8+ cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '8'] <- "macrophage"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '9'] <- "B cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '10'] <- "Pro-Natural killer cell (Pro-NK)"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '11'] <- "Naive CD4+/CD8+ cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '12'] <- "Monocyte derived dendritic cell"
pbmc$curated_celltypes[pbmc$RNA_snn_res.0.5 == '13'] <- "Megakaryocyte progenitor cell"

DimPlot(pbmc, group.by = "curated_celltypes", label = T)

##################################################

#Celldex and SingleR

#BiocManager::install("celldex")
#BiocManager::install("SingleR")
#imm.reference <- celldex::DatabaseImmuneCellExpressionData()
#pbmc_data_mat <- GetAssayData(pbmc, assay = "RNA", slot = "Data")
#mapping <- SingleR::SingleR(test = pbmc_data_mat, ref = imm.reference, labels = imm.reference$label.main)
#pbmc[["SingleR_labels"]] <- mapping$labels
#DimPlot(pbmc, group.by = "SingleR_labels", label = T)