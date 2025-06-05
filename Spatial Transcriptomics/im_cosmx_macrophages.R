library(Seurat)


#######
# Macrophages initial clustering
#######

Idents(object = im_cosmx) <- im_cosmx@meta.data$seurat_clusters

im_cosmx_macs <- subset(im_cosmx, idents = c("6","8"))

im_cosmx_macs <- SCTransform(im_cosmx_macs, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_macs <- FindVariableFeatures(im_cosmx_macs, selection.method = "vst", nfeatures = 100)

im_cosmx_macs <- RunPCA(im_cosmx_macs, npcs = 30, features = VariableFeatures(im_cosmx_macs))


ElbowPlot(im_cosmx_macs, ndims = 30, reduction = "pca")


im_cosmx_macs <- RunUMAP(im_cosmx_macs, dims = 1:30)
im_cosmx_macs <- FindNeighbors(im_cosmx_macs, reduction = "pca", dims = 1:30)
im_cosmx_macs <- FindClusters(im_cosmx_macs, resolution = 0.5)

DimPlot(im_cosmx_macs, reduction = "umap", label = T, label.size = 5, pt.size = 2)


all.markers.macs <- FindAllMarkers(im_cosmx_macs)

all.markers.macs.reduced <- all.markers.macs %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)




#######
# Macrophages clean clustering
#######

im_cosmx_macs_2 <- subset(im_cosmx_macs, idents = c("1","3","6","8","9","10"), invert = TRUE)

im_cosmx_macs_2 <- SCTransform(im_cosmx_macs_2, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_macs_2 <- FindVariableFeatures(im_cosmx_macs_2, selection.method = "vst", nfeatures = 50)

im_cosmx_macs_2 <- RunPCA(im_cosmx_macs_2, npcs = 30, features = VariableFeatures(im_cosmx_macs_2))


ElbowPlot(im_cosmx_macs_2, ndims = 30, reduction = "pca")


im_cosmx_macs_2 <- RunUMAP(im_cosmx_macs_2, dims = 1:30)
im_cosmx_macs_2 <- FindNeighbors(im_cosmx_macs_2, reduction = "pca", dims = 1:30)
im_cosmx_macs_2 <- FindClusters(im_cosmx_macs_2, resolution = 0.6)

DimPlot(im_cosmx_macs_2, reduction = "umap", pt.size = 1, cols = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88"))



im_cosmx_macs_2 <- Seurat::ScaleData(im_cosmx_macs_2, features = rownames(im_cosmx_macs_2))




all.markers.macs.2 <- FindAllMarkers(im_cosmx_macs_2)

all.markers.macs.reduced.2 <- all.markers.macs.2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)


FeaturePlot(im_cosmx_macs_2, features = "CD274", pt.size = 1)

VlnPlot(im_cosmx_macs_2, features = "CD163", pt.size = 0)

saveRDS(im_cosmx_macs_2, ".../base_macs.rds")



im_cosmx_macs_2 <- readRDS(".../base_macs.rds")

mac.features <- c("CCL18","IL2RA","GPNMB","CXCL12","C1QB","CCL5","MMP12","MMP9","SPP1",
                  "S100A8","S100A9","IL1B","DUSP1","CXCL10","CXCL9","CD44","MT2A","CCL8","MT1X","CCL2","ITGAX",
                  "IDO1","CD74","CIITA","HLA-DQB1/2","HSPA1A/B","THBS1","IGFBP5","COL3A1","TIMP1")

DoHeatmap(im_cosmx_macs_2, features = mac.features, disp.min = -2.5, disp.max = 2, group.colors = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88"), angle = 0) +
  scale_fill_distiller(palette ="RdBu", direction = -1)

FeaturePlot(im_cosmx_macs_2, features = "CD74", pt.size = 1)



