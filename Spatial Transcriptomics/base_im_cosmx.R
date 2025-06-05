library(Seurat)
library(ggpubr)

##########
# Import and cluster
##########

im_cosmx <- readRDS(".../base_cleaned_IM.rds")

im_cosmx <- SCTransform(im_cosmx, assay = "RNA", clip.range = c(-10, 10))

im_cosmx <- FindVariableFeatures(im_cosmx, selection.method = "vst", nfeatures = 800)

im_cosmx <- RunPCA(im_cosmx, npcs = 50, features = VariableFeatures(im_cosmx))

ElbowPlot(im_cosmx, ndims = 50, reduction = "pca")

im_cosmx <- RunUMAP(im_cosmx, dims = 1:50)
im_cosmx <- FindNeighbors(im_cosmx, reduction = "pca", dims = 1:50)
im_cosmx <- FindClusters(im_cosmx, resolution = 0.5) 

DimPlot(im_cosmx, reduction = "umap", label = T, label.size = 5)

im_cosmx@meta.data$clustering.base <- im_cosmx@active.ident
im_cosmx@meta.data$clustering.high.res <- im_cosmx@active.ident
#im_cosmx<- SetIdent(im_cosmx, value = im_cosmx@meta.data$seurat_clusters)

features <- c("CD3E","CD3G","CD3D","MS4A1","CD19","CD4",
              "CD8A","CD68","CD14","CD163","SPP1","MKI67",
              "EBNA2","LMP1","CD79A","KRT8",
              "KRT18","JCHAIN")
library(viridis)

fig.1.features <- rev(c("GZMB","COTL1","CD3D","NKG7","LYZ",
                        "C1QC","CXCL9","IGHG1","IGHG2","XBP1","MZB1",
                        "CXCL8","HCAR2/3","TIMP1","COL6A2","COL1A1",
                        "COL3A1","FN1","RPMS1/A73","COL4A1","IGFBP7",
                        "SPARCL1","HSP90AA1","HSPA1A/B","HSPB1",
                        "KRT14","KRT17","KRT6A/B/C","KRT5"))

DotPlot(object = im_cosmx, features = fig.1.features, idents = c("T cells","Macrophages",
                                                           "Plasma/B cells","CXCL8+",
                                                           "Fibroblasts","EBV-infected B cells",
                                                           "Endothelial cells","Glandular cells",
                                                           "Epithelial cells")) + theme_bw() + 
  rotate_x_text(45) + coord_flip() + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") # 5.8 x 7



DoHeatmap(object = im_cosmx)

mac.features <- c("CCL18","IL2RA","GPNMB","CXCL12","C1QB","CCL5","MMP12","MMP9","SPP1",
                  "S100A8","S100A9","IL1B","DUSP1","CXCL10","CXCL9","CD44","MT2A","CCL8","MT1X","CCL2","ITGAX",
                  "IDO1","CD74","CIITA","HLA-DQB1/2","HSPA1A/B","THBS1","IGFBP5","COL3A1","TIMP1")

DoHeatmap(im_cosmx, features = mac.features, 
          disp.min = -2.5, disp.max = 2,
          group.colors = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88"), angle = 0) +
  scale_fill_distiller(palette ="RdBu", direction = -1)


saveRDS(im_cosmx,".../base_cleaned_IM_clustered.rds")

im_cosmx <- readRDS(".../base_cleaned_IM_clustered.rds")

##########
# Annotating clusters
##########

all.markers <- FindAllMarkers(im_cosmx)

all.markers.reduced <- all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.reduced, ".../all_markers_reduced_2.csv", row.names = F)

cluster.1.markers <- FindMarkers(im_cosmx, ident.1 = 1, logfc.threshold = 0.05)
cluster.9.markers <- FindMarkers(im_cosmx, ident.1 = 9, ident.2 = c(0), logfc.threshold = 0.05)
cluster.11.markers <- FindMarkers(im_cosmx, ident.1 = 11, logfc.threshold = 0.05)

new_names <- c("Unclassified", "T cells", "Unclassified", "Macrophages", "Plasma/B cells", 
               "CXCL8+", "Plasma/B cells", "Fibroblasts",
               "Plasma/B cells","EBV-infected B cells","EBV-infected B cells",
               "Endothelial cells","Glandular cells","Epithelial cells",
               "Plasma/B cells","EBV-infected B cells")

names(new_names) <- levels(im_cosmx)

im_cosmx <- RenameIdents(object = im_cosmx, new_names)

im_cosmx@meta.data$clustering.annotations <- im_cosmx@active.ident

im_cosmx<- SetIdent(im_cosmx, value = im_cosmx@meta.data$seurat_clusters)


pdf(file = ".../umap_2.pdf",   # The directory you want to save the file in
    width = 5.2, # The width of the plot in inches
    height = 5)

DimPlot(im_cosmx, reduction = "umap", label = T, label.size = 5,
        cols = c("#B8BAC0","#078FD3","#DC4B00","#94CE95","#006178","#FE9075","#ABDFEF","#FFD7B4","#C44C68","black")) + theme(legend.position = "none")

# Step 3: Run dev.off() to create the file!
dev.off()


FeaturePlot(object = im_cosmx, features = "KRT16")


##########
# EBV expression
##########

ebv.features <- c("BALF2","BARF1","BcLF1","BCRF1","BGLF4",
                  "BHRF1","BLLF1","BNLF2a","BNLF2b","BNRF1",
                  "BRLF1","EBER1","EBER2","EBNA1","EBNA2",
                  "EBNA3A","EBNA3B","EBNA3C","EBNA-LP","RPMS1/A73",
                  "LMP1","LMP2A","LMP2A/B","BZLF1")

ebv.features.reduced <- ebv.features[c(12:15,21:22)]

FeaturePlot(object = im_cosmx, features = ebv.features.reduced, ncol = 6)

FeaturePlot(object = im_cosmx, features = c("EBER1","EBER2",
                                            "EBNA1","EBNA2",
                                            "EBNA3A","EBNA-LP"), min.cutoff = "q70", 
            raster = TRUE, pt.size = 2, label.size = 0, ncol = 2) # 6 x 7 


FeaturePlot(object = im_cosmx, features = c("LMP1","LMP2A",
                                            "RPMS1/A73","BZLF1",
                                            "BHRF1","BNLF2b"), min.cutoff = "q70", 
            raster = TRUE, pt.size = 2, label.size = 0, ncol = 3) 


VlnPlot(im_cosmx, features = ebv.features.reduced)

##########
# Checking spatial location of clusters
##########

# y needs to be flipped

cluster <- 4
fov <- 54  #17 18 19 57

cluster.4.location <- data.frame(x1 = im_cosmx@meta.data[["x_FOV_px"]][im_cosmx@meta.data$fov == fov & im_cosmx@meta.data$seurat_clusters == cluster],
                                    y1 = im_cosmx@meta.data[["y_FOV_px"]][im_cosmx@meta.data$fov == fov & im_cosmx@meta.data$seurat_clusters == cluster])

ggpubr::ggscatter(data = cluster.4.location, x = "x1", y = "y1")

cluster.location <- data.frame(x1 = im_cosmx@meta.data[["x_FOV_px"]][im_cosmx@meta.data$fov == fov],
                                 y1 = im_cosmx@meta.data[["y_FOV_px"]][im_cosmx@meta.data$fov == fov])

ggpubr::ggscatter(data = cluster.location, x = "x1", y = "y1")


##########
# Epithelial cell re-clustering
##########

im_cosmx_epi <- subset(im_cosmx, idents = c(12))

#im_cosmx_epi <- subset(im_cosmx_epi, idents = 2, invert = TRUE)

im_cosmx_epi <- SCTransform(im_cosmx_epi, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_epi <- FindVariableFeatures(im_cosmx_epi, selection.method = "vst", nfeatures = 200)

#var.genes <- c(VariableFeatures(im_cosmx_epi),ebv.features)
#var.genes <- unique(var.genes[!(var.genes %in% c("EBER1","EBER2"))])

im_cosmx_epi <- RunPCA(im_cosmx_epi, npcs = 50, features = VariableFeatures(im_cosmx_epi))


ElbowPlot(im_cosmx_epi, ndims = 50, reduction = "pca")


im_cosmx_epi <- RunUMAP(im_cosmx_epi, dims = 1:50)
im_cosmx_epi <- FindNeighbors(im_cosmx_epi, reduction = "pca", dims = 1:50)
im_cosmx_epi <- FindClusters(im_cosmx_epi, resolution = 0.3)

DimPlot(im_cosmx_epi, reduction = "umap", label = T, label.size = 5, pt.size = 2, cols = c("#031D44","#70A288","#DAB785"))

DimPlot(im_cosmx_epi, reduction = "umap", label = T, label.size = 5, pt.size = 2, group.by = "donor")

DimPlot(im_cosmx_epi, reduction = "umap", label = T, label.size = 5, pt.size = 2, group.by = "fov")



FeaturePlot(object = im_cosmx_epi, features = "EBER1", min.cutoff = "q25", pt.size = 2)

FeaturePlot(object = im_cosmx_epi, features = "EBER1", pt.size = 2)

FeaturePlot(object = im_cosmx_epi, features = ebv.features.reduced, ncol = 6)


DotPlot(im_cosmx_epi, features= ebv.features[!(ebv.features %in% c("EBER1","EBER2"))]) + theme_bw() + 
  rotate_x_text(45) + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma")






all.markers.epi <- FindAllMarkers(im_cosmx_epi)

all.markers.epi.reduced <- all.markers.epi %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.1)

cluster.2.markers.epi <- FindMarkers(im_cosmx, ident.1 = 2, logfc.threshold = 0.05)
cluster.3.markers.epi <- FindMarkers(im_cosmx, ident.1 = 3, logfc.threshold = 0.05)
cluster.5.markers.epi <- FindMarkers(im_cosmx, ident.1 = 5, logfc.threshold = 0.05)

VlnPlot(im_cosmx_epi, idents = 2, features = ebv.features[c(1:12)])

VlnPlot(im_cosmx_epi, features = ebv.features[c(1:12)])

VlnPlot(im_cosmx_epi, features = ebv.features[c(13:24)])


DotPlot(im_cosmx_epi, features= c("CD44","KRT14","KRT13","KRT4")) + rotate_x_text(45)

saveRDS(im_cosmx_epi, ".../base_cleaned_IM_epi_clustered.rds")

im_cosmx_epi <- readRDS(".../base_cleaned_IM_epi_clustered.rds")

DotPlot(im_cosmx_epi, features= ebv.features) + rotate_x_text(45)

cluster <- c(0,1,2)
fov <- 10  #17 18 19 57
marker <- "LMP2A"

cluster.4.location <- data.frame(x1 = im_cosmx_epi@meta.data[["x_FOV_px"]][im_cosmx_epi@meta.data$fov == fov & im_cosmx_epi@meta.data$seurat_clusters %in% cluster],
                                 y1 = im_cosmx_epi@meta.data[["y_FOV_px"]][im_cosmx_epi@meta.data$fov == fov & im_cosmx_epi@meta.data$seurat_clusters %in% cluster],
                                 col1 = im_cosmx_epi@meta.data$seurat_clusters[im_cosmx_epi@meta.data$fov == fov & im_cosmx_epi@meta.data$seurat_clusters %in% cluster],
                                 exp = as.matrix(im_cosmx_epi@assays$RNA@data)[marker,][im_cosmx_epi@meta.data$fov == fov & im_cosmx_epi@meta.data$seurat_clusters %in% cluster])

ggpubr::ggscatter(data = cluster.4.location, x = "x1", y = "y1", color = "col1", palette = c("#031D44","#70A288","#DAB785"))


new_names <- c("Columnar", "Polygonal", "Squamous")
names(new_names) <- levels(im_cosmx_epi)

im_cosmx_epi <- RenameIdents(object = im_cosmx_epi, new_names)


FeaturePlot(im_cosmx_epi, features = c("KRT4"))



DotPlot(im_cosmx_epi, features= c("CD44","KRT14","KRT13","KRT4")) + theme_bw() + 
  rotate_x_text(45) + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma")


##############
# Look in EBV+ only

Idents(object = im_cosmx_epi) <- "ebv_positive"

im_cosmx_epi_pos <- subset(im_cosmx_epi, idents = "TRUE")

im_cosmx_epi_pos@meta.data$ebv_positive

Idents(object = im_cosmx_epi_pos) <- "seurat_clusters"

DotPlot(im_cosmx_epi_pos, features= ebv.features[!(ebv.features %in% c("EBER1","EBER2"))]) + theme_bw() + 
  rotate_x_text(45) + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma")

DotPlot(im_cosmx_epi_pos, features= c("EBER1","EBER2")) + theme_bw() + 
  rotate_x_text(45) + xlab("") + ylab("") + 
  theme(axis.text = element_text(size = 14, colour = "black")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) 


###################################
# Epithelial cell analysis for paper
###################################

# How many cells are EBV? EBER1 & EBER2 > 1
cut.off <- 1
table(as.matrix(im_cosmx_epi@assays$SCT@data)["EBER1",] > cut.off & as.matrix(im_cosmx_epi@assays$SCT@data)["EBER2",] > cut.off)

im_cosmx_epi@meta.data$ebv_positive <- as.matrix(im_cosmx_epi@assays$SCT@data)["EBER1",] > cut.off & as.matrix(im_cosmx_epi@assays$SCT@data)["EBER2",] > cut.off

overall_percent_positive <- data.frame(class = c("Negative","Positive"),
                                    percent_positive = c((table(im_cosmx_epi@meta.data$ebv_positive)[1] / sum(table(im_cosmx_epi@meta.data$ebv_positive)))*100, (table(im_cosmx_epi@meta.data$ebv_positive)[2] / sum(table(im_cosmx_epi@meta.data$ebv_positive)))*100))

ggpubr::ggbarplot(overall_percent_positive, x = "class", y = "percent_positive", fill = "class") + scale_y_continuous(expand = c(0, 0)) + xlab("") + ylab("Percentage Positive") + theme(legend.position = "none")



percent_positive_by_diff <- data.frame(class = c("Col-Negative","Col-Positive","Poly-Negative","Poly-Positive","Squa-Negative","Squa-Positive"),
                                       percent_positive = c((table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[1] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100, 
                                                            (table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[2] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100,
                                                            (table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[1] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100, 
                                                            (table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[2] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100,
                                                            (table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[1] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100, 
                                                            (table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[2] / sum(table(im_cosmx_epi@meta.data$ebv_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100))

percent_positive_by_diff_sub <- percent_positive_by_diff[c(5,6),]

ggpubr::ggbarplot(percent_positive_by_diff_sub, x = "class", 
                  y = "percent_positive", fill = "class", palette = "jama") + 
  scale_y_continuous(expand = c(0, 0)) + xlab("") + 
  ylab("Percentage Positive") + theme(legend.position = "none") + rotate_x_text(45)


# How many cells are lytic and latent?

cut.off <- 0
im_cosmx_epi@meta.data$lytic_positive <- as.matrix(im_cosmx_epi@assays$SCT@data)["BZLF1",] > cut.off | as.matrix(im_cosmx_epi@assays$SCT@data)["BRLF1",] > cut.off

lytic_positive_by_diff <- data.frame(class = c("Col-Negative","Col-Positive","Poly-Negative","Poly-Positive","Squa-Negative","Squa-Positive"),
                                       percent_positive = c((table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[1] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100, 
                                                            (table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[2] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100,
                                                            (table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[1] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100, 
                                                            (table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[2] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100,
                                                            (table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[1] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100, 
                                                            (table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[2] / sum(table(im_cosmx_epi@meta.data$lytic_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100))

lytic_positive_by_diff_sub <- lytic_positive_by_diff[c(5,6),]

ggpubr::ggbarplot(lytic_positive_by_diff_sub, x = "class", y = "percent_positive", fill = "class") + scale_y_continuous(expand = c(0, 0)) + xlab("") + 
  ylab("Percentage Positive") + theme(legend.position = "none") + rotate_x_text(45)






im_cosmx_epi@meta.data$latent_positive <- as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA2",] > cut.off | as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA3A",] > cut.off | 
                                          as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA3B",] > cut.off | as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA3C",] > cut.off |
                                          as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA-LP",] > cut.off | as.matrix(im_cosmx_epi@assays$SCT@data)["LMP1",] > cut.off |
                                          as.matrix(im_cosmx_epi@assays$SCT@data)["LMP2A",] > cut.off | as.matrix(im_cosmx_epi@assays$SCT@data)["LMP2A/B",] > cut.off

latent_positive_by_diff <- data.frame(class = c("Col-Negative","Col-Positive","Poly-Negative","Poly-Positive","Squa-Negative","Squa-Positive"),
                                     percent_positive = c((table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[1] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100, 
                                                          (table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])[2] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 0])))*100,
                                                          (table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[1] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100, 
                                                          (table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])[2] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 1])))*100,
                                                          (table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[1] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100, 
                                                          (table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])[2] / sum(table(im_cosmx_epi@meta.data$latent_positive[im_cosmx_epi@meta.data$seurat_clusters == 2])))*100))

latent_positive_by_diff_sub <- latent_positive_by_diff[c(5,6),]

ggpubr::ggbarplot(latent_positive_by_diff_sub, x = "class", y = "percent_positive", fill = "class") + scale_y_continuous(expand = c(0, 0)) + xlab("") + 
  ylab("Percentage Positive") + theme(legend.position = "none") + rotate_x_text(45)


# DEG of EBV+ and EBV- epithelial cells

im_cosmx_epi@meta.data$epi_annotation <- im_cosmx_epi@active.ident

im_cosmx_epi <- SetIdent(im_cosmx_epi, value = "ebv_positive")

ebv_pos_vs_neg_epi <- FindMarkers(im_cosmx_epi, ident.1 = TRUE, ident.2 = FALSE, logfc.threshold = 0.05)

write.csv(x = ebv_pos_vs_neg_epi, file = "D:/IM_Tonsil_Paper/CosMx/ebv_pos_epi_vs_ebv_neg_epi.csv", sep = ",", quote = F, row.names = T, col.names = T)

ebv_pos_vs_neg_epi <- read.csv(".../ebv_pos_epi_vs_ebv_neg_epi.csv")

library(EnhancedVolcano)

keyvals <- rep("darkgrey", times = nrow(ebv_pos_vs_neg_epi))
names(keyvals) <- "All"
keyvals[ebv_pos_vs_neg_epi$X %in% c("EBER2","EBER1")] <- "red"
names(keyvals)[ebv_pos_vs_neg_epi$X %in% c("EBER2","EBER1")] <- "EBV"
keyvals[ebv_pos_vs_neg_epi$X %in% c("IGHG1","IGHG2")] <- "blue"
names(keyvals)[ebv_pos_vs_neg_epi$X %in% c("IGHG1","IGHG2")] <- "IgG"

EnhancedVolcano::EnhancedVolcano(ebv_pos_vs_neg_epi,
                                 lab = ebv_pos_vs_neg_epi$X,
                                 x = 'avg_log2FC',
                                 y = 'p_val',
                                 pCutoff = 10e-3,
                                 FCcutoff = 0.5,
                                 pointSize = 3.0,
                                 labSize = 4,
                                 colCustom = keyvals,
                                 selectLab = c("EBER2","EBER1","IGHG1","IGHG2","MMP12")) + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 12, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12, face = "plain"),
        axis.ticks.length=unit(.2, "cm"))


# col = c('darkgrey', 'darkgrey', 'darkgrey',"darkgrey")

# Co-expression analysis

cor(as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA1",],  as.matrix(im_cosmx_epi@assays$SCT@data)["EBNA2",])

# Monocle

library(monocle3)
library(Seurat)
library(ggplot2)
source("D:/NPC/T_cell_paper/R_Scripts/asCellDataSet.R")

im_cosmx_epi_cds <- as.cell_data_set(im_cosmx_epi)

im_cosmx_epi_cds <- estimate_size_factors(im_cosmx_epi_cds) # necessary for running plot_gene_in_pseudotime

fData(im_cosmx_epi_cds)$gene_short_name <- rownames(fData(im_cosmx_epi_cds))
head(fData(im_cosmx_epi_cds))

head(counts(im_cosmx_epi_cds))

# Retrieve clustering information from Seurat object

recreate.partitions <- c(rep(1, length(im_cosmx_epi_cds@colData@rownames)))
names(recreate.partitions) <- im_cosmx_epi_cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

im_cosmx_epi_cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- im_cosmx_epi@active.ident
im_cosmx_epi_cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

im_cosmx_epi_cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- im_cosmx_epi@reductions$umap@cell.embeddings

# Plot

cluster.before.traj <- plot_cells(im_cosmx_epi_cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                  group_label_size = 5, cell_size = 2) + theme(legend.position = "right") + scale_color_manual(values = c("#031D44","#70A288","#DAB785"))

cluster.before.traj

# Learn trajectory - adjust partitions

im_cosmx_epi_cds@clusters@listData$UMAP$partitions <- im_cosmx_epi_cds@clusters@listData$UMAP$clusters
levels(im_cosmx_epi_cds@clusters@listData$UMAP$partitions) <- c("1","2","3") # Needs to be 1,2,3,etc. for whatever reason
levels(im_cosmx_epi_cds@clusters@listData$UMAP$partitions)

im_cosmx_epi_cds <- learn_graph(im_cosmx_epi_cds, use_partition = FALSE)

im_cosmx_epi_cds <- order_cells(im_cosmx_epi_cds, reduction_method = "UMAP")

plot_cells(im_cosmx_epi_cds, color_cells_by = "pseudotime", label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 2)



ebv.genes <- ebv.features[c(17:24)]
im_cosmx_epi_cds_ds <- im_cosmx_epi_cds[rowData(im_cosmx_epi_cds)$gene_short_name %in% ebv.genes,
                                        colData(im_cosmx_epi_cds)$ebv_positive %in% c("TRUE")]
im_cosmx_epi_cds_ds <- order_cells(im_cosmx_epi_cds_ds)

plot_genes_in_pseudotime(im_cosmx_epi_cds_ds,
                         color_cells_by="pseudotime",
                         min_expr=0.01, cell_size = 1)

saveRDS(im_cosmx_epi_cds, ".../epithelium_monocle.rds")



