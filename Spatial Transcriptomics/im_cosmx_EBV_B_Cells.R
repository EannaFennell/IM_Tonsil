library(Seurat)

#######
# B cells initial clustering
#######

im_cosmx_only_ebv_b <- subset(im_cosmx_ebv_b_cells_2, idents = c("0","7"))


im_cosmx_only_ebv_b <- SCTransform(im_cosmx_only_ebv_b, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_only_ebv_b <- FindVariableFeatures(im_cosmx_only_ebv_b, selection.method = "vst", nfeatures = 150)

im_cosmx_only_ebv_b <- RunPCA(im_cosmx_only_ebv_b, npcs = 10, features = ebv.features)

ElbowPlot(im_cosmx_only_ebv_b, ndims = 10, reduction = "pca")

im_cosmx_only_ebv_b <- RunUMAP(im_cosmx_only_ebv_b, dims = 1:10)
im_cosmx_only_ebv_b <- FindNeighbors(im_cosmx_only_ebv_b, reduction = "pca", dims = 1:10)
im_cosmx_only_ebv_b <- FindClusters(im_cosmx_only_ebv_b, resolution = 0.4)

DimPlot(im_cosmx_only_ebv_b, reduction = "umap", label = T, label.size = 5, pt.size = 2)

all.markers.only.ebv.b <- FindAllMarkers(im_cosmx_only_ebv_b)

all.markers.only.ebv.b.reduced <- all.markers.only.ebv.b %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

ebv.features <- c("BALF2","BARF1","BcLF1","BCRF1","BGLF4",
                                    "BHRF1","BLLF1","BNLF2a","BNLF2b","BNRF1",
                                    "BRLF1","EBER1","EBER2","EBNA1","EBNA2",
                                    "EBNA3A","EBNA3B","EBNA3C","EBNA-LP","RPMS1/A73",
                                    "LMP1","LMP2A","LMP2A/B","BZLF1")

DotPlot(im_cosmx_only_ebv_b, features = ebv.features[!(ebv.features %in% c("EBER1","EBER2"))], dot.scale = 5) + 
  coord_flip() + 
  rotate_x_text(90) + xlab("") + ylab("") +
  viridis::scale_colour_viridis(option = "inferno") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1.1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  grids(linetype = "dashed") + ylab("")

DotPlot(im_cosmx_only_ebv_b, features = "EBNA2") + rotate_x_text(90)


im_cosmx_only_ebv_b@meta.data$clustering.base <- im_cosmx_only_ebv_b@active.ident

new_names <- c("Early","Immediate Early","Latent","Latent","Latent","Early","Early","Latent","Immediate Early")
names(new_names) <- levels(im_cosmx_only_ebv_b)
im_cosmx_only_ebv_b <- RenameIdents(object = im_cosmx_only_ebv_b, new_names)

Idents(im_cosmx_only_ebv_b) <- im_cosmx_only_ebv_b@meta.data$clustering.base


DimPlot(im_cosmx_only_ebv_b, reduction = "umap", label = T, label.size = 5, pt.size = 2)



saveRDS(im_cosmx_only_ebv_b,".../only_EBV_B_cell_2.rds")

im_cosmx_only_ebv_b <- readRDS(".../only_EBV_B_cell_2.rds")





