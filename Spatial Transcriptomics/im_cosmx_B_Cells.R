library(Seurat)

#######
# B cells initial clustering
#######

im_cosmx_ebv_b_cells <- subset(im_cosmx, idents = c("3","4","9","10","13","14","15","16","19"))

im_cosmx_ebv_b_cells <- SCTransform(im_cosmx_ebv_b_cells, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_ebv_b_cells <- FindVariableFeatures(im_cosmx_ebv_b_cells, selection.method = "vst", nfeatures = 100)

#var.genes <- c(VariableFeatures(im_cosmx_ebv_b_cells),ebv.features)
#var.genes <- unique(var.genes[!(var.genes %in% c("EBER1","EBER2"))])

im_cosmx_ebv_b_cells <- RunPCA(im_cosmx_ebv_b_cells, npcs = 30, features = VariableFeatures(im_cosmx_ebv_b_cells))


ElbowPlot(im_cosmx_ebv_b_cells, ndims = 30, reduction = "pca")


im_cosmx_ebv_b_cells <- RunUMAP(im_cosmx_ebv_b_cells, dims = 1:30)
im_cosmx_ebv_b_cells <- FindNeighbors(im_cosmx_ebv_b_cells, reduction = "pca", dims = 1:30)
im_cosmx_ebv_b_cells <- FindClusters(im_cosmx_ebv_b_cells, resolution = 0.5)

DimPlot(im_cosmx_ebv_b_cells, reduction = "umap", label = T, label.size = 5, pt.size = 2)


all.markers.ebv.b.cells <- FindAllMarkers(im_cosmx_ebv_b_cells)

all.markers.ebv.b.cells.reduced <- all.markers.ebv.b.cells %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

#######
# EBV-infected B cells - second clustering with cleaner cells
#######

im_cosmx_ebv_b_cells_2 <- subset(im_cosmx_ebv_b_cells, idents = c("1","3","4","6","7","8","9","10","11"))

im_cosmx_ebv_b_cells_2 <- SCTransform(im_cosmx_ebv_b_cells_2, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_ebv_b_cells_2 <- FindVariableFeatures(im_cosmx_ebv_b_cells_2, selection.method = "vst", nfeatures = 100)

im_cosmx_ebv_b_cells_2 <- RunPCA(im_cosmx_ebv_b_cells_2, npcs = 30, features = VariableFeatures(im_cosmx_ebv_b_cells_2))

ElbowPlot(im_cosmx_ebv_b_cells_2, ndims = 30, reduction = "pca")

im_cosmx_ebv_b_cells_2 <- RunUMAP(im_cosmx_ebv_b_cells_2, dims = 1:30)
im_cosmx_ebv_b_cells_2 <- FindNeighbors(im_cosmx_ebv_b_cells_2, reduction = "pca", dims = 1:30)
im_cosmx_ebv_b_cells_2 <- FindClusters(im_cosmx_ebv_b_cells_2, resolution = 0.4)

DimPlot(im_cosmx_ebv_b_cells_2, reduction = "umap", label = T, label.size = 5, pt.size = 2)

all.markers.ebv.b.cells.2 <- FindAllMarkers(im_cosmx_ebv_b_cells_2)

all.markers.ebv.b.cells.reduced.2 <- all.markers.ebv.b.cells.2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

im_cosmx_ebv_b_cells_2@meta.data$clustering.base <- im_cosmx_ebv_b_cells_2@active.ident

new_names <- c("EBV+ B cells","IgG+ Plasma cells","IgA+ Plasma cells",
               "Naive cells","IgM+ Plasma cells","IgM+ Plasma cells",
               "Mature B cells","LMP1+ EBV+ B cells","TIMP1+ Plasma cells")
names(new_names) <- levels(im_cosmx_ebv_b_cells_2)
im_cosmx_ebv_b_cells_2 <- RenameIdents(object = im_cosmx_ebv_b_cells_2, new_names)

im_cosmx_ebv_b_cells_2@meta.data$clustering.base.annotation <- im_cosmx_ebv_b_cells_2@active.ident

im_cosmx_ebv_b_cells_2 <- SetIdent(im_cosmx_ebv_b_cells_2, value = im_cosmx_ebv_b_cells_2@meta.data$SCT_snn_res.0.4)

DimPlot(im_cosmx_ebv_b_cells_2, reduction = "umap", label = T, label.size = 5, pt.size = 2)

cluster.5.markers <- FindMarkers(im_cosmx_ebv_b_cells_2, ident.1 = 5, logfc.threshold = 0.05)
cluster.4.markers <- FindMarkers(im_cosmx_ebv_b_cells_2, ident.1 = 4, logfc.threshold = 0.05)

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("LMP1"), pt.size = 2)

lmp1_vs_ebv <- FindMarkers(im_cosmx_ebv_b_cells_2, ident.1 = "LMP1+ EBV+ B cells", ident.2 = "EBV+ B cells", logfc.threshold = 0)
write.csv(lmp1_vs_ebv,".../lmp1_vs_ebv_dge.csv", row.names = T)

lmp1_vs_naive <- FindMarkers(im_cosmx_ebv_b_cells_2, ident.1 = "LMP1+ EBV+ B cells", ident.2 = "Naive cells", logfc.threshold = 0)

lmp1_vs_plasma <- FindMarkers(im_cosmx_ebv_b_cells_2, ident.1 = "LMP1+ EBV+ B cells", ident.2 = c("IgM+ Plasma cells","IgG+ Plasma cells","IgA+ Plasma cells"), logfc.threshold = 0)




DimPlot(im_cosmx_ebv_b_cells_2, reduction = "umap", pt.size = 1.2,
        cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40")) 


#F9CF40

#######
# EBV-infected B cells - ebv positive per cluster
#######

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("EBER1"), pt.size = 2, min.cutoff = "q50")

cut.off <- c(3.5,2)
table(as.matrix(im_cosmx_ebv_b_cells_2@assays$SCT@data)["EBER1",] > cut.off[1] & as.matrix(im_cosmx_ebv_b_cells_2@assays$SCT@data)["EBER2",] > cut.off[2])

im_cosmx_ebv_b_cells_2@meta.data$ebv_positive <- as.matrix(im_cosmx_ebv_b_cells_2@assays$SCT@data)["EBER1",] > cut.off[1] & as.matrix(im_cosmx_ebv_b_cells_2@assays$SCT@data)["EBER2",] > cut.off[2]


FeaturePlot(im_cosmx_ebv_b_cells_2, features = "ebv_positive", pt.size = 1, order = TRUE)

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("CD274"), pt.size = 2)

DotPlot(im_cosmx_ebv_b_cells_2, features = c("IGHG1","IGHG2","IGHM","IGHA1","IGHD","JCHAIN",
                                             "MS4A1","MZB1","IGKC","CD19","XBP1","IRF4","CD79A",
                                             "CD24","CD38","IL7R","MYC","CXCR5","CD40","CD40LG","CCR7",
                                             "RAG1","IFNG","CXCR4","BCL2","TIMP1","TNFRSF9","ACKR4","LAIR1","LTB","TCL1A")) + rotate_x_text(90)

DotPlot(im_cosmx_ebv_b_cells_2, features = c("IGHG1","IGHG2","IGHM","IGHA1","IGHD","JCHAIN",
                                             "MS4A1","MZB1","IGKC","CD19","XBP1","IRF4","CD79A",
                                             "CD24","CD38","IL7R","MYC","CXCR5","CD40","CD40LG","CCR7",
                                             "RAG1","IFNG","CXCR4","BCL2","TIMP1","TNFRSF9","ACKR4","LAIR1","LTB","TCL1A"), cols = c("#F9CF40","#74C242"), dot.scale = 5) + 
  rotate_x_text(90) + xlab("") + ylab("") +
  viridis::scale_colour_viridis(option = "viridis") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1.1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  grids(linetype = "dashed") + ylab("")




saveRDS(im_cosmx_ebv_b_cells_2,".../base_B_cell.rds")

im_cosmx_ebv_b_cells_2 <- readRDS(".../base_B_cell.rds")

#######
# EBV-infected B cells - lytic module score
#######

lytic_genes <- list(c("BALF2","BARF1","BcLF1","BCRF1","BGLF4",
                      "BHRF1","BLLF1","BNLF2a","BNLF2b","BNRF1",
                      "BRLF1","RPMS1/A73","BZLF1"))
im_cosmx_ebv_b_cells_2 <- AddModuleScore(
  object = im_cosmx_ebv_b_cells_2,
  features = lytic_genes,
  ctrl = 5,
  name = 'lytic_score'
)

FeaturePlot(im_cosmx_ebv_b_cells_2, features = "lytic_score1", min.cutoff = "q2")

VlnPlot(im_cosmx_ebv_b_cells_2, features = "lytic_score1", sort = "increasing", pt.size = 0,
        cols = c("#364CA2","#74C242","#B89365","#AD3496","#4CC7EF","#F9CF40","#18A58A","#F04037")) + theme(legend.position = "none") + labs(title="",
                                                                                                                                      x ="", y = "EBV Lytic Score")

VlnPlot(im_cosmx_ebv_b_cells_2, features = "LMP1", sort = "increasing", pt.size = 0,
        cols = c("#74C242","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#364CA2","#F9CF40")) + theme(legend.position = "none") + labs(title="",
                                                                                                                                      x ="", y = "LMP1 Expression")


VlnPlot(im_cosmx_ebv_b_cells_2, features = "LMNA", sort = "increasing", pt.size = 0,
        cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40")) + theme(legend.position = "none") + labs(title="",
                                                                                                                              x ="", y = "LMNA Expression") + coord_flip()
                                                                                                                                      
im_cosmx_ebv_b_cells_2 <- SetIdent(im_cosmx_ebv_b_cells_2, value = im_cosmx_ebv_b_cells_2@meta.data$clustering.base.annotation)




library(viridis)

DotPlot(im_cosmx_ebv_b_cells_2, features =ebv.features[!(ebv.features %in% c("EBER1","EBER2"))], cols = c("#F9CF40","#74C242"), dot.scale = 6) + 
  rotate_x_text(90) + xlab("") + ylab("") +
  viridis::scale_colour_viridis(option = "magma") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=1.1) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  grids(linetype = "dashed") + ylab("") 


#######
# LMP1 and LMP2A/LMP2B co-expression
#######

im_cosmx_ebv_b_cells_2@active.assay <- "SCT"

rownames(im_cosmx_ebv_b_cells_2@assays$SCT@scale.data)[rownames(im_cosmx_ebv_b_cells_2@assays$SCT@scale.data) == "LMP2A/B"] <- "LMP2B"

im_cosmx_ebv_b_cells_2 <- readRDS(".../base_B_cell.rds")

FeatureScatter(im_cosmx_ebv_b_cells_2, feature1 = "LMP1", feature2 = "LMP2A", slot = "scale.data",
               cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40"))

FeatureScatter(im_cosmx_ebv_b_cells_2, feature1 = "LMP1", feature2 = "LMP2B", slot = "scale.data",
               cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40")) + ylab("LMP2A/B")




FeatureScatter(im_cosmx_ebv_b_cells_2, feature1 = "EBNA2", feature2 = "EBNA3C", slot = "scale.data",
               cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40"))


#######
# MIF, T-bet, CD11c and CXCR3
#######

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("TBX21"), pt.size = 0.5)
VlnPlot(im_cosmx_ebv_b_cells_2, features = c("TBX21"))

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("ITGAX"), pt.size = 0.5)

FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("CXCR3"), pt.size = 0.5)


FeaturePlot(im_cosmx_ebv_b_cells_2, features = c("MIF"), pt.size = 0.5)



im_cosmx_ebv_b_cells_2@active.assay <- "RNA"

VlnPlot(im_cosmx_ebv_b_cells_2, features = c("MIF"), idents = c("0","7"), pt.size = 0) + xlab("") + ggdist::theme_ggdist() + theme(legend.position = "none",
                                                                                                                                   axis.text.x = element_text(size = 18, colour = "black"),
                                                                                                                                   axis.text.y = element_text(size = 18, colour = "black"),
                                                                                                                                   axis.title.y = element_text(size = 20, colour = "black"),
                                                                                                                                   axis.ticks.length=unit(.25, "cm"),
                                                                                                                                   axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                   axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                   axis.ticks = element_line(color="black")) + 
  stat_compare_means(comparisons = list(c("0","7"))) + ggtitle("") + ylab("MIF Expression") + scale_fill_manual(values = c("#74C242","#364CA2"))



VlnPlot(im_cosmx_ebv_b_cells_2, features = c("MIF"), pt.size = 0, cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A","#B89365","#F04037","#74C242","#F9CF40"))
