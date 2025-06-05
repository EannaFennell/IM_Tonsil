library(Seurat)
library(dplyr)
library(viridis)

#######
# T cells - initial clustering
#######

im_cosmx_t_cells <- subset(im_cosmx, idents = c("1"))

im_cosmx_t_cells <- SCTransform(im_cosmx_t_cells, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_t_cells <- FindVariableFeatures(im_cosmx_t_cells, selection.method = "vst", nfeatures = 100)

#var.genes <- c(VariableFeatures(im_cosmx_t_cells),ebv.features)
#var.genes <- unique(var.genes[!(var.genes %in% c("EBER1","EBER2"))])

im_cosmx_t_cells <- RunPCA(im_cosmx_t_cells, npcs = 30, features = VariableFeatures(im_cosmx_t_cells))


ElbowPlot(im_cosmx_t_cells, ndims = 30, reduction = "pca")


im_cosmx_t_cells <- RunUMAP(im_cosmx_t_cells, dims = 1:30)
im_cosmx_t_cells <- FindNeighbors(im_cosmx_t_cells, reduction = "pca", dims = 1:30)
im_cosmx_t_cells <- FindClusters(im_cosmx_t_cells, resolution = 0.5)


DimPlot(im_cosmx_t_cells, reduction = "umap", label = T, label.size = 5, pt.size = 2)

all.markers.t.cells <- FindAllMarkers(im_cosmx_t_cells)

all.markers.t.cells.reduced <- all.markers.t.cells %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.1)

FeaturePlot(object = im_cosmx_t_cells, features = ebv.features.reduced, ncol = 6)


DotPlot(im_cosmx_t_cells, features = c("CD3E","CD4","CD8A","IL27RA","FOXP3","NKG7"))


#######
# T cells - clean clustering
#######

im_cosmx_t_cells_2 <- subset(im_cosmx_t_cells, idents = c("1","5","7","10"), invert = TRUE)

im_cosmx_t_cells_2 <- SCTransform(im_cosmx_t_cells_2, assay = "RNA", clip.range = c(-10, 10))

im_cosmx_t_cells_2 <- FindVariableFeatures(im_cosmx_t_cells_2, selection.method = "vst", nfeatures = 50)

#var.genes <- c(VariableFeatures(im_cosmx_t_cells),ebv.features)
#var.genes <- unique(var.genes[!(var.genes %in% c("EBER1","EBER2"))])

im_cosmx_t_cells_2 <- RunPCA(im_cosmx_t_cells_2, npcs = 30, features = VariableFeatures(im_cosmx_t_cells))

ElbowPlot(im_cosmx_t_cells_2, ndims = 30, reduction = "pca")

im_cosmx_t_cells_2 <- RunUMAP(im_cosmx_t_cells_2, dims = 1:30)
im_cosmx_t_cells_2 <- FindNeighbors(im_cosmx_t_cells_2, reduction = "pca", dims = 1:30)
im_cosmx_t_cells_2 <- FindClusters(im_cosmx_t_cells_2, resolution = 0.5)


DimPlot(im_cosmx_t_cells_2, reduction = "umap", label = T, label.size = 5, pt.size = 1, cols = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E"))

all.markers.t.cells.2 <- FindAllMarkers(im_cosmx_t_cells_2)

all.markers.t.cells.2.reduced <- all.markers.t.cells.2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.1)

write.csv(all.markers.t.cells.2, ".../t_cell_subsets_dge.csv")



saveRDS(im_cosmx_t_cells_2, ".../base_T_cells.rds")

im_cosmx_t_cells_2 <- readRDS(".../base_T_cells.rds")
im_cosmx_t_cells_3 <- readRDS(".../base_T_cells_cd4_cd8.rds")
im_cosmx_t_cells_2 <- im_cosmx_t_cells_3

DimPlot(im_cosmx_t_cells_2, reduction = "umap", label = T, label.size = 5, pt.size = 2)


all.markers.t.cells.2.reduced.cluster <- all.markers.t.cells.2.reduced[all.markers.t.cells.2.reduced$cluster == 6,]

genes_per_cluster <- unique(c("H4C3","STMN1",
                       "CCL5","ITGAE","GZMA","TIMP1","NKG7",
                       "GNLY",
                       "LTB","MAF","IL7R","TCF7","SPOCK2",
                       "CCL4/L1/L2","CCL3/L1/L3","IFNG","IL10","LAG3",
                       "GZMB","TNFRSF9","HSP90B1","IL2RA","PRF1",
                       "XCL1/2","LIF","MYC"))

DotPlot(im_cosmx_t_cells_2, features = genes_per_cluster) + rotate_x_text(45) + 
  geom_point(mapping = aes_string(size = 'pct.exp'), colour = "black" , stroke = 1, shape = 21) +
  scale_color_viridis(option="magma") + coord_flip() + xlab("") + ylab("")

FeaturePlot(im_cosmx_t_cells_2, features = c("LAG3"), pt.size = 2)
FeaturePlot(im_cosmx_t_cells_2, features = c("PDCD1"), pt.size = 2)
FeaturePlot(im_cosmx_t_cells_2, features = c("HAVCR2"), pt.size = 2)

new.cluster.ids <- c("CD8","CD8","CD8","CD4","CD8","CD8","CD8")
names(new.cluster.ids) <- levels(im_cosmx_t_cells_2)
im_cosmx_t_cells_2 <- RenameIdents(im_cosmx_t_cells_2, new.cluster.ids)
DimPlot(im_cosmx_t_cells_2, reduction = "umap", label = TRUE, pt.size = 1, cols = c("black","red")) + NoLegend()

Idents(object = im_cosmx_t_cells_2) <- "seurat_clusters"

#######
# T cells - in niche or not
#######

t_cell_function_niche <- as.data.frame(cbind(as.character(im_cosmx_t_cells_2@meta.data$seurat_clusters),im_cosmx@meta.data$in_niche[im_cosmx@meta.data$cell_id %in% im_cosmx_t_cells_2@meta.data[["cell_id"]]]))
colnames(t_cell_function_niche) <- c("Cluster","Niche")

results <- c(table(t_cell_function_niche$Cluster[t_cell_function_niche$Niche == "LMP1+"]) / table(t_cell_function_niche$Cluster),
             table(t_cell_function_niche$Cluster[t_cell_function_niche$Niche == "EBV+"]) / table(t_cell_function_niche$Cluster),
             table(t_cell_function_niche$Cluster[t_cell_function_niche$Niche == "None"]) / table(t_cell_function_niche$Cluster),
             table(t_cell_function_niche$Cluster[t_cell_function_niche$Niche == "LMP1+ EBV+"]) / table(t_cell_function_niche$Cluster))

t_cells_in_niche <- as.data.frame(cbind(rep(c(0:6), times = 4),rep(c("LMP1+","EBV+","None","LMP1+ EBV+"), each = 7), results))
colnames(t_cells_in_niche) <- c("Cluster","Niche","Enrichment")
t_cells_in_niche$Enrichment <- as.numeric(t_cells_in_niche$Enrichment)

t_cells_in_niche$Enrichment <- t_cells_in_niche$Enrichment * 100
t_cells_in_niche$realtiveEnrichment <- t_cells_in_niche$Enrichment - 25



ggbarplot(t_cells_in_niche, x = "Cluster", y = "Enrichment", fill = "Niche")

ggbarplot(t_cells_in_niche, x = "Niche", y = "Enrichment", fill = "Cluster")


ggplot(t_cells_in_niche, aes(fill=Niche, y=realtiveEnrichment, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", colour = "black") + theme_classic() + 
  scale_fill_manual(values = c("None" = "#d62728", "LMP1+ EBV+" = "#2ca02c", "EBV+" = "#ff7f0e", "LMP1+" = "#1f77b4")) + 
  ylab("Relative Enrichment [%]") + xlab("") + 
  theme(axis.ticks.length=unit(.25, "cm"),
        text = element_text(size = 16))


#######
# DGE of T cell subsets
#######


all.markers.t.cells.cluster.5 <- FindMarkers(im_cosmx_t_cells_2, ident.1 = "5")

all.markers.t.cells.cluster.5.reduced <- all.markers.t.cells.cluster.5 %>%
  dplyr::filter(avg_log2FC > 0.1)

write.csv(all.markers.t.cells.cluster.5, ".../cluster_5_dge.csv")



#######
# T cell DEG between in niche or not
#######

im_cosmx_t_cells_2@meta.data$in_niche <- im_cosmx@meta.data$in_niche[im_cosmx@meta.data$cell_id %in% im_cosmx_t_cells_2@meta.data[["cell_id"]]]

im_cosmx_t_cells_3 <- subset(im_cosmx_t_cells_2, idents = c("0","1","2","4","5","6"))

im_cosmx_t_cells_3 <- SetIdent(im_cosmx_t_cells_3, value = "in_niche")

all.markers.t.cells.in.niche <- FindMarkers(im_cosmx_t_cells_3, ident.1 = "LMP1+")

all.markers.t.cells.in.niche.reduced <- all.markers.t.cells.in.niche %>%
  dplyr::filter(avg_log2FC > 0.05)

all.markers.t.cells.in.niche.reduced.lmp1 <- all.markers.t.cells.in.niche.reduced[all.markers.t.cells.in.niche.reduced$cluster == "LMP1+",]

write.csv(all.markers.t.cells.in.niche.reduced, ".../lmp1_t_cell_dge.csv")





library(EnhancedVolcano)

keyvals <- rep("darkgrey", times = nrow(all.markers.t.cells.in.niche))
names(keyvals) <- "All"

chemokines <- c("KLRB1","LEP","HSP90AB1","PFN1","HSP90AA1","ENO1","HLA-DRB","LIF")
mhc <- c("HLA-DRB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DQA1", "CD74")
targets <- c("HSP90AB1", "HSP90AA1", "HSP90B1", "HSPA1A/B", "HSPB1")

keyvals[rownames(lmp1.macs.niche) %in% chemokines] <- "red"
names(keyvals)[rownames(lmp1.macs.niche) %in% chemokines] <- "Chemokines (Receptors)"
keyvals[rownames(lmp1.macs.niche) %in% mhc] <- "blue"
names(keyvals)[rownames(lmp1.macs.niche) %in% mhc] <- "Antigen Presentation"
keyvals[rownames(lmp1.macs.niche) %in% targets] <- "green"
names(keyvals)[rownames(lmp1.macs.niche) %in% targets] <- "Oxidative stress"

#ebv_pos_vs_neg_epi[ebv_pos_vs_neg_epi$X == "CCL22", 'p_val_adj'] <- 0

EnhancedVolcano::EnhancedVolcano(all.markers.t.cells.in.niche,
                                 lab = rownames(all.markers.t.cells.in.niche),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 pCutoff = 10e-3,
                                 FCcutoff = 0.25,
                                 pointSize = 3.0,
                                 labSize = 6, 
                                 selectLab = chemokines, 
                                 colCustom = keyvals) + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, face = "plain"),
        axis.ticks.length=unit(.2, "cm"))





#######
# All cells DEG between in niche or not
#######

im_cosmx <- SetIdent(im_cosmx, value = "in_niche")

all.markers.in.niche <- FindAllMarkers(im_cosmx)

all.markers.in.niche.reduced <- all.markers.in.niche %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

all.markers.in.niche.reduced.lmp1 <- all.markers.in.niche.reduced[all.markers.in.niche.reduced$cluster == "LMP1+",]

#######
# T cells hypoxia markers
#######

library(escape)

GS.hallmark <- getGeneSets(library = "H")

glycoysis <- list(GS.hallmark[["HALLMARK_GLYCOLYSIS"]]@geneIds[GS.hallmark[["HALLMARK_GLYCOLYSIS"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])
hypoxia <- list(GS.hallmark[["HALLMARK_HYPOXIA"]]@geneIds[GS.hallmark[["HALLMARK_HYPOXIA"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])
ox_phor <- list(GS.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]@geneIds[GS.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])
fat_acid_met <- list(GS.hallmark[["HALLMARK_FATTY_ACID_METABOLISM"]]@geneIds[GS.hallmark[["HALLMARK_FATTY_ACID_METABOLISM"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])
angio <- list(GS.hallmark[["HALLMARK_ANGIOGENESIS"]]@geneIds[GS.hallmark[["HALLMARK_ANGIOGENESIS"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])
ros <- list(GS.hallmark[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]@geneIds[GS.hallmark[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])


im_cosmx_t_cells_2 <- SetIdent(im_cosmx_t_cells_2, value = "seurat_clusters")
im_cosmx_t_cells_2 <- SetIdent(im_cosmx_t_cells_2, value = "in_niche")

im_cosmx_t_cells_2 <- AddModuleScore(
  object = im_cosmx_t_cells_2,
  features = hypoxia,
  ctrl = 5,
  name = 'hypoxia'
)

VlnPlot(object = im_cosmx_t_cells_2, features = 'hypoxia1')
RidgePlot(object = im_cosmx_t_cells_2, features = 'hypoxia1')

im_cosmx_t_cells_2 <- AddModuleScore(
  object = im_cosmx_t_cells_2,
  features = glycoysis,
  ctrl = 5,
  name = 'glycoysis'
)

VlnPlot(object = im_cosmx_t_cells_2, features = 'glycoysis1', pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("Glycolysis Module Score")

RidgePlot(object = im_cosmx_t_cells_2, features = 'glycoysis1')

im_cosmx_t_cells_2 <- AddModuleScore(
  object = im_cosmx_t_cells_2,
  features = ox_phor,
  ctrl = 5,
  name = 'ox_phor'
)

VlnPlot(object = im_cosmx_t_cells_2, features = 'ox_phor1', pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("Ox Phos Module Score")

RidgePlot(object = im_cosmx_t_cells_2, features = 'ox_phor1')

im_cosmx_t_cells_2 <- AddModuleScore(
  object = im_cosmx_t_cells_2,
  features = fat_acid_met,
  ctrl = 5,
  name = 'fat_acid_met'
)

# 4.5 x 3 - landscape

VlnPlot(object = im_cosmx_t_cells_2, features = 'fat_acid_met1', pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("Fatty Acid Metabolism Module Score")



RidgePlot(object = im_cosmx_t_cells_2, features = 'fat_acid_met1')

inflam_response <- list(GS.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]]@geneIds[GS.hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]]@geneIds %in% rownames(im_cosmx_t_cells_2)])

im_cosmx_t_cells_2 <- AddModuleScore(
  object = im_cosmx_t_cells_2,
  features = inflam_response,
  ctrl = 5,
  name = 'inflam_response'
)

VlnPlot(object = im_cosmx_t_cells_2, features = 'inflam_response1', pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("Inflammatory Response Module Score")

RidgePlot(object = im_cosmx_t_cells_2, features = 'inflam_response1')




VlnPlot(object = im_cosmx_t_cells_2, features = c('ICOS'), pt.size = 0) + 
    scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("ICOS Expression")

VlnPlot(object = im_cosmx_t_cells_2, features = c('LAG3'), pt.size = 0)+ 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("LAG3 Expression")

VlnPlot(object = im_cosmx_t_cells_2, features = c('HAVCR2'), pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("HAVCR2 Expression")

VlnPlot(object = im_cosmx_t_cells_2, features = c('TIGIT'), pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("TIGIT Expression")

VlnPlot(object = im_cosmx_t_cells_2, features = c('PDCD1'), pt.size = 0) + 
  scale_fill_manual(values = c("#6CCFF6","#98CE00","#FFB30F","#9A031E","#7ADFBB","#BD9391","#211A1E")) + theme(legend.position = "none") + ggtitle("") + xlab("") + ylab("PDCD1 Expression")

#######
# GZMB+ volcano
#######

ido1.markers.macs.2 <- Seurat::FindMarkers(im_cosmx_, ident.1 = "6")

ido1.markers.macs.reduced.2 <- ido1.markers.macs.2 %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(ido1.markers.macs.2, ".../ido1_dge.csv")

cluster.5.t.cells <- read.csv(".../cluster_5_dge.csv", row.names = 1)

# plot macrophage dge

library(EnhancedVolcano)

keyvals <- rep("darkgrey", times = nrow(lmp1.macs.niche))
names(keyvals) <- "All"

chemokines <- c("CXCL9", "CXCL10", "CCL19", "CCL22", "CCR5","IL2RG","IL4R")
mhc <- c("HLA-DRB", "HLA-DPA1", "HLA-DPB1", "HLA-DRA", "HLA-DQA1", "CD74")
targets <- c("HSP90AB1", "HSP90AA1", "HSP90B1", "HSPA1A/B", "HSPB1")

keyvals[rownames(lmp1.macs.niche) %in% chemokines] <- "red"
names(keyvals)[rownames(lmp1.macs.niche) %in% chemokines] <- "Chemokines (Receptors)"
keyvals[rownames(lmp1.macs.niche) %in% mhc] <- "blue"
names(keyvals)[rownames(lmp1.macs.niche) %in% mhc] <- "Antigen Presentation"
keyvals[rownames(lmp1.macs.niche) %in% targets] <- "green"
names(keyvals)[rownames(lmp1.macs.niche) %in% targets] <- "Oxidative stress"

#ebv_pos_vs_neg_epi[ebv_pos_vs_neg_epi$X == "CCL22", 'p_val_adj'] <- 0

EnhancedVolcano::EnhancedVolcano(cluster.5.t.cells,
                                 lab = rownames(cluster.5.t.cells),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 pCutoff = 10e-3,
                                 FCcutoff = 0.25,
                                 pointSize = 3.0,
                                 labSize = 6) + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, face = "plain"),
        axis.ticks.length=unit(.2, "cm"))

#######
# XCL1/2+ volcano
#######

cluster.6.t.cells <- read.csv(".../cluster_6_dge.csv", row.names = 1)

# plot macrophage dge

library(EnhancedVolcano)

EnhancedVolcano::EnhancedVolcano(cluster.6.t.cells,
                                 lab = rownames(cluster.6.t.cells),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 pCutoff = 10e-3,
                                 FCcutoff = 0.25,
                                 pointSize = 3.0,
                                 labSize = 6) + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 15, face = "plain"),
        axis.text.y = element_text(color = "black", size = 15, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 15, face = "plain"),
        axis.title.y = element_text(color = "black", size = 15, face = "plain"),
        axis.ticks.length=unit(.2, "cm"))

#######
# GZMB+ 
#######

library(ggpubr)

genesets <- read.table(".../GO_cluster_6.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

genesets <- genesets[order(genesets$Adjusted.P.value),]

genesets$log10 <- -log10(genesets$Adjusted.P.value)

# genesets.reduced <- genesets[c(1,11,13),]
genesets.reduced <- genesets[c(1,3,4,5,6),]

genesets.reduced$Term <- as.factor(genesets.reduced$Term)
genesets.reduced$Term <- factor(genesets.reduced$Term, 
                                levels = rev(c("Cellular Response To Tumor Necrosis Factor (GO:0071356)",
                                               "Positive Regulation Of T Cell Migration (GO:2000406)",
                                               "Response To Interleukin-1 (GO:0070555)",
                                               "Apoptotic Process (GO:0006915)",
                                               "Cytokine-Mediated Signaling Pathway (GO:0019221)")))

ggbarplot(genesets.reduced, x = "Term", y = "log10", fill = "#68CCC9") + coord_flip() + xlab("") + ggdist::theme_ggdist() + theme(legend.position = "none",
                                                                                                                                  axis.text.x = element_text(size = 18, colour = "black"),
                                                                                                                                  axis.text.y = element_text(size = 12, colour = "black"),
                                                                                                                                  axis.title.x = element_text(size = 20, colour = "black"),
                                                                                                                                  axis.title.y = element_text(size = 20, colour = "black"),
                                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                  axis.ticks = element_line(color="black"))









#######
# CD8/CD4 ratio
#######

cd8_cd4_ratio <- data.frame(Ratio = c((table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 1"])[1] / table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 1"])[2]),
table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 2"])[1] / table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 2"])[2],
table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 3"])[1] / table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 3"])[2],
table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 4"])[1] / table(im_cosmx_t_cells_3@active.ident[im_cosmx_t_cells_3@meta.data$donor == "IM Tonsil 4"])[2]),
Group = rep("CD8/CD4 Ratio"))

ggdotplot(cd8_cd4_ratio, x = "Group", y = "Ratio", fill = "black", color = "black", size = 1.05) +
  xlab("")   + ylab("CD8/CD4 Ratio") +
  geom_hline(aes(yintercept = 1),
             linetype = 2) + ggdist::theme_ggdist() + theme(panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(),
                                                            legend.position = "none",
                                                            axis.text.x = element_text(size = 4, colour = "black"),
                                                            axis.text.y = element_text(size = 18, colour = "black"),
                                                            axis.title.y = element_text(size = 20, colour = "black"),
                                                            axis.title.x = element_text(size = 20, colour = "black"),
                                                            axis.ticks.length=unit(.25, "cm"),
                                                            axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                            axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                            axis.ticks = element_line(color="black"))  + rotate_x_text(60)





