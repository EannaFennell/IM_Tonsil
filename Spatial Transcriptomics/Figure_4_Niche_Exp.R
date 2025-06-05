library(Seurat)
library(ggpubr)
library(RANN)

im_cosmx <- readRDS(".../base_cleaned_IM_clustered_with_cd4_cd8_m1_m2_ebv_b_cell_with_niches.rds")

im_cosmx_macs_2 <- readRDS(".../base_macs.rds")

###########################
# Add B cell subsets back into the base cosmx seurat object
###########################

im_cosmx_macs_2@meta.data$macrophage_clusters <- im_cosmx_macs_2@active.ident

mac_function_niche <- as.data.frame(cbind(as.character(im_cosmx_macs_2@meta.data$macrophage_clusters),im_cosmx@meta.data$in_niche[im_cosmx@meta.data$cell_id %in% im_cosmx_macs_2@meta.data[["cell_id"]]]))
colnames(mac_function_niche) <- c("Cluster","Niche")

results <- c(table(mac_function_niche$Cluster[mac_function_niche$Niche == "LMP1+"]) / table(mac_function_niche$Cluster),
table(mac_function_niche$Cluster[mac_function_niche$Niche == "EBV+"]) / table(mac_function_niche$Cluster),
table(mac_function_niche$Cluster[mac_function_niche$Niche == "None"]) / table(mac_function_niche$Cluster),
table(mac_function_niche$Cluster[mac_function_niche$Niche == "LMP1+ EBV+"]) / table(mac_function_niche$Cluster))

macs_in_niche <- as.data.frame(cbind(rep(c(0:7), times = 4),rep(c("LMP1+","EBV+","None","LMP1+ EBV+"), each = 8), results))
colnames(macs_in_niche) <- c("Cluster","Niche","Enrichment")
macs_in_niche$Enrichment <- as.numeric(macs_in_niche$Enrichment)


macs_in_niche$Niche <- factor(macs_in_niche$Niche, levels = c("None", "LMP1+ EBV+", "EBV+", "LMP1+"))

p <- ggbarplot(macs_in_niche, x = "Cluster", y = "Enrichment", fill = "Niche")

p + scale_fill_manual(values = c("None" = "#d62728", "LMP1+ EBV+" = "#2ca02c", "EBV+" = "#ff7f0e", "LMP1+" = "#1f77b4"))


macs_in_niche$Enrichment <- macs_in_niche$Enrichment * 100

macs_in_niche$realtiveEnrichment <- macs_in_niche$Enrichment - 25


ggbarplot(macs_in_niche, x = "Cluster", y = "Enrichment", fill = "Niche")

ggbarplot(macs_in_niche, x = "Niche", y = "Enrichment", fill = "Cluster")



ggplot(macs_in_niche, aes(fill=Niche, y=realtiveEnrichment, x=Cluster)) + 
  geom_bar(position="dodge", stat="identity", colour = "black") + theme_classic() + 
  scale_fill_manual(values = c("None" = "#d62728", "LMP1+ EBV+" = "#2ca02c", "EBV+" = "#ff7f0e", "LMP1+" = "#1f77b4")) + 
  ylab("Relative Enrichment [%]") + xlab("") + 
  theme(axis.ticks.length=unit(.25, "cm"),
        text = element_text(size = 16))


###########################
# DGE macrophages
###########################

all.markers.macs.2 <- FindAllMarkers(im_cosmx_macs_2)

all.markers.macs.reduced.2 <- all.markers.macs.2 %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.macs.reduced.2, ".../macrophage_subclusters_dge.csv")

###########################
# DGE macrophages - CXCL10
###########################

cxcl10.markers.macs.2 <- Seurat::FindMarkers(im_cosmx_macs_2, ident.1 = "4")

cxcl10.markers.macs.reduced.2 <- cxcl10.markers.macs.2 %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(cxcl10.markers.macs.2, ".../cxcl10_dge.csv")


###########################
# DGE macrophages - SPP1
###########################

spp1.markers.macs.2 <- Seurat::FindMarkers(im_cosmx_macs_2, ident.1 = "2")

spp1.markers.macs.reduced.2 <- spp1.markers.macs.2 %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(spp1.markers.macs.2, ".../spp1_dge.csv")

###########################
# DGE macrophages - IDO1
###########################

ido1.markers.macs.2 <- Seurat::FindMarkers(im_cosmx_macs_2, ident.1 = "6")

ido1.markers.macs.reduced.2 <- ido1.markers.macs.2 %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(ido1.markers.macs.2, ".../ido1_dge.csv")

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

EnhancedVolcano::EnhancedVolcano(cxcl10.markers.macs.2,
                                 lab = rownames(cxcl10.markers.macs.2),
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


###########################
# DGE macrophages - in niche etc
###########################

im_cosmx_macs_2@meta.data$in_niche <- im_cosmx@meta.data$in_niche[im_cosmx@meta.data$cell_id %in% im_cosmx_macs_2@meta.data[["cell_id"]]]

im_cosmx_macs_2 <- SetIdent(im_cosmx_macs_2, value = "in_niche")

all.markers.macs.niche <- FindAllMarkers(im_cosmx_macs_2)

all.markers.macs.niche.reduced <- all.markers.macs.niche %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.macs.niche.reduced, ".../macrophage_niche_dge.csv")


lmp1.macs.niche <- FindMarkers(im_cosmx_macs_2, ident.1 = "LMP1+")


write.csv(lmp1.macs.niche, ".../macrophage_niche_dge_lmp1_vs_others.csv")


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

EnhancedVolcano::EnhancedVolcano(lmp1.macs.niche,
                                 lab = rownames(lmp1.macs.niche),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj',
                                 pCutoff = 10e-3,
                                 FCcutoff = 0.25,
                                 pointSize = 3.0,
                                 labSize = 4,
                                 colCustom = keyvals,
                                 selectLab = c(chemokines,mhc,targets)) + theme_classic() + 
  theme(axis.text.x = element_text(color = "black", size = 12, face = "plain"),
        axis.text.y = element_text(color = "black", size = 12, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 12, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12, face = "plain"),
        axis.ticks.length=unit(.2, "cm"))

# selectLab = c(chemokines,mhc,targets,"CRIP1","CCL5")

# Pathway

library(ggpubr)

genesets <- read.table(".../lmp1_niche_pathway.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

genesets <- genesets[order(genesets$Adjusted.P.value),]

genesets$log10 <- -log10(genesets$Adjusted.P.value)


genesets.reduced <- genesets[c(1,3,4,6,7,8),]

genesets.reduced$Term <- as.factor(genesets.reduced$Term)
genesets.reduced$Term <- factor(genesets.reduced$Term, 
                                levels = rev(c("Cellular Response To Cytokine Stimulus (GO:0071345)",
                                               "Positive Regulation Of Cytokine Production (GO:0001819)",
                                               "Inflammatory Response (GO:0006954)",
                                               "Positive Regulation Of Protein Phosphorylation (GO:0001934)",
                                               "Positive Regulation Of Macromolecule Metabolic Process (GO:0010604)",
                                               "Interleukin-4-Mediated Signaling Pathway (GO:0035771)")))



ggbarplot(genesets.reduced, x = "Term", y = "log10", fill = "#68CCC9") + coord_flip() + xlab("") + ggdist::theme_ggdist() + theme(legend.position = "none",
                                                                                                                                  axis.text.x = element_text(size = 18, colour = "black"),
                                                                                                                                  axis.text.y = element_text(size = 12, colour = "black"),
                                                                                                                                  axis.title.x = element_text(size = 20, colour = "black"),
                                                                                                                                  axis.title.y = element_text(size = 20, colour = "black"),
                                                                                                                                  axis.ticks.length=unit(.25, "cm"),
                                                                                                                                  axis.line.x = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                  axis.line.y = element_line(color = "black", linewidth = rel(0.5)),
                                                                                                                                  axis.ticks = element_line(color="black"))

###########################
# GSEA macrophages - in niche etc
###########################

library(escape)

GS.hallmark <- getGeneSets(library = "H")

glycoysis <- list(GS.hallmark[["HALLMARK_GLYCOLYSIS"]]@geneIds[GS.hallmark[["HALLMARK_GLYCOLYSIS"]]@geneIds %in% rownames(im_cosmx)])
hypoxia <- list(GS.hallmark[["HALLMARK_HYPOXIA"]]@geneIds[GS.hallmark[["HALLMARK_HYPOXIA"]]@geneIds %in% rownames(im_cosmx)])
ox_phor <- list(GS.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]@geneIds[GS.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]]@geneIds %in% rownames(im_cosmx)])
fat_acid_met <- list(GS.hallmark[["HALLMARK_FATTY_ACID_METABOLISM"]]@geneIds[GS.hallmark[["HALLMARK_FATTY_ACID_METABOLISM"]]@geneIds %in% rownames(im_cosmx)])
angio <- list(GS.hallmark[["HALLMARK_ANGIOGENESIS"]]@geneIds[GS.hallmark[["HALLMARK_ANGIOGENESIS"]]@geneIds %in% rownames(im_cosmx)])
ros <- list(GS.hallmark[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]@geneIds[GS.hallmark[["HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]]@geneIds %in% rownames(im_cosmx)])


im_cosmx_macs_2@meta.data$in_niche_clusters <- paste0(im_cosmx_macs_2@meta.data$seurat_clusters, " ", im_cosmx_macs_2@meta.data$in_niche)

im_cosmx_macs_2 <- SetIdent(im_cosmx_macs_2, value = "in_niche")
im_cosmx_macs_2 <- SetIdent(im_cosmx_macs_2, value = "seurat_clusters")
im_cosmx_macs_2 <- SetIdent(im_cosmx_macs_2, value = "in_niche_clusters")

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = glycoysis,
  ctrl = 5,
  name = 'glycoysis'
)

RidgePlot(object = im_cosmx_macs_2, features = 'glycoysis1')
VlnPlot(object = im_cosmx_macs_2, features = 'glycoysis1')

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = hypoxia,
  ctrl = 5,
  name = 'hypoxia'
)

VlnPlot(object = im_cosmx_macs_2, features = 'hypoxia1')

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = ox_phor,
  ctrl = 5,
  name = 'ox_phor'
)

VlnPlot(object = im_cosmx_macs_2, features = 'ox_phor1')
RidgePlot(object = im_cosmx_macs_2, features = 'ox_phor1')

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = fat_acid_met,
  ctrl = 5,
  name = 'fat_acid_met'
)

VlnPlot(object = im_cosmx_macs_2, features = 'fat_acid_met1')
RidgePlot(object = im_cosmx_macs_2, features = 'fat_acid_met1')

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = angio,
  ctrl = 5,
  name = 'angio'
)

RidgePlot(object = im_cosmx_macs_2, features = 'angio1')

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = ros,
  ctrl = 5,
  name = 'ros'
)

VlnPlot(object = im_cosmx_macs_2, features = 'ros1')
RidgePlot(object = im_cosmx_macs_2, features = 'ros1')


VlnPlot(object = im_cosmx_macs_2, features = 'HIF1A')


###########################
# Macrophages expression
###########################

im_cosmx_macs_2 <- AddModuleScore(
  object = im_cosmx_macs_2,
  features = list(c("CIITA","HLA-DRB","HLA-DPA1","HLA-DRA","HLA-DQB1/2","HLA-DPB1")),
  ctrl = 5,
  name = 'Class_2_score'
)

im_cosmx_macs_2@meta.data$in_niche <- as.character(im_cosmx@meta.data$in_niche[im_cosmx@meta.data$cell_id %in% im_cosmx_macs_2@meta.data[["cell_id"]]])


# MHC Class I
VlnPlot(im_cosmx_macs_2, features = c("TAP1","TAP2","MHC I"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8"), log = TRUE)

# MHC Class II
VlnPlot(im_cosmx_macs_2, features = c("CIITA","HLA-DRB","HLA-DPA1","HLA-DRA","HLA-DQB1/2","HLA-DPB1"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8"), log = TRUE)

# Proliferative
VlnPlot(im_cosmx_macs_2, features = c("CCND1","RGCC","MKI67"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8"), log = TRUE)

# Other
VlnPlot(im_cosmx_macs_2, features = c("CSF2RA","CSF1R","CD274","GAS6","IDO1"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8"), log = TRUE)

# Suppressive
VlnPlot(im_cosmx_macs_2, features = c("CHI3L1","GPNMB","MMP9"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8"), log = TRUE)


# Class 2 score
VlnPlot(im_cosmx_macs_2, features = c("Class_2_score1"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8","black"), log = TRUE, split.by = "in_niche") + xlab("") + ylab("MHC Class II score") + 
  theme(legend.position = "none") + labs(title = "")


VlnPlot(im_cosmx_macs_2, features = c("TAP2"), layer = "data", pt.size = 0, 
        idents = c("2","4","6"), cols = c("#F2903F","#C87DB4","#2FBFD8","black"), log = TRUE, split.by = "in_niche") + xlab("") + ylab("MHC Class II score") + labs(title = "")





###########################
# Pathway analysis of clusters
###########################

cluster3 <- read.table(".../Cluster_3_Pathway.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

cluster4 <- read.table(".../Cluster_4_Pathway.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

cluster5 <- read.table(".../Cluster_5_Pathway.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

cluster7 <- read.table(".../Cluster_7_Pathway.txt", header = TRUE,
                       fill = TRUE, sep = "\t")

pathway_heatmap <- as.data.frame(phonTools::zeros(4,25)) + 1
colnames(pathway_heatmap) <- unique(c(cluster3[c(1,2,3,6,7,10,11),"Term"],cluster4[c(1,2,3,4,6,7,12),"Term"],cluster5[c(1,2,3,5,6,7,11),"Term"],cluster7[c(1,3,5,7,9,12,14),"Term"]))
rownames(pathway_heatmap) <- c("SPP1+","IL1B+","CXCL10+","IDO1+")

clus.3.inter <- cluster3[cluster3$Term %in% colnames(pathway_heatmap), c(1,3)]
pathway_heatmap[1,clus.3.inter$Term] <- clus.3.inter$P.value

clus.4.inter <- cluster4[cluster4$Term %in% colnames(pathway_heatmap), c(1,3)]
pathway_heatmap[2,clus.4.inter$Term] <- clus.4.inter$P.value

clus.5.inter <- cluster5[cluster5$Term %in% colnames(pathway_heatmap), c(1,3)]
pathway_heatmap[3,clus.5.inter$Term] <- clus.5.inter$P.value

clus.7.inter <- cluster7[cluster7$Term %in% colnames(pathway_heatmap), c(1,3)]
pathway_heatmap[4,clus.7.inter$Term] <- clus.7.inter$P.value

pathway_heatmap <- -log10(pathway_heatmap)

#pathway_heatmap_plot <- log2(pathway_heatmap)
#pathway_heatmap_plot[pathway_heatmap_plot == -Inf] <- 0

library(circlize)
library(ComplexHeatmap)

pathway_heatmap[pathway_heatmap > 10] <- 10

ComplexHeatmap::Heatmap(as.matrix(pathway_heatmap), col = viridis::viridis(100), rect_gp = gpar(col = "grey50", lwd = 2), show_column_dend = FALSE, 
                        show_column_names = TRUE, heatmap_legend_param = list(title = "-log10(P value)", 
                                                                               title_gp = gpar(fontsize = 11, fontface = "bold"), 
                                                                               legend_height = unit(6, "cm"), title_position = "leftcenter-rot"))


DimPlot(im_cosmx_macs_2, reduction = "umap", pt.size = 1, cols = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88"))

FeaturePlot(im_cosmx_macs_2, feature = c("CD274"), pt.size = 1)


