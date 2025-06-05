library(Seurat)
library(RANN)
library(ggpubr)

im_cosmx <- readRDS(".../base_cleaned_IM_clustered.rds")

im_cosmx_ebv_b_cells_2 <- readRDS(".../base_B_cell.rds")

###########################
# Add B cell subsets back into the base cosmx seurat object
###########################

im_cosmx@meta.data$clustering.annotations.ebv.b.cells <- as.character(im_cosmx@meta.data$clustering.annotations)
im_cosmx@meta.data$clustering.annotations.ebv.b.cells[im_cosmx@meta.data$cell_id %in% im_cosmx_ebv_b_cells_2@meta.data[["cell_id"]]] <- as.character(im_cosmx_ebv_b_cells_2@meta.data$clustering.base.annotation)

table(im_cosmx@meta.data$clustering.annotations.ebv.b.cells)

###########################
# Import T cells and break into CD4s and CD8s, put them back into base cosmx seurat object
###########################

im_cosmx_t_cells_2 <- readRDS("...base_T_cells.rds")

DimPlot(im_cosmx_t_cells_2, reduction = "umap", label = T, label.size = 5, pt.size = 2)

DotPlot(im_cosmx_t_cells_2, features = c("CD3E","CD4","CD8A","CD8B","IL27RA",
                                         "FOXP3","NKG7","PFN1","GZMK","GZMA",
                                         "GZMH","PRF1","KLRB1","GNLY","TCF7",
                                         "IL7R","CCR7","IFNG","CD2")) + rotate_x_text(45) + xlab("")

im_cosmx_t_cells_2@meta.data$clustering.base.pre.cd4.cd8 <- im_cosmx_t_cells_2@active.ident
new_names <- c("CD8+ T cells","CD8+ T cells","CD8+ T cells","CD4+ T cells","CD8+ T cells","CD8+ T cells","CD8+ T cells")
names(new_names) <- levels(im_cosmx_t_cells_2)
im_cosmx_t_cells_2 <- RenameIdents(object = im_cosmx_t_cells_2, new_names)

saveRDS(im_cosmx_t_cells_2,".../base_T_cells_cd4_cd8.rds")

im_cosmx_t_cells_2 <- readRDS(".../base_T_cells_cd4_cd8.rds")

length(unname(Idents(im_cosmx_t_cells_2)))

im_cosmx@meta.data$clustering.annotations.ebv.b.cells <- as.character(im_cosmx@meta.data$clustering.annotations.ebv.b.cells)
im_cosmx@meta.data$clustering.annotations.ebv.b.cells[im_cosmx@meta.data$cell_id %in% im_cosmx_t_cells_2@meta.data[["cell_id"]]] <- as.character(unname(Idents(im_cosmx_t_cells_2)))

table(im_cosmx@meta.data$clustering.annotations.ebv.b.cells)

###########################
# Import macrophages and break into M1s and M2s, put them back into base cosmx seurat object
###########################

im_cosmx_macs_2 <- readRDS(".../base_macs.rds")

DimPlot(im_cosmx_macs_2, reduction = "umap", pt.size = 1, cols = c("#6295CB","#5CA53F","#F2903F","#EC5D6A","#C87DB4","#EBE747","#2FBFD8","#2BDD88"))

mac.features <- c("CCL18","IL2RA","GPNMB","CXCL12","C1QB","CCL5","MMP12","MMP9","SPP1",
                  "S100A8","S100A9","IL1B","DUSP1","CXCL10","CXCL9","CD44","MT2A","CCL8","MT1X","CCL2","ITGAX",
                  "IDO1","CD74","CIITA","HLA-DQB1/2","HSPA1A/B","THBS1","IGFBP5","COL3A1","TIMP1")

DotPlot(im_cosmx_macs_2, features = mac.features) + rotate_x_text(45) + xlab("")


mac.features <- c("CD14","CD68","CD163","STAT1","TLR4","CD86","CD80","MRC1","CD209",
                  "STAT3","STAT6","IL10","IL6","TLR8","ARG1","IL6R",
                  "CXCL9","CXCL10","CCL5","CXCL11","CXCL16","CIITA","HLA-DRA","HLA-DRB")

DotPlot(im_cosmx_macs_2, features = mac.features) + rotate_x_text(45) + xlab("")


im_cosmx_macs_2@meta.data$clustering.base.pre.m1.m2 <- im_cosmx_macs_2@active.ident
new_names <- c("M0","M1","M2","M1","M2","M2","M1","M1")
names(new_names) <- levels(im_cosmx_macs_2)
im_cosmx_macs_2 <- RenameIdents(object = im_cosmx_macs_2, new_names)

saveRDS(im_cosmx_macs_2,".../base_macs_m1_m2.rds")

im_cosmx@meta.data$clustering.annotations.ebv.b.cells <- as.character(im_cosmx@meta.data$clustering.annotations.ebv.b.cells)
im_cosmx@meta.data$clustering.annotations.ebv.b.cells[im_cosmx@meta.data$cell_id %in% im_cosmx_macs_2@meta.data[["cell_id"]]] <- as.character(unname(Idents(im_cosmx_macs_2)))

saveRDS(im_cosmx, ".../base_cleaned_IM_clustered_with_cd4_cd8_m1_m2_ebv_b_cell.rds")

im_cosmx <-readRDS(".../base_cleaned_IM_clustered_with_cd4_cd8_m1_m2_ebv_b_cell.rds")

###########################
# Calculate which cells are in EBV+ / LMP1+ B cell niches
###########################

Idents(im_cosmx) <- im_cosmx@meta.data$clustering.annotations.ebv.b.cells

im_cosmx@meta.data$in_niche <- rep("None", times = nrow(im_cosmx@meta.data))
im_cosmx@meta.data$in_ebv_niche <- rep("None", times = nrow(im_cosmx@meta.data))
im_cosmx@meta.data$in_lmp1_niche <- rep("None", times = nrow(im_cosmx@meta.data))

# LMP1+

dists <- RANN::nn2(im_cosmx@meta.data[,c("x_slide_mm","y_slide_mm")],
                   query = im_cosmx@meta.data[im_cosmx@meta.data$clustering.annotations.ebv.b.cells == "LMP1+ EBV+ B cells",c("x_slide_mm","y_slide_mm")], 
                   radius = 0.1, k = 200)

cutoff <- 0.05 # in mm, so 0.05 is 50um

for(i in 1:nrow(dists[["nn.dists"]])){
  im_cosmx@meta.data$in_lmp1_niche[dists[["nn.idx"]][i,dists[["nn.dists"]][i,] < cutoff]] <- "LMP1+"
}

# EBV+

dists <- RANN::nn2(im_cosmx@meta.data[,c("x_slide_mm","y_slide_mm")],
                   query = im_cosmx@meta.data[im_cosmx@meta.data$clustering.annotations.ebv.b.cells == "EBV+ B cells",c("x_slide_mm","y_slide_mm")], 
                   radius = 0.1, k = 200)

for(i in 1:nrow(dists[["nn.dists"]])){
  im_cosmx@meta.data$in_ebv_niche[dists[["nn.idx"]][i,dists[["nn.dists"]][i,] < cutoff]] <- "EBV+"
}


# Combine niches

im_cosmx@meta.data$in_niche <- paste0(im_cosmx@meta.data$in_lmp1_niche," ", im_cosmx@meta.data$in_ebv_niche)
im_cosmx@meta.data$in_niche[im_cosmx@meta.data$in_niche == "None None"] <- "None"
im_cosmx@meta.data$in_niche[im_cosmx@meta.data$in_niche == "LMP1+ None"] <- "LMP1+"
im_cosmx@meta.data$in_niche[im_cosmx@meta.data$in_niche == "None EBV+"] <- "EBV+"

# Save RDS file

saveRDS(im_cosmx, ".../base_cleaned_IM_clustered_with_cd4_cd8_m1_m2_ebv_b_cell_with_niches.rds")

im_cosmx <- readRDS(".../base_cleaned_IM_clustered_with_cd4_cd8_m1_m2_ebv_b_cell_with_niches.rds")

###########################
# Split into subsetted macrophage and T cell subsets for plotting
###########################

# I saved the pdfs from this as 4.5 x 4 in landscape

chemokines <- c("SPP1")

median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

# M1 

im_cosmx_m1 <- subset(im_cosmx, idents = c("M1"))

im_cosmx_m1@meta.data$in_niche <- as.factor(im_cosmx_m1@meta.data$in_niche)
im_cosmx_m1@meta.data$in_niche <- factor(im_cosmx_m1@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

#im_cosmx_m1 <- ScaleData(im_cosmx_m1, assay = "SCT")

Seurat::VlnPlot(im_cosmx_m1, features = "SPP1", split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_m1) <- im_cosmx_m1@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_m1, features = "CD274", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T) + 
  stat_compare_means(comparisons = list(c("LMP1+ EBV+","LMP1+"), c("EBV+","LMP1+"),c("None","LMP1+"))) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) + theme(legend.position = "none") + xlab("")

table(im_cosmx_m1@meta.data$in_niche)

all.markers.m1.niches <- FindAllMarkers(im_cosmx_m1)

all.markers.m1.niches.reduced <- all.markers.m1.niches %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.m1.niches.reduced, ".../m1.csv")

# M2 


im_cosmx_m2 <- subset(im_cosmx, idents = c("M2"))

im_cosmx_m2@meta.data$in_niche <- as.factor(im_cosmx_m2@meta.data$in_niche)
im_cosmx_m2@meta.data$in_niche <- factor(im_cosmx_m2@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

#im_cosmx_m2 <- ScaleData(im_cosmx_m2, assay = "SCT")

im_cosmx_m2@active.assay <- "RNA"
im_cosmx_m2@active.assay <- "SCT"

Seurat::VlnPlot(im_cosmx_m2, features = "SPP1", split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_m2) <- im_cosmx_m2@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_m2, features = "CD274", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T)  + 
  stat_compare_means(comparisons = list(c("LMP1+ EBV+","LMP1+"),c("EBV+","LMP1+"),c("None","LMP1+"))) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) + theme(legend.position = "none") + xlab("")

table(im_cosmx_m2@meta.data$in_niche)


all.markers.m2.niches <- FindAllMarkers(im_cosmx_m2)

all.markers.m2.niches.reduced <- all.markers.m2.niches %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.m2.niches.reduced, ".../m2.csv")

# CD4

im_cosmx_cd4 <- subset(im_cosmx, idents = c("CD4+ T cells"))

im_cosmx_cd4@meta.data$in_niche <- as.factor(im_cosmx_cd4@meta.data$in_niche)
im_cosmx_cd4@meta.data$in_niche <- factor(im_cosmx_cd4@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

Seurat::VlnPlot(im_cosmx_cd4, features = chemokines[3], split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_cd4) <- im_cosmx_cd4@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_cd4, features = "CXCR4", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T)  + 
  stat_compare_means(comparisons = list(c("LMP1+ EBV+","LMP1+"),c("EBV+","LMP1+"),c("None","LMP1+"))) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) + theme(legend.position = "none") + xlab("")

table(im_cosmx_cd4@meta.data$in_niche)

all.markers.cd4.niches <- FindAllMarkers(im_cosmx_cd4)

all.markers.cd4.niches.reduced <- all.markers.cd4.niches %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.cd4.niches.reduced, ".../cd4.csv")

# CD8

im_cosmx_cd8 <- subset(im_cosmx, idents = c("CD8+ T cells"))

im_cosmx_cd8@meta.data$in_niche <- as.factor(im_cosmx_cd8@meta.data$in_niche)
im_cosmx_cd8@meta.data$in_niche <- factor(im_cosmx_cd8@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

Seurat::VlnPlot(im_cosmx_cd8, features = chemokines[4], split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_cd8) <- im_cosmx_cd8@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_cd8, features = "CXCR4", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T)  + 
  stat_compare_means(comparisons = list(c("LMP1+ EBV+","LMP1+"),c("EBV+","LMP1+"),c("None","LMP1+")))  +
                       stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95) + theme(legend.position = "none") + xlab("")

table(im_cosmx_cd8@meta.data$in_niche)

all.markers.cd8.niches <- FindAllMarkers(im_cosmx_cd8)

all.markers.cd8.niches.reduced <- all.markers.cd8.niches %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.05)

write.csv(all.markers.cd8.niches.reduced, ".../cd8.csv")


# M0

im_cosmx_cxcl8 <- subset(im_cosmx, idents = c("M0"))

im_cosmx_cxcl8@meta.data$in_niche <- as.factor(im_cosmx_cxcl8@meta.data$in_niche)
im_cosmx_cxcl8@meta.data$in_niche <- factor(im_cosmx_cxcl8@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

#im_cosmx_m2 <- ScaleData(im_cosmx_m2, assay = "SCT")

Seurat::VlnPlot(im_cosmx_cxcl8, features = chemokines[4], split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_cxcl8) <- im_cosmx_cxcl8@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_cxcl8, features = "CX3CR1", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T)  + 
  stat_compare_means(comparisons = list(c("EBV+","LMP1+"),c("None","LMP1+"))) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95)


# Fibroblasts

im_cosmx_fibroblasts <- subset(im_cosmx, idents = c("Endothelial cells"))

im_cosmx_fibroblasts@meta.data$in_niche <- as.factor(im_cosmx_fibroblasts@meta.data$in_niche)
im_cosmx_fibroblasts@meta.data$in_niche <- factor(im_cosmx_fibroblasts@meta.data$in_niche, levels = c("None","EBV+","LMP1+ EBV+","LMP1+"))

#im_cosmx_m2 <- ScaleData(im_cosmx_m2, assay = "SCT")

Seurat::VlnPlot(im_cosmx_fibroblasts, features = chemokines[4], split.by = 'in_niche', slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"))

Idents(im_cosmx_fibroblasts) <- im_cosmx_fibroblasts@meta.data$in_niche

Seurat::VlnPlot(im_cosmx_fibroblasts, features = "CXCL12", slot = "scale.data", cols = c("#364CA2","#4CC7EF","#AD3496","#18A58A"), y.max = 10, log = T)  + 
  stat_compare_means(comparisons = list(c("EBV+","LMP1+"),c("None","LMP1+"))) +
  stat_summary(fun.y = median, geom='point', size = 15, colour = "white", shape = 95)









