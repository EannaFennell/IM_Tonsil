library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(MISSILe)
library(ggpubr)
library(FastPG)
library(grid)
library(ComplexHeatmap)
library(ggsignif)
library(RColorBrewer)

############################### 
# Import Latest MISSILe Object
############################### 

saveRDS(IM_TMA, ".../TMA_Select_Final.rds")
IM_TMA <- readRDS(".../IDO1_ME.rds")

############################### 
# Gene Lists
###############################
all.genes.TMA <- colnames(IM_TMA@Expression[["MultiIHC"]]@counts)[c(2:52)]
lineageMarkers <- colnames(IM_TMA@Expression$MultiIHC@counts)[c(3,4,6,9,10,12,13,14,16,17,21,22,23,26,34,39,42,44)]

############################### 
# Clustering
###############################
IM_TMA@current.identity$ActiveIdents <- "Final_Annotations"
table(IM_TMA@meta.data$Final_Annotations)

MISSILe::markerCorrelation(MISSILeObject = IM_TMA, expSet = "counts", markers = all.genes.TMA)

# Make T cell cluster from recluster
length(levels(IM_TMA@meta.data$Final_Annotations))
levels(IM_TMA@meta.data$Final_Annotations)
T_annotation <- c("B cells","CD4+ T cells","CD8+ T cells","Dendritic cells","Endothelium","Epithelium","GC B cells","LMP1+ B cells",
       "M1 macrophages","M2 macrophages", "Myeloid","Neutrophils","Plasma cells","T cells","Smooth muscle","T regs","Undefined")
length(T_annotation)
table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
levels(IM_TMA@meta.data$Final_Annotations) <- T_annotation
table(IM_TMA@meta.data$Final_Annotations)

MISSILe::BubblePlot(MISSILeObject = IM_TMA, markers = lineageMarkers,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Final_Annotations",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5"))+ grids(linetype = "dashed")

MISSILe::plotClusterIntensities(MISSILeObject = IM_TMA, markers = lineageMarkers,
                                clusteringName = "Final_Annotations", enrichmentLimits = c(-1.5,1.5)) 

plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "B cells", fill = rainbow(55), yscale = c(0,120))
quantile(x = IM_TMA@Expression$MultiIHC@counts$LMP1[IM_TMA@meta.data$Final_Annotations == "B cells"], probs = c(0,0.1,0.5,0.7,0.8,0.9,0.95))

#########################
# Clustering to match FF
#########################

table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ Macrophages"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "M2 macrophages"] <- "CD68+ CD163+ Macrophages"
table(IM_TMA@meta.data$Final_Annotations)

saveRDS(IM_TMA, ".../TMA_Select_Final_After_Match.rds")

############################### 
# tSNE, Downsample function from tonsilAnalysis.R
############################### 

colnames(IM_TMA@Expression$MultiIHC@counts)
IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
levels(IM_TMA@meta.data$Final_Annotations)
orderMarkers_IM_TMA <- list(c(39),                       #Bcells
                            c(4,34,13),                  #CD4+Tcells
                            c(10,6,26),                  #M2s
                            c(10,26,44,13),              #M1s
                            c(17,34),                    #CD8+Tcells
                            c(10,16,13,44),              #DCs
                            c(3),                        #Endothelium
                            c(12),                       #Epithelium
                            c(39,23),                    #GC B
                            c(39,42, 13),                #LMP1+ B cells
                            c(26),                       #Myeloid
                            c(44),                       #Neutrophils
                            c(14),                       #Plasma
                            c(9),                        #Smooth muscle
                            c(34),                       #T cells
                            c(4,34,21),                  #T regs
                            c(1))                        #Unknown

a <- downSampleVisualisation(clusters = IM_TMA@meta.data$Final_Annotations,
                             downsamplePercent = 0.2,
                             countsMatrix = IM_TMA@Expression$MultiIHC@counts,
                             markers = orderMarkers_IM_TMA)

a <- a[!(a$membership %in% c("Unknown")),]

reducedExp <- IM_TMA@Expression$MultiIHC@counts[as.numeric(a$index),lineageMarkers]

rtsne_out <- Rtsne::Rtsne(reducedExp, verbose = TRUE, perplexity = 30, num_threads = 14)

rtsne_plot <- data.frame(x = rtsne_out$Y[,1],
                         y = rtsne_out$Y[,2],
                         Clusters = a$membership)
rtsne_plot$Clusters <- as.factor(rtsne_plot$Clusters)

colorIn <- rtsne_plot$Clusters

centroids.data <- zeros(length(levels(colorIn)),2)
colnames(centroids.data) <- c("centroids_x","centroids_y")
levelsTo <- levels(colorIn)
for(i in 1:length(unique(colorIn))){
  x_coords <- rtsne_plot$x[colorIn == levelsTo[i]]
  y_coords <- rtsne_plot$y[colorIn == levelsTo[i]]
  centroids.data[i,1] <- mean(x_coords)
  centroids.data[i,2] <- mean(y_coords)
  rm(x_coords, y_coords)
}

qualpalette <- qualpalr::qualpal(17, "pretty")   # pretty, pretty_dark, rainbow, pastels

ggplot() + geom_point(data = rtsne_plot, aes(x=x, y=y, color=Clusters), alpha = 0.7) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_bw() + labs(y= "tSNE 2", x = "tSNE 1") +
  theme_classic() + scale_color_manual(values = qualpalette$hex) #+
geom_text(data = as.data.frame(centroids.data), mapping = aes(x = centroids_x,
                                                              y = centroids_y,
                                                              label = 1:length(unique(colorIn))),
          color = "black", size = 6, fontface = 2)


############################################
### Phenotype of LMP1 an B cell clusters ###
#############################################

mean(IM_TMA@Expression$MultiIHC@counts$Ki67[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ki67[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ki67[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ki67[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD44[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD44[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD44[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD44[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$PDL1[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$PDL1[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$PDL1[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$PDL1[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$IDO1[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IDO1[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IDO1[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IDO1[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$ICOS[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$ICOS[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$ICOS[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$ICOS[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD30[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD30[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD30[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD30[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$TIGIT[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIGIT[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIGIT[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIGIT[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CTLA4[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CTLA4[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CTLA4[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CTLA4[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$Ox40[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ox40[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ox40[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Ox40[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$IFNG[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IFNG[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IFNG[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$IFNG[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD69[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD69[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD69[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD69[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD40[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$T.bet[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$T.bet[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$T.bet[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$T.bet[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD40L[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40L[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40L[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD40L[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$VISTA[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$VISTA[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$VISTA[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$VISTA[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$DDR1[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$DDR1[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$DDR1[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$DDR1[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$TIM3[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIM3[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIM3[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TIM3[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$TCF7[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TCF7[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TCF7[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TCF7[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$Col6[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Col6[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Col6[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$Col6[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$LAG3[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$LAG3[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$LAG3[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$LAG3[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$BCL2[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$BCL2[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$BCL2[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$BCL2[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD39[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD39[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD39[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD39[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$TP63[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TP63[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TP63[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$TP63[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD20[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD20[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD20[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD20[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

mean(IM_TMA@Expression$MultiIHC@counts$CD38[IM_TMA@meta.data$Final_Annotations == "B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD38[IM_TMA@meta.data$Final_Annotations == "GC B cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD38[IM_TMA@meta.data$Final_Annotations == "Plasma cells"])
mean(IM_TMA@Expression$MultiIHC@counts$CD38[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"])

#####################
# Expression Boxplots
#####################

cell_types_of_interest <- c("LMP1+ B cells", "Plasma cells", "GC B cells", "B cells")

cell_type_order <- c( "B cells", "GC B cells", "Plasma cells", "LMP1+ B cells")

cell_type_colors <- c(
  "B cells" = "black",
  "GC B cells" = "#606060",
  "Plasma cells" = "#CCCCCC",
  "LMP1+ B cells" = "#CF2626")

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)

#ki67
#%Pos

table(IM_TMA@meta.data$Ki67_Positive)
quantile(IM_TMA@Expression$MultiIHC@counts$Ki67[IM_TMA@meta.data$allRegions == "1"], probs = c(0.1,0.6,0.7,0.75,0.8,0.817,0.85,0.905,0.98,0.984,0.99)) 
IM_TMA@meta.data$Ki67_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
IM_TMA@meta.data$Ki67_Manual[which(IM_TMA@Expression$MultiIHC@counts$Ki67 > 15)] <- "1"
table(IM_TMA@meta.data$Ki67_Manual)

Ki67_All_B <- IM_TMA@meta.data$Ki67_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
Ki67_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "Ki67 Expression" = Ki67_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])
unique(Ki67_All_B_data_frame$Sample)

Ki67_percentage_positive <- aggregate(`Ki67.Expression` ~ Sample + `Cell.Type`, 
                                      data = Ki67_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(Ki67_percentage_positive, 
          x = "Cell.Type", 
          y = "Ki67.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for Ki67",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
 theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                       p.adjust.method = "bonferroni", 
                                                                       ref.group = "B cells", 
                                                                       label = "p.adj.signif", tip.length = 0)


# ggbarplot(Ki67_percentage_positive, x = "Cell.Type", y = "Ki67.Expression", fill = "Cell.Type", palette = cell_type_colors, size = 0, add = "jitter",
#           position = position_dodge(0.8), order = cell_type_order,
#           add.params = list(alpha = 0.7), color = "blue", title = "Percentage of B cells positive for Ki67") +
#   ylab("Percentage of Cells Positive for Ki67") +
#   xlab("Cell Type") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#CD44
#%Pos
CD44_All_B <- IM_TMA@meta.data$CD44_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
CD44_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "CD44 Expression" = CD44_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

CD44_percentage_positive <- aggregate(`CD44.Expression` ~ Sample + `Cell.Type`, 
                                      data = CD44_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(CD44_percentage_positive, 
          x = "Cell.Type", 
          y = "CD44.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for CD44",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                    p.adjust.method = "bonferroni",
                                                                    ref.group = "B cells",
                                                                    label = "p.adj.signif", tip.length = 0)


#CD20
#%Pos
CD20_All_B <- IM_TMA@meta.data$CD20_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
CD20_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "CD20 Expression" = CD20_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

CD20_percentage_positive <- aggregate(`CD20.Expression` ~ Sample + `Cell.Type`, 
                                      data = CD20_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(CD20_percentage_positive, 
          x = "Cell.Type", 
          y = "CD20.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for CD20",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                      p.adjust.method = "bonferroni",
                                                                      ref.group = "B cells",
                                                                      label = "p.adj.signif", tip.length = 0)


#CD38
#%Pos
CD38_All_B <- IM_TMA@meta.data$CD38_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
CD38_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "CD38 Expression" = CD38_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

CD38_percentage_positive <- aggregate(`CD38.Expression` ~ Sample + `Cell.Type`, 
                                      data = CD38_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(CD38_percentage_positive, 
          x = "Cell.Type", 
          y = "CD38.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for CD38",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                      p.adjust.method = "bonferroni",
                                                                      ref.group = "B cells",
                                                                      label = "p.adj.signif", tip.length = 0)


#PDL1
#%Pos
PDL1_All_B <- IM_TMA@meta.data$PDL1_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
PDL1_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "PDL1 Expression" = PDL1_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

PDL1_percentage_positive <- aggregate(`PDL1.Expression` ~ Sample + `Cell.Type`, 
                                      data = PDL1_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(PDL1_percentage_positive, 
          x = "Cell.Type", 
          y = "PDL1.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for PDL1",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                      p.adjust.method = "bonferroni",
                                                                      ref.group = "B cells",
                                                                      label = "p.adj.signif", tip.length = 0)

#IDO1
#%Pos
IDO1_All_B <- IM_TMA@meta.data$IDO1_Manual[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
IDO1_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "IDO1 Expression" = IDO1_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

IDO1_percentage_positive <- aggregate(`IDO1.Expression` ~ Sample + `Cell.Type`, 
                                      data = IDO1_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(IDO1_percentage_positive, 
          x = "Cell.Type", 
          y = "IDO1.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for IDO1",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                      p.adjust.method = "bonferroni",
                                                                      ref.group = "B cells",
                                                                      label = "p.adj.signif", tip.length = 0)

#CD30
#%Pos
CD30_All_B <- IM_TMA@meta.data$CD30_Positive[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest]
CD30_All_B_data_frame <- data.frame("Cell Type" = IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest], "CD30 Expression" = CD30_All_B, "Sample" = IM_TMA@meta.data$allRegions[IM_TMA@meta.data$Final_Annotations %in% cell_types_of_interest])

CD30_percentage_positive <- aggregate(`CD30.Expression` ~ Sample + `Cell.Type`, 
                                      data = CD30_All_B_data_frame,
                                      function(x) mean(x == 1) * 100)

# Plotting
ggboxplot(CD30_percentage_positive, 
          x = "Cell.Type", 
          y = "CD30.Expression", 
          fill = "Cell.Type",
          palette = cell_type_colors, order = cell_type_order) +
  labs(title = "Percentage of B cells positive for CD30",
       x = "Cell Type",
       y = "Percentage Positive") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_pwc(method = "wilcox.test", 
                                                                      p.adjust.method = "bonferroni",
                                                                      ref.group = "B cells",
                                                                      label = "p.adj.signif", tip.length = 0)
############################
# IDO1 Threshold = 23.11 or 5

# table(IM_TMA@meta.data[["IDO1_Positive"]])
# table(IM_TMA@meta.data[["IDO1_Manual"]])
# 
# IM_TMA@meta.data$IDO1_Manual[which(IM_TMA@Expression$MultiIHC@counts$IDO1 > 5)] <- "1"
# table(IM_TMA@meta.data$IDO1_Manual)


############################
# Spatial Location of LMP1
############################

spatialPhenotype(MISSILeObject = IM_TMA, region = 1,
                 phenotypes = c("GC B cells", "Epithelium","LMP1+ B cells"), 
                 clustering = "Final_Annotations",
                 colours = c("yellow","black", "red"), pt.size = 0.5)

#################################################
# Cell Type Abundance Per Sample Part of Whole
#################################################

for(i in 1:length(unique(IM_TMA@meta.data$allRegions))){
  
  IM_Cell_Type_Abun <- as.data.frame(cbind(levels(IM_TMA@meta.data$Final_Annotations),cbind(rep(unique(IM_TMA@meta.data$SampleID)[i], times = length(levels(IM_TMA@meta.data$Final_Annotations))),table(IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]]) / length(IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]])))) 
  colnames(IM_Cell_Type_Abun) <- c("Cell.type","Patient","Abundance")
  IM_Cell_Type_Abun$Abundance <- as.numeric(IM_Cell_Type_Abun$Abundance)*100
  
  if(i == 1){
    IM_Cell_Type_Abun_all <- IM_Cell_Type_Abun
  } else{
    IM_Cell_Type_Abun_all <- rbind(IM_Cell_Type_Abun_all,IM_Cell_Type_Abun)
  }
  
}


IM_Cell_Type_Abun <- MISSILe::plotCellAbundanceGrouped(MISSILeObject = IM_TMA, Idents = "Final_Annotations",
                                                       cluster = unique(IM_TMA@meta.data$Final_Annotations), condition = "SampleID")

ggbarplot(IM_Cell_Type_Abun_all_LMP1s, x="Patient", y="Abundance", 
          fill = "Patient", 
          color = "Patient", 
          palette = pallette_2,
          width = 0.85)

# "#FF0000","#00FF00", "#0000FF", "#FFFF00",  "#808000","#FF00FF", "#00FFFF", "#FFA500", "#800080", 
#"#008000", "#C0C0C0", "#008080", "#800000", "#000080", "#808080", "black"

pallette_2 <- c("#92CC71","#B66CCC","#6CB4C7","#C96C69","#E3D0C4","#C89E6D","#706EC4","#BDC7DE","#E7C6DE","#D3E8C7","#B39FD8","#CA80A4","#7ABB9B","#CCC681","#7395C3","#B2E4E2","#73CA6F")

##########################################
# Run neighbourhoods 
##########################################

IM_TMA <- readRDS(".../IDO1_ME.rds")

IM_TMA@current.identity[["ActiveIdents"]]
IM_TMA@current.identity[["ActiveIdents"]] <- "Final_Annotations"

IM_TMA <- MISSILe::cellNeighbourhoods(MISSILeObject = IM_TMA, numOfCells = 10, kMeans = 8, neighbourhoodName = "CellularNeighbourhood8_Select")
MISSILe::neighbourhoodEnrichment(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select") 

spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 1, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 2, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 3, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 4, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 5, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 6, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 7, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 8, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))

IM_NeighFrequ <- neighbourhoodFrequencies(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", comparison = "allRegions")
ggboxplot(IM_NeighFrequ, x = "Region", y = "Frequency", facet.by = "Neighbourhood")

ggbarplot(IM_NeighFrequ, x="Region", y="Frequency", 
          fill = "Neighbourhood", 
          color = "Neighbourhood", 
          palette =  c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"),
          width = 0.85)



cell_type_enrichment <- MISSILe::neighbourhoodCellAbundance(MISSILeObject = IM_TMA, cellType = "Endothelium", neighbourhoodName = "CellularNeighbourhood8_Select")

cell_type_enrichment_long <- gather(cell_type_enrichment, condition, measurement, N1:N8, factor_key=TRUE)

cell_type_enrichment_long <- cell_type_enrichment_long[!(cell_type_enrichment_long$condition %in% c("N6")),]

cell_type_enrichment_long$condition <- as.factor(cell_type_enrichment_long$condition)

cell_type_enrichment_long$condition <- factor(cell_type_enrichment_long$condition, levels = c("N4","N1","N7","N2","N8","N3","N5"))

ggpubr::ggboxplot(cell_type_enrichment_long, x = "condition", y = "measurement") + stat_compare_means(comparisons = list(c("N1","N8"))) + xlab("Neighbourhood")



cell_type_enrichment <- MISSILe::neighbourhoodCellAbundance(MISSILeObject = IM_TMA, cellType = "Smooth muscle", neighbourhoodName = "CellularNeighbourhood8_Select")

cell_type_enrichment_long <- gather(cell_type_enrichment, condition, measurement, N1:N8, factor_key=TRUE)

ggpubr::ggboxplot(cell_type_enrichment_long, x = "condition", y = "measurement")

##########################
# LMP1 and M1 M2 DC correlation 
##########################

IM_Cell_Type_Abun_all_M1s <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "M1 macrophages",]
IM_Cell_Type_Abun_all_M2s <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "M2 macrophages",]
IM_Cell_Type_Abun_all_DCs <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "Dendritic cells",]
IM_Cell_Type_Abun_all_LMP1s <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "LMP1+ B cells",]

M1_Abun <- IM_Cell_Type_Abun_all_M1s$Abundance
M2_Abun <- IM_Cell_Type_Abun_all_M2s$Abundance
DC_Abun <- IM_Cell_Type_Abun_all_DCs$Abundance
LMP1_Abun <- IM_Cell_Type_Abun_all_LMP1s$Abundance

cor.test(M1_Abun, LMP1_Abun, method=c("spearman"))
cor.test(DC_Abun, LMP1_Abun, method=c("spearman"))
cor.test(M2_Abun, LMP1_Abun, method=c("spearman"))

####
# PDL1 IDO1 M1 LMP1 correaltion
###

for(i in 1:length(unique(IM_TMA@meta.data$allRegions))){
  
  IM_Cell_Type_Abun <- as.data.frame(cbind(levels(IM_TMA@meta.data$PDL1_IDO1_M1),cbind(rep(unique(IM_TMA@meta.data$SampleID)[i], times = length(levels(IM_TMA@meta.data$PDL1_IDO1_M1))),table(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]]) / length(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]])))) 
  colnames(IM_Cell_Type_Abun) <- c("Cell.type","Patient","Abundance")
  IM_Cell_Type_Abun$Abundance <- as.numeric(IM_Cell_Type_Abun$Abundance)*100
  
  if(i == 1){
    IM_Cell_Type_Abun_all <- IM_Cell_Type_Abun
  } else{
    IM_Cell_Type_Abun_all <- rbind(IM_Cell_Type_Abun_all,IM_Cell_Type_Abun)
  }
  
}

IM_Cell_Type_Abun_all_LMP1s <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "LMP1+ B cells",]
IM_Cell_Type_Abun_all_PDL1IDO1 <- IM_Cell_Type_Abun_all[IM_Cell_Type_Abun_all$Cell.type == "CD68+ PDL1+ IDO1+ Macrophages",]

PDL1IDO1_Abun <- IM_Cell_Type_Abun_all_PDL1IDO1$Abundance
LMP1_Abun <- IM_Cell_Type_Abun_all_LMP1s$Abundance

cor.test(PDL1IDO1_Abun, LMP1_Abun, method=c("spearman"))

##############################
# CD8:CD4 ratio
##############################

IM_Cell_Type_Abun_all_CD8s <- IM_Cell_Type_Abun[IM_Cell_Type_Abun$Phenotype == "CD8+ T cells",]
IM_Cell_Type_Abun_all_CD4s <- IM_Cell_Type_Abun[IM_Cell_Type_Abun$Phenotype == "CD4+ T cells",]

CD8_CD4_Ratio <- data.frame(Group = rep("CD8/CD4 Ratio", times = 8),
                               Ratio = IM_Cell_Type_Abun_all_CD8s$Frequency/IM_Cell_Type_Abun_all_CD4s$Frequency)

ggdotplot(CD8_CD4_Ratio, x = "Group", y = "Ratio", fill = "black", color = "black", size = 1.05) +
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



#############################
# M1s in LMP1 Neighbourhood
#############################

#LMP1+ cells#
IM_TMA@meta.data$LMP1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$LMP1[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"] <- "LMP1+ B cells"
table(IM_TMA@meta.data$LMP1)

IM_TMA@meta.data$LMP1 <- as.factor(IM_TMA@meta.data$LMP1)
levels(IM_TMA@meta.data$LMP1)
IM_TMA@current.identity[["ActiveIdents"]] <- "LMP1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$LMP1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$LMP1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$LMP1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "LMP1+ B cells",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)





#M1s only#
IM_TMA@meta.data$M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$M1[IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ Macrophages"
table(IM_TMA@meta.data$M1)

IM_TMA@meta.data$M1 <- as.factor(IM_TMA@meta.data$M1)
levels(IM_TMA@meta.data$M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)


#PDL1 only ASSIGN#
IM_TMA@meta.data$PDL1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$PDL1_M1[IM_TMA@meta.data[["PDL1_Positive"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ PDL1+ Macrophages"
table(IM_TMA@meta.data$PDL1_M1)

IM_TMA@meta.data$PDL1_M1 <- as.factor(IM_TMA@meta.data$PDL1_M1)
levels(IM_TMA@meta.data$PDL1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "PDL1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$PDL1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$PDL1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$PDL1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$PDL1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ PDL1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)
#PDL1 only Manual#
table(imTonsil@meta.data$PDL1_Positive)
quantile(IM_TMA@Expression$MultiIHC@counts$PDL1[IM_TMA@meta.data$allRegions == "1"], probs = c(0.1,0.6,0.7,0.75,0.8,0.822,0.9,0.91,0.98,0.984,0.99))
IM_TMA@meta.data$PDL1_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
IM_TMA@meta.data$PDL1_Manual[which(IM_TMA@Expression$MultiIHC@counts$PDL1 > 1.4)] <- "1" #75% because Qupath 25%
table(IM_TMA@meta.data$PDL1_Manual)

IM_TMA@meta.data$PDL1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$PDL1_M1[IM_TMA@meta.data[["PDL1_Manual"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ PDL1+ Macrophages"
table(IM_TMA@meta.data$PDL1_M1)

IM_TMA@meta.data$PDL1_M1 <- as.factor(IM_TMA@meta.data$PDL1_M1)
levels(IM_TMA@meta.data$PDL1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "PDL1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$PDL1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$PDL1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$PDL1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$PDL1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ PDL1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)


#IDO1 only ASSIGN#
IM_TMA@meta.data$IDO1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$IDO1_M1[IM_TMA@meta.data[["IDO1_Positive"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ IDO1+ Macrophages"
table(IM_TMA@meta.data$IDO1_M1)

IM_TMA@meta.data$IDO1_M1 <- as.factor(IM_TMA@meta.data$IDO1_M1)
levels(IM_TMA@meta.data$IDO1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "IDO1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$IDO1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$IDO1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$IDO1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$IDO1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ IDO1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)
#IDO1 only manual#
table(IM_TMA@meta.data$IDO1_Positive)
quantile(IM_TMA@Expression$MultiIHC@counts$IDO1[IM_TMA@meta.data$allRegions == "1"], probs = c(0.1,0.6,0.7,0.75,0.8,0.817,0.85,0.905,0.98,0.984,0.99)) #old 90.05% BECAUSE 9% POS = 35.26
IM_TMA@meta.data$IDO1_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
IM_TMA@meta.data$IDO1_Manual[which(IM_TMA@Expression$MultiIHC@counts$IDO1 > 23.11)] <- "1"
table(IM_TMA@meta.data$IDO1_Manual)

IM_TMA@meta.data$IDO1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$IDO1_M1[IM_TMA@meta.data[["IDO1_Manual"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ IDO1+ Macrophages"
table(IM_TMA@meta.data$IDO1_M1)

IM_TMA@meta.data$IDO1_M1 <- as.factor(IM_TMA@meta.data$IDO1_M1)
levels(IM_TMA@meta.data$IDO1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "IDO1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$IDO1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$IDO1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$IDO1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$IDO1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ IDO1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)


#PDL1 and IDO1 assign#
IM_TMA@meta.data$PDL1_IDO1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@meta.data[["PDL1_Positive"]] == 1 & IM_TMA@meta.data[["IDO1_Positive"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ PDL1+ IDO1+ Macrophages"
table(IM_TMA@meta.data$PDL1_IDO1_M1)

IM_TMA@meta.data$PDL1_IDO1_M1 <- as.factor(IM_TMA@meta.data$PDL1_IDO1_M1)
levels(IM_TMA@meta.data$PDL1_IDO1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "PDL1_IDO1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$PDL1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$PDL1_IDO1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ PDL1+ IDO1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)


#PDL1 and IDO1 manual#
IM_TMA@meta.data$PDL1_IDO1_M1 <- as.character(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@meta.data[["PDL1_Manual"]] == 1 & IM_TMA@meta.data[["IDO1_Manual"]] == 1 & IM_TMA@meta.data$Final_Annotations == "M1 macrophages"] <- "CD68+ PDL1+ IDO1+ Macrophages"
table(IM_TMA@meta.data$PDL1_IDO1_M1)

IM_TMA@meta.data$PDL1_IDO1_M1 <- as.factor(IM_TMA@meta.data$PDL1_IDO1_M1)
levels(IM_TMA@meta.data$PDL1_IDO1_M1)
IM_TMA@current.identity[["ActiveIdents"]] <- "PDL1_IDO1_M1"

neighbourhoodAbundance <- as.data.frame(phonTools::zeros(length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)),length(levels(IM_TMA@meta.data$PDL1_M1))))
colnames(neighbourhoodAbundance) <- levels(IM_TMA@meta.data$PDL1_IDO1_M1)

for(i in c(1:length(unique(IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID)))){
  
  neighbourhoodAbundance[i,] <- (as.numeric(table(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == i])) / length(IM_TMA@meta.data$PDL1_IDO1_M1[IM_TMA@neighbourhoods$generalNeighbourhood$CellularNeighbourhood8_Select@neighbourhoodID == i]))*100
  
}

neighbourhoodAbundance$Neighbourhoods <- paste0("N - ", c(1:8))

#neighbourhoodAbundance <- neighbourhoodAbundance[!(neighbourhoodAbundance$Neighbourhoods == "N - 5"),]



neighbourhoodAbundance_long <- gather(neighbourhoodAbundance, CellTypes, Abundance, `B cells`:`Undefined`, factor_key=TRUE)

neighbourhoodAbundance_long_subset <- neighbourhoodAbundance_long[neighbourhoodAbundance_long$CellTypes == "CD68+ PDL1+ IDO1+ Macrophages",]

ggbarplot(neighbourhoodAbundance_long_subset, x = "Neighbourhoods", y = "Abundance", fill = "Neighbourhoods" ,
          palette = qualpalr::qualpal(length(unique(neighbourhoodAbundance$Neighbourhoods)), "pretty_dark")$hex) + ylab("Cell Abundance [%]") + xlab("") + theme_bw() + theme(panel.grid.major = element_blank(),
                                                                                                                                                                              panel.grid.minor = element_blank(),
                                                                                                                                                                              legend.position = "none") + ggtitle(neighbourhoodAbundance_long_subset$CellTypes)
#####################################
# Microenvironments 
#####################################

IM_TMA@meta.data$Final_Annotations_character <- as.character(IM_TMA@meta.data$Final_Annotations)
cell_types <- levels(IM_TMA@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_TMA@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_TMA@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_TMA <- calculateMicroenvironment2(MISSILeObject = IM_TMA, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_TMA@meta.data$Final_Annotations_character)

###########################
### Downsample for ME ###
###########################

IM_TMA@meta.data$Final_Annotations_ME <- as.character(IM_TMA@meta.data$Final_Annotations)
numBcells <- length(IM_TMA@meta.data$Final_Annotations_ME[IM_TMA@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_TMA@meta.data$Final_Annotations_ME[IM_TMA@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_TMA@meta.data$Final_Annotations_ME[IM_TMA@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_TMA@meta.data$Final_Annotations_ME)

IM_TMA <- calculateMicroenvironment2(MISSILeObject = IM_TMA, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_TMA, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)



plotMicroenvironment(MISSILeObject = IM_TMA, phenotype = "M2 Macr", neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)


######################
# Sample by sample
######################

unique(IM_TMA@meta.data$allRegions)
IM_Core1 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "1")
IM_Core2 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "2")
IM_Core3 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "3")
IM_Core4 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "4")
IM_Core5 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "5")
IM_Core6 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "6")
IM_Core7 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "7")
IM_Core8 <- splitMISSILe(MISSILeObject = IM_TMA, data = "allRegions", criteria = "8")

####IMCORE1

table(IM_Core1@meta.data$Final_Annotations)
IM_Core1@meta.data$Final_Annotations_character <- as.character(IM_Core1@meta.data$Final_Annotations)
cell_types <- levels(IM_Core1@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core1@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core1@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core1 <- calculateMicroenvironment2(MISSILeObject = IM_Core1, phenotypeColumn = "Final_Annotations_character", 
                                     cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core1@meta.data$Final_Annotations_character)

IM_Core1@meta.data$Final_Annotations_ME <- as.character(IM_Core1@meta.data$Final_Annotations)
numBcells <- length(IM_Core1@meta.data$Final_Annotations_ME[IM_Core1@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core1@meta.data$Final_Annotations_ME[IM_Core1@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core1@meta.data$Final_Annotations_ME[IM_Core1@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core1@meta.data$Final_Annotations_ME)

IM_Core1 <- calculateMicroenvironment2(MISSILeObject = IM_Core1, phenotypeColumn = "Final_Annotations_ME", 
                                     cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core1, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)


####IMCORE2

table(IM_Core2@meta.data$Final_Annotations)
IM_Core2@meta.data$Final_Annotations_character <- as.character(IM_Core2@meta.data$Final_Annotations)
cell_types <- levels(IM_Core2@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core2@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core2@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core2 <- calculateMicroenvironment2(MISSILeObject = IM_Core2, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core2@meta.data$Final_Annotations_character)

IM_Core2@meta.data$Final_Annotations_ME <- as.character(IM_Core2@meta.data$Final_Annotations)
numBcells <- length(IM_Core2@meta.data$Final_Annotations_ME[IM_Core2@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core2@meta.data$Final_Annotations_ME[IM_Core2@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core2@meta.data$Final_Annotations_ME[IM_Core2@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core2@meta.data$Final_Annotations_ME)

IM_Core2 <- calculateMicroenvironment2(MISSILeObject = IM_Core2, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core2, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)

#####IMCORE3

table(IM_Core3@meta.data$Final_Annotations)
IM_Core3@meta.data$Final_Annotations_character <- as.character(IM_Core3@meta.data$Final_Annotations)
cell_types <- levels(IM_Core3@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core3@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core3@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core3 <- calculateMicroenvironment2(MISSILeObject = IM_Core3, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core3@meta.data$Final_Annotations_character)

IM_Core3@meta.data$Final_Annotations_ME <- as.character(IM_Core3@meta.data$Final_Annotations)
numBcells <- length(IM_Core3@meta.data$Final_Annotations_ME[IM_Core3@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core3@meta.data$Final_Annotations_ME[IM_Core3@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core3@meta.data$Final_Annotations_ME[IM_Core3@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core3@meta.data$Final_Annotations_ME)

IM_Core3 <- calculateMicroenvironment2(MISSILeObject = IM_Core3, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core3, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)

#### IMCORE4

IM_Core4@meta.data$Final_Annotations_character <- as.character(IM_Core4@meta.data$Final_Annotations)
cell_types <- levels(IM_Core4@meta.data$Final_Annotations)
funcMarkers <- colnames(IM_Core4@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core4@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core4 <- calculateMicroenvironment2(MISSILeObject = IM_Core4, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

IM_Core4@meta.data$Final_Annotations_ME <- as.character(IM_Core4@meta.data$Final_Annotations)
numBcells <- length(IM_Core4@meta.data$Final_Annotations_ME[IM_Core4@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core4@meta.data$Final_Annotations_ME[IM_Core4@meta.data$Final_Annotations_ME == "LMP1+ B cells"])
randomNums <- sample(1:numBcells,numVirusCells)
IM_Core4@meta.data$Final_Annotations_ME[IM_Core4@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core4@meta.data$Final_Annotations_ME)

IM_Core4 <- calculateMicroenvironment2(MISSILeObject = IM_Core4, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core4, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)


#### IMCORE5

table(IM_Core5@meta.data$Final_Annotations)
IM_Core5@meta.data$Final_Annotations_character <- as.character(IM_Core5@meta.data$Final_Annotations)
cell_types <- levels(IM_Core5@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core5@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core5@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core5 <- calculateMicroenvironment2(MISSILeObject = IM_Core5, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core5@meta.data$Final_Annotations_character)

IM_Core5@meta.data$Final_Annotations_ME <- as.character(IM_Core5@meta.data$Final_Annotations)
numBcells <- length(IM_Core5@meta.data$Final_Annotations_ME[IM_Core5@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core5@meta.data$Final_Annotations_ME[IM_Core5@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core5@meta.data$Final_Annotations_ME[IM_Core5@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core5@meta.data$Final_Annotations_ME)

IM_Core5 <- calculateMicroenvironment2(MISSILeObject = IM_Core5, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core5, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)


### IMCORE6

table(IM_Core6@meta.data$Final_Annotations)
IM_Core6@meta.data$Final_Annotations_character <- as.character(IM_Core6@meta.data$Final_Annotations)
cell_types <- levels(IM_Core6@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core6@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core6@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core6 <- calculateMicroenvironment2(MISSILeObject = IM_Core6, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core6@meta.data$Final_Annotations_character)

IM_Core6@meta.data$Final_Annotations_ME <- as.character(IM_Core6@meta.data$Final_Annotations)
numBcells <- length(IM_Core6@meta.data$Final_Annotations_ME[IM_Core6@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core6@meta.data$Final_Annotations_ME[IM_Core6@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core6@meta.data$Final_Annotations_ME[IM_Core6@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core6@meta.data$Final_Annotations_ME)

IM_Core6 <- calculateMicroenvironment2(MISSILeObject = IM_Core6, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core6, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)


### IMSORE7

table(IM_Core7@meta.data$Final_Annotations)
IM_Core7@meta.data$Final_Annotations_character <- as.character(IM_Core7@meta.data$Final_Annotations)
cell_types <- levels(IM_Core7@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core7@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core7@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core7 <- calculateMicroenvironment2(MISSILeObject = IM_Core7, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core7@meta.data$Final_Annotations_character)

IM_Core7@meta.data$Final_Annotations_ME <- as.character(IM_Core7@meta.data$Final_Annotations)
numBcells <- length(IM_Core7@meta.data$Final_Annotations_ME[IM_Core7@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core7@meta.data$Final_Annotations_ME[IM_Core7@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core7@meta.data$Final_Annotations_ME[IM_Core7@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core7@meta.data$Final_Annotations_ME)

IM_Core7 <- calculateMicroenvironment2(MISSILeObject = IM_Core7, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core7, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)

### IMcore8

table(IM_Core8@meta.data$Final_Annotations)
IM_Core8@meta.data$Final_Annotations_character <- as.character(IM_Core8@meta.data$Final_Annotations)
cell_types <- levels(IM_Core8@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core8@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core8@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core8 <- calculateMicroenvironment2(MISSILeObject = IM_Core8, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core8@meta.data$Final_Annotations_character)

IM_Core8@meta.data$Final_Annotations_ME <- as.character(IM_Core8@meta.data$Final_Annotations)
numBcells <- length(IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core8@meta.data$Final_Annotations_ME)

IM_Core8 <- calculateMicroenvironment2(MISSILeObject = IM_Core8, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core8, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)

### IMCore8

table(IM_Core8@meta.data$Final_Annotations)
IM_Core8@meta.data$Final_Annotations_character <- as.character(IM_Core8@meta.data$Final_Annotations)
cell_types <- levels(IM_Core8@meta.data$Final_Annotations)

funcMarkers <- colnames(IM_Core8@Expression$MultiIHC@counts)[c(2,5,7,8,11,15,18,19,20,24,25,27,28,29,30,31,32,33,36,37,40,41,43,45,47,48,49,50,51,52)]
IM_Core8@current.identity$ActiveIdents <- "Final_Annotations_character"

IM_Core8 <- calculateMicroenvironment2(MISSILeObject = IM_Core8, phenotypeColumn = "Final_Annotations_character", 
                                       cellOfInterest = c("LMP1+ B cells"), functionalMarkers = funcMarkers)  

table(IM_Core8@meta.data$Final_Annotations_character)

IM_Core8@meta.data$Final_Annotations_ME <- as.character(IM_Core8@meta.data$Final_Annotations)
numBcells <- length(IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "B cells"])
numVirusCells <- length(IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "LMP1+ B cells"])

randomNums <- sample(1:numBcells,numVirusCells)

IM_Core8@meta.data$Final_Annotations_ME[IM_Core8@meta.data$Final_Annotations_ME == "B cells"][randomNums] <- "Sample B cells"
table(IM_Core8@meta.data$Final_Annotations_ME)

IM_Core8 <- calculateMicroenvironment2(MISSILeObject = IM_Core8, phenotypeColumn = "Final_Annotations_ME", 
                                       cellOfInterest = c("Sample B cells"), functionalMarkers = funcMarkers) 

plotMicroenvironment(MISSILeObject = IM_Core8, phenotype = cell_types, neighbourhoodName = "LMP1+ B cells",
                     cellComparisons = c("LMP1+ B cells","Sample B cells"), colours = c("blue","red"), returnDF = F)



################
# Spatial Expression
###############

spatialExpression(MISSILeObject = IM_TMA, region = 2, marker = "Ki67")







###########
#Functions
###########

downSampleVisualisation <- function(clusters = NULL, downsamplePercent = NULL, countsMatrix = NULL, markers = NULL){
  
  clusterIDs <- levels(clusters)
  
  memIdx <- as.data.frame(as.character(clusters))
  memIdx$index <- rownames(memIdx)
  colnames(memIdx) <- c("membership","index")
  
  if(nrow(memIdx) != nrow(countsMatrix)){
    stop("Stop")
  }
  
  if(length(markers) != length(clusterIDs)){
    stop("Stop 2")
  }
  
  for(i in 1:length(clusterIDs)){
    
    clusterCells <- length(memIdx$membership[memIdx$membership == clusterIDs[i]])
    
    proportion <- round(clusterCells * downsamplePercent)
    
    downDF <- memIdx[memIdx$membership == clusterIDs[i],]
    
    countsMatrixTemp <- countsMatrix[memIdx$membership == clusterIDs[i],]
    
    if(length(markers[[i]]) == 1){
      downDF <- downDF[order(-countsMatrixTemp[,markers[[i]]]),]
    } else if(length(markers[[i]]) == 2){
      downDF <- downDF[order(-countsMatrixTemp[,markers[[i]][1]],-countsMatrixTemp[,markers[[i]][2]]),]
    }
    
    if(i == 1){
      downSampled <- downDF[c(1:proportion),]
    } else {
      tempDownSampled <- downDF[c(1:proportion),]
      downSampled <- rbind(downSampled,tempDownSampled)
      
    }
    
  }
  
  downSampled <- downSampled[ order(as.numeric(downSampled$index)), ]
  
  return(downSampled)
  
}
