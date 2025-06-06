library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(MISSILe)
library(ggpubr)
library(FastPG)
library(grid)
library(ComplexHeatmap)

############################### 
# Import Latest MISSILe Object
############################### 

IM_TMA_All_Cells <- readRDS(".../IM_Missile_Neighbourhoods.rds")
IM_TMA <- readRDS(".../TMA_Select_Final.rds")

############################### 
# Import TRUE/FALSE mask
############################### 

Mask_InOut <- importSeparateFiles(inputDirectory = "D:/PRINCE/IMTonsil/fcs_masked", multipleTMAs = FALSE, filetype = "csv")
colnames(Mask_InOut[[1]])
for(i in 2:length(Mask_InOut)){colnames(Mask_InOut[[i]]) <- colnames(Mask_InOut[[1]])}

############################### 
# Make MISSILe from TRUE/FALSE mask and remove the same cells as before
############################### 

Mask_InOut_MISSILe <- MISSILe::createMISSILe(listFCS = Mask_InOut, markerChannels = c(15:66),
                                             spatialData = c(3,4), extraMetaDataColumns = c(1:2,5:14,67))
colnames(Mask_InOut_MISSILe@Expression$MultiIHC@counts)
pixelToDiameter <- function(pxArea, pixelSize = 0.5045)
{onePx <- pixelSize * pixelSize
cellArea <- pxArea/onePx
approximateDiameter <- sqrt((4*cellArea)/pi)
return(approximateDiameter)}
Mask_InOut_MISSILe@meta.data$cellDiameter <- pixelToDiameter(pxArea = Mask_InOut_MISSILe@meta.data$Size, pixelSize = 1)
quantile(Mask_InOut_MISSILe@meta.data$cellDiameter,na.rm = T,probs = c(0.01,0.25,0.5,0.75,0.95,0.995))
diaCut0ffMax <- 35
diaCut0ffMin <- 7
Mask_InOut_MISSILe@Expression$MultiIHC@counts <- Mask_InOut_MISSILe@Expression$MultiIHC@counts[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Mask_InOut_MISSILe@Expression$MultiIHC@spatial.data <- Mask_InOut_MISSILe@Expression$MultiIHC@spatial.data[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Mask_InOut_MISSILe@meta.data <- Mask_InOut_MISSILe@meta.data[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
saveRDS(Mask_InOut_MISSILe, ".../Mask_InOut_MISSILe.rds")

############################### 
# Copy TRUE/FALSE column to existing object and split MISSILe to keep only true 
############################### 

IM_TMA_All_Cells@meta.data$CellsInclude <- Mask_InOut_MISSILe@meta.data$InOut
IM_TMA <- splitMISSILe(MISSILeObject = IM_TMA_All_Cells, data = "CellsInclude", criteria = "True")

saveRDS(IM_TMA_All_Cells,".../IM_TMA_All_Cells.rds")
saveRDS(IM_TMA,".../TMA_Select_Final.rds")

IM_TMA <- readRDS(".../IM_TMA_Select_Cells.rds")

############################### 
# Import TRUE/FALSE mask v2
############################### 

Mask_InOut <- importSeparateFiles(inputDirectory = "D:/PRINCE/IMTonsil/Analysis/TMA_Analysis/fcs_masked_gc", multipleTMAs = FALSE, filetype = "csv")
Mask_InOut <- Mask_InOut[[1]]
#for(i in 2:length(Mask_InOut)){colnames(Mask_InOut[[i]]) <- colnames(Mask_InOut[[1]])}

############################### 
# Make MISSILe from TRUE/FALSE mask and remove the same cells as before v2
############################### 

nrow(IM_TMA@meta.data[IM_TMA@meta.data$allRegions == 2,])

colnames(Mask_InOut)

IM_TMA@meta.data$GC <- rep("False", times = nrow(IM_TMA@meta.data))

toInsert <- merge(IM_TMA@meta.data[IM_TMA@meta.data$allRegions == 2,],Mask_InOut, by.x = "Cell_ID", by.y = "Cell_ID")

ncol(IM_TMA@meta.data[IM_TMA@meta.data$allRegions == 2,]) == ncol(toInsert)

toInsert <- toInsert[,c(1:76,144)]

toInsert <-toInsert[,c(2,1,3:77)]

colnames(toInsert)[77] <- "GC"

colnames(toInsert) <- colnames(IM_TMA@meta.data)

IM_TMA@meta.data[IM_TMA@meta.data$allRegions == 2,] <- toInsert

table(IM_TMA@meta.data$GC)

saveRDS(IM_TMA,".../IM_TMA_Select_Cells_with_GC_mask.rds")

############################### 
# Clustering
############################### 

IM_TMA <- readRDS(".../TMA_Select_Final_Annotations.rds")
saveRDS(IM_TMA, ".../TMA_Select_Final.rds")
all.genes.TMA <- colnames(IM_TMA@Expression[["MultiIHC"]]@counts)

#################################
# Make GC B cell cluster
#################################

table(IM_TMA@meta.data$CD21_Positive)
quantile(IM_TMA@Expression$MultiIHC@counts$CD21, probs = c(0.94,0.94789,0.95)) #AP threshold ~57.051
IM_TMA@meta.data$CD21_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
table(IM_TMA@meta.data$CD21_Manual)
IM_TMA@meta.data$CD21_Manual[which(IM_TMA@Expression$MultiIHC@counts$CD21 > 65)] <- "1"    
table(IM_TMA@meta.data$CD21_Manual)

IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$PhenographClusteringAnnotated3)
table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$CD21_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "B cells"] <- "GC B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$CD21_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"] <- "GC B cells"
table(IM_TMA@meta.data$Final_Annotations)

spatialPhenotype(MISSILeObject = IM_TMA, region = 8, phenotypes = "GC B cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")

##########################
# B cells in Epithelium
##########################

table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD20_Positive == 1] <- "B cells"
table(IM_TMA@meta.data$Final_Annotations)

saveRDS(IM_TMA, ".../TMA_Select_Final_Annotations.rds")


#################
# LMP1
#################

table(IM_TMA@meta.data$LMP1_Positive)
quantile(IM_TMA@Expression$MultiIHC@counts$LMP1, probs = c(0.93,0.932,0.932758517,0.933)) #AP threshold ~31.49
IM_TMA@meta.data$LMP1_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
table(IM_TMA@meta.data$LMP1_Manual)
IM_TMA@meta.data$LMP1_Manual[which(IM_TMA@Expression$MultiIHC@counts$LMP1 > 65)] <- "1"    
table(IM_TMA@meta.data$LMP1_Manual)

table(IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$allRegions == 8 & IM_TMA@meta.data$LMP1_Manual == 1])

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "B cells"] <- "GC B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"] <- "GC B cells"
table(IM_TMA@meta.data$Final_Annotations)

saveRDS(IM_TMA, ".../TMA_Select_Final_Annotations.rds")

#########################
# Recluster
#########################

Recluster <- splitMISSILe(MISSILeObject = IM_TMA, data = "Final_Annotations", criteria = c("B cells", "LMP1+ B cells", "CD8+ T cells"))
table(Recluster@meta.data$Final_Annotations)

Recluster <- MISSILe::clusterMISSILe(MISSILeObject = Recluster, markers = c("CD20","LMP1","CD21","CD8","CD3e","CD45"), numNeighbours = 30, expSet = "counts", clusteringName = "Reclusterer_30") 
MISSILe::BubblePlot(MISSILeObject = Recluster, markers = lineageMarkers,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Reclusterer_30",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5")) + grids(linetype = "dashed")
table(Recluster@meta.data$Reclusterer_30)

Recluster_annotations <- c("Undefined","LMP1+ B cells","B cells","Undefined","B cells","Undefined","GC B cells","Undefined", #7
                           "GC B cells", "Undefined", "CD4+ T cells", "LMP1+ B cells","B cells","CD8+ T cells","GC B cells", #14
                           "Undefined","Undefined","CD8+ T cells","LMP1+ B cells","CD8+ T cells","B cells","Undefined","Undefined", #22
                           "B cells", "B cells", "Recluster","Undefined","Recluster","CD4+ T cells","GC B cells","Undefined","Undefined", #31
                           "Recluster", "Recluster","Recluster","Recluster","Recluster","Recluster","Recluster", "Recluster")

length(levels(Recluster@meta.data$Reclusterer_30))
length(Recluster_annotations)

Recluster@meta.data$Final_Annotations <- Recluster@meta.data$Reclusterer_30
Recluster@meta.data$Final_Annotations <- as.factor(Recluster@meta.data$Final_Annotations)
levels(Recluster@meta.data$Final_Annotations) <- Recluster_annotations
table(Recluster@meta.data$Final_Annotations)

#Again

Recluster_2 <- splitMISSILe(MISSILeObject = Recluster, data = "Final_Annotations", criteria = c("Recluster"))
table(Recluster_2@meta.data$Final_Annotations)

Recluster_2 <- MISSILe::clusterMISSILe(MISSILeObject = Recluster_2, markers = c("CD20","LMP1","CD21","CD8","CD3e","CD45","CD11c","CD14","CD163","CD31","CD4","CD68"), 
                                       numNeighbours = 30, expSet = "counts", clusteringName = "Reclusterer_30_2") 
MISSILe::BubblePlot(MISSILeObject = Recluster_2, markers = lineageMarkers,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Reclusterer_30_2",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5")) + grids(linetype = "dashed")
table(Recluster_2@meta.data$Reclusterer_30_2)

Recluster_2_annotations <- c("CD8+ T cells","LMP1+ B cells","GC B cells","B cells","CD4+ T cells","B cells","Undefined","Undefined", #7
                             "CD8+ T cells", "Undefined", "CD8+ T cells", "CD8+ T cells","M1 macrophages","Dendritic cells","CD8+ T cells", #14
                             "Myeloid","Undefined","GC B cells","Dendritic cells","B cells","M2 macrophages","CD8+ T cells","Endothelium", #22
                             "M2 macrophages", "M2 macrophages", "CD4+ T cells")

length(levels(Recluster_2@meta.data$Reclusterer_30_2))
length(Recluster_2_annotations)

Recluster_2@meta.data$Final_Annotations <- Recluster_2@meta.data$Reclusterer_30_2
Recluster_2@meta.data$Final_Annotations <- as.factor(Recluster_2@meta.data$Final_Annotations)
levels(Recluster_2@meta.data$Final_Annotations) <- Recluster_2_annotations
table(Recluster_2@meta.data$Final_Annotations)

saveRDS(Recluster_2, ".../Recluster_2_23.10.rds")

COI <- c("B cells", "LMP1+ B cells", "CD8+ T cells")
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% COI] <- as.character(Recluster@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations %in% COI] <- as.character(Recluster_2@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Recluster"] <- as.character(Recluster_2@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
levels(IM_TMA@meta.data$Final_Annotations)

MISSILe::BubblePlot(MISSILeObject = IM_TMA, markers = lineageMarkers,
                    threshold.percent = 0.2, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Final_Annotations",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5")) + grids(linetype = "dashed")

spatialPhenotype(MISSILeObject = IM_TMA, region = 4, phenotypes = "LMP1+ B cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")

##############################################
# Final adjustments

##############################################

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$CD4_Positive == 1 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$CD20_Positive == 0] <- "CD4+ T cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$CD8_Positive == 1 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$CD20_Positive == 0] <- "CD8+ T cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$CD20_Positive == 1 &IM_TMA@meta.data$LMP1_Manual == 0] <- "B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$LMP1_Manual == 1] <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "GC B cells" & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "GC B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "B cells"

spatialPhenotype(MISSILeObject = IM_TMA, region = 8, phenotypes = c("LMP1+ B cells","GC B cells","B cells"), clustering = "Final_Annotations", backgroundCol = "black", colours = c("red", "yellow","blue"))
spatialPhenotype(MISSILeObject = IM_TMA, region = 1, phenotypes = "LMP1+ B cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")
spatialPhenotype(MISSILeObject = IM_TMA, region = 1, phenotypes = "CD8+ T cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")
spatialPhenotype(MISSILeObject = IM_TMA, region = 5, phenotypes = "Endothelium", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")

plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "B cells", fill = rainbow(54), yscale = c(0,100))

table(IM_TMA@meta.data$Final_Annotations)

############################
############################
# FINAL
############################
############################ 
# table(IM_TMA@meta.data$LMP1_Positive)
# quantile(IM_TMA@Expression$MultiIHC@counts$LMP1, probs = c(0.93,0.932,0.932758517,0.933)) #AP threshold ~31.49
# IM_TMA@meta.data$LMP1_Manual <- rep("0", times = nrow(IM_TMA@meta.data))
# table(IM_TMA@meta.data$LMP1_Manual)
# IM_TMA@meta.data$LMP1_Manual[which(IM_TMA@Expression$MultiIHC@counts$LMP1 > 65)] <- "1"    
# table(IM_TMA@meta.data$LMP1_Manual)
# 
# table(IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$allRegions == 8 & IM_TMA@meta.data$LMP1_Manual == 1])
# 
# # IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "B cells"] <- "GC B cells"
# # IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"] <- "GC B cells"
# # table(IM_TMA@meta.data$Final_Annotations)
# 
# IM_TMA@Expression$MultiIHC@counts$CD31 <- IM_TMA@Expression$MultiIHC@counts$CD31 -50 
# IM_TMA@Expression$MultiIHC@counts$CD31[IM_TMA@Expression$MultiIHC@counts$CD31 < 0] <- 0
# IM_TMA <- assignPositivity(MISSILeObject = IM_TMA, markers = "CD31", replicateTimes = 5)

saveRDS(IM_TMA, ".../TMA_Select_Final_Annotations.rds")


IM_TMA <- readRDS(".../TMA_Select_Final_Annotations.rds")

#######################
# B cells
#######################
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "B cells", fill = rainbow(54), yscale = c(0,100))
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & 
                                     IM_TMA@meta.data$CD21_Manual == 1] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" &                                  
                                     IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD4_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1]  <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & 
                                     IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & 
                                     IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & 
                                     IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1] <- "Myeloid"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" &                                  
                                     IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD4_Positive == 1 & IM_TMA@meta.data$CD3e_Positive == 1]  <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$CD20_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 &
                                     (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1)]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & 
                                     IM_TMA@meta.data$LMP1_Manual == 1 & 
                                     IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "B cells", fill = rainbow(54), yscale = c(0,100))

#######################
# CD4 cells
#######################

IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "CD4+ T cells", fill = rainbow(54), yscale = c(0,100))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD21_Manual == 1] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD4_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD4_Positive == 0 & (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1) ]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells" 

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD4+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Manual == 1 ] <- "B cells" 

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "CD4+ T cells", fill = rainbow(54), yscale = c(0,100))

#######################
# CD8 cells
#######################
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "CD8+ T cells", fill = rainbow(54), yscale = c(0,100))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$CD21_Manual == 1 ] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD4_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0 ] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1 ] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD4_Positive == 1 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 ] <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells" 

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "CD8+ T cells" & IM_TMA@meta.data$CD3e_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 ] <- "B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "CD8+ T cells", fill = rainbow(54), yscale = c(0,100))

#######################
# DC 
#######################
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Dendritic cells", fill = rainbow(54), yscale = c(0,150))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD20_Positive == 1 & IM_TMA@meta.data$CD21_Manual == 1 ] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0 ] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1 ] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$CD4_Positive == 1] <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1) ]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells" 

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Dendritic cells" & IM_TMA@meta.data$CD11c_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 ] <- "B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Dendritic cells", fill = rainbow(54), yscale = c(0,200))

#######################
# Endothelium 
#######################
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Endothelium", fill = rainbow(54), yscale = c(0,100))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$aSMA_Positive == 1 ] <- "Smooth muscle"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD20_Positive == 1 & IM_TMA@meta.data$CD21_Manual == 1 ] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0 ] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1 ] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$CD4_Positive == 1] <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1) ]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells" 

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Endothelium" & IM_TMA@meta.data$CD31_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 ] <- "B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Endothelium", fill = rainbow(54), yscale = c(0,100))

#######################
# Epithelium
#######################
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Epithelium", fill = rainbow(54), yscale = c(0,100))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD20_Positive == 1 & IM_TMA@meta.data$CD21_Manual == 1 ] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0 ] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1 ] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$CD4_Positive == 1] <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1) ]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "LMP1+ B cells" 

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Epithelium" & IM_TMA@meta.data$PanCytokeratin_Positive == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 ] <- "B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "Epithelium", fill = rainbow(54), yscale = c(0,100))

#######################
# LMP1
#######################
IM_TMA@meta.data$Final_Annotations <- as.character(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "LMP1+ B cells", fill = rainbow(54), yscale = c(0,100))

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD8_Positive == 0 & IM_TMA@meta.data$GrB_Positive == 0 & IM_TMA@meta.data$FoxP3_Positive == 0 & IM_TMA@meta.data$CD3e_Positive == 1 ] <- "Recluster"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 0 ] <- "M1 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD68_Positive == 1 & IM_TMA@meta.data$CD163_Positive == 1 ] <- "M2 macrophages"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD68_Positive == 0 & IM_TMA@meta.data$CD163_Positive == 0 & IM_TMA@meta.data$CD14_Positive == 1 ] <- "Myeloid"                                                                   

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & IM_TMA@meta.data$CD4_Positive == 1] <- "CD4+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD3e_Positive == 1 & (IM_TMA@meta.data$CD8_Positive == 1 | IM_TMA@meta.data$GrB_Positive == 1) ]  <- "CD8+ T cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 & IM_TMA@meta.data$CD21_Manual == 1 ] <- "GC B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells" & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$LMP1_Manual == 0 & IM_TMA@meta.data$CD20_Positive == 1 & IM_TMA@meta.data$CD21_Manual == 0 ] <- "B cells"

IM_TMA@meta.data$Final_Annotations <- as.factor(IM_TMA@meta.data$Final_Annotations)
table(IM_TMA@meta.data$Final_Annotations)
plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "LMP1+ B cells", fill = rainbow(54), yscale = c(0,100))

MISSILe::BubblePlot(MISSILeObject = IM_TMA, markers = lineageMarkers,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Final_Annotations",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5")) + grids(linetype = "dashed")


###########################
# Put all LMP1+ B cells into B cells then take out the ones of certain thresholds
###########################
IM_TMA@current.identity[["ActiveIdents"]] 
IM_TMA@current.identity[["ActiveIdents"]] <- "Final_Annotations"
table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "LMP1+ B cells"] <- "B cells"                                  
table(IM_TMA@meta.data$Final_Annotations)

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "1" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 79]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "2" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 90 & IM_TMA@meta.data$GC == "False"]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "3" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 90 & IM_TMA@meta.data$CD21_Manual == "0"]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "4" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 90]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "5" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 20]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "6" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 20]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "7" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 20]  <- "LMP1+ B cells"
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "8" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 20]  <- "LMP1+ B cells"

IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "2" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 10 & IM_TMA@meta.data$CD21_Manual == "0"]  <- "LMP1+ B cells"
# IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "B cells" & IM_TMA@meta.data$allRegions == "2" & IM_TMA@Expression$MultiIHC@counts$LMP1 > 100 & IM_TMA@meta.data$CD21_Manual == "0" & (IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "1" | IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "2" |IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "3" |IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "4" | IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "5" | IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8_Select"]]@neighbourhoodID == "6")]  <- "LMP1+ B cells"

IM_TMA@meta.data[["GC"]]

spatialPhenotype(MISSILeObject = IM_TMA, region = 1, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 2, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 3, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 4, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 5, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 6, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 7, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)
spatialPhenotype(MISSILeObject = IM_TMA, region = 8, 
                 phenotypes = c("LMP1+ B cells","B cells"), clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green","blue"), pt.size = 0.25)

IM_TMA <- MISSILe::cellNeighbourhoods(MISSILeObject = IM_TMA, numOfCells = 10, kMeans = 8, neighbourhoodName = "CellularNeighbourhood8_Select")
MISSILe::neighbourhoodEnrichment(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select") 

spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 1, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 2, pt.size = 0.05, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"), ignoreNeighbourhoods = "7")
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 3, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 4, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 5, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 6, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 7, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))
spatialNeighbourhood(MISSILeObject = IM_TMA, neighbourhoodName = "CellularNeighbourhood8_Select", region = 8, pt.size = 0.25, colours = c("#DFD561", "#CE7666", "#CA6C70", "#DDA4B3", "#4AA4DE", "#23304E", "#2A7C39", "#295A43"))


plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", clusterID = "B cells", fill = rainbow(60), yscale = c(0,130))


################
# ReImage clustering
################
table(IM_TMA@meta.data$Final_Annotations)
IM_TMA@meta.data$Final_Annotations[IM_TMA@meta.data$Final_Annotations == "Recluster" & IM_TMA@meta.data$CD4_Positive == 1 & IM_TMA@meta.data$CD3e_Positive == 1] <- "CD4+ T cells"
table(IM_TMA@meta.data$Final_Annotations)

MISSILe::BubblePlot(MISSILeObject = IM_TMA, markers = all.genes.TMA,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Final_Annotations",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5")) + grids(linetype = "dashed")

MISSILe::plotClusterIntensities(MISSILeObject = IM_TMA,markers = all.genes.TMA,clusteringName = "Final_Annotations",expSet = "counts")


#####################################
# CHECK
#####################################

spatialPhenotype(MISSILeObject = IM_TMA, region = 4, phenotypes = c("LMP1+ B cells","GC B cells","B cells"), clustering = "Final_Annotations", backgroundCol = "black", colours = c("red", "yellow","blue"))
spatialPhenotype(MISSILeObject = IM_TMA, region = 1, phenotypes = "LMP1+ B cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")
spatialPhenotype(MISSILeObject = IM_TMA, region = 1, phenotypes = "CD8+ T cells", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")
spatialPhenotype(MISSILeObject = IM_TMA, region = 5, phenotypes = "Endothelium", clustering = "Final_Annotations", backgroundCol = "black", colours = "red")



lineageMarkers <- colnames(IM_TMA@Expression$MultiIHC@counts)[c(3,4,6,9,10,12,13,14,16,17,21,22,23,26,34,35,39,42,44)]
MISSILe::BubblePlot(MISSILeObject = IM_TMA, markers = lineageMarkers,
                    threshold.percent = 0, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "Final_Annotations",
                    transposePlot = TRUE, colours = c("#132B43","#2871b5"))

MISSILe::plotClusterIntensities(MISSILeObject = IM_TMA, markers = lineageMarkers,
                                clusteringName = "Final_Annotations")

############################### 
# tSNE, Downsample function from tonsilAnalysis.R
############################### 

colnames(IM_TMA@Expression$MultiIHC@counts)
levels(IM_TMA@meta.data$Final_Annotations)
orderMarkers_IM_TMA <- list(c(39),                       #Bcells
                            c(4,34,13),                  #CD4+Tcells
                            c(17,34),                    #CD8+Tcells
                            c(10,16,13,44),              #DCs
                            c(3),                        #Endothelium
                            c(12),                       #Epithelium
                            c(39,23),                    #GC B
                            c(39,42, 13),                #LMP1+ B cells
                            c(10,26,44,13),              #M1s
                            c(10,6,26),                  #M2s
                            c(26),                       #Myeloid
                            c(44),                       #Neutrophils
                            c(14),                       #Plasma
                            c(9),                        #Smooth muscle
                            c(4,34,21),                  #T regs
                            c(1))                        #Undefined

a <- downSampleVisualisation(clusters = IM_TMA@meta.data$Final_Annotations,
                             downsamplePercent = 0.2,
                             countsMatrix = IM_TMA@Expression$MultiIHC@counts,
                             markers = orderMarkers_IM_TMA)

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

qualpalette <- qualpalr::qualpal(16, "pretty")   # pretty, pretty_dark, rainbow, pastels

ggplot() + geom_point(data = rtsne_plot, aes(x=x, y=y, color=Clusters), alpha = 0.7) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_bw() + labs(y= "tSNE 2", x = "tSNE 1") +
  theme_classic() + scale_color_manual(values = qualpalette$hex) #+
geom_text(data = as.data.frame(centroids.data), mapping = aes(x = centroids_x,
                                                              y = centroids_y,
                                                              label = 1:length(unique(colorIn))),
          color = "black", size = 6, fontface = 2)

############################### 
# UMAP not working
############################### 

UMAPplot(MISSILeObject = IM_TMA, labels = T, 
         orderMarkers = orderMarkers_IM_TMA, 
         downsample = T, 
         downsamplePercent = 0.1, 
         clusteringName = "PhenographClusteringAnnotated3",
         colPal = qualpalette)


############################### 
# Expression Correlation 
############################### 

all.genes.TMA <- colnames(IM_TMA@Expression[["MultiIHC"]]@counts)[c(2:52)]
markerCorrelation(MISSILeObject = IM_TMA, expSet = "counts", markers = all.genes, corrMethod = "pearson")

IM_TMA <- readRDS(".../IM_TMA_Select_Cells.rds")

####################
# Assign Positivity 
####################

all.genes.TMA <- colnames(IM_TMA@Expression[["MultiIHC"]]@counts)[c(2:52)]
IM_TMA <- assignPositivity(MISSILeObject = IM_TMA, markers = all.genes.TMA, replicateTimes = 5)

####################
# Spatial Phenotype 
####################

pallette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                       "#8c564b", "yellow", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                       "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5") 
                       
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 1, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 2, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 3, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 4, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 5, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 6, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 7, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 8, clustering = "Final_Annotations", pt.size = 0.25, colours = pallette)

MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 1, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 2, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 3, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 4, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 5, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 6, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 7, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")
MISSILe::spatialPhenotype(MISSILeObject = IM_TMA, region = 8, clustering = "Final_Annotations", pt.size = 0.25, phenotypes = "Undefined")

########################################
# Cell Type Abundance Per Sample Part of Whole
########################################
for(i in 1:length(unique(IM_TMA@meta.data$allRegions))){
  
  IM_Cell_Type_Abun <- as.data.frame(cbind(levels(IM_TMA@meta.data$PhenographClusteringAnnotated3),cbind(rep(unique(IM_TMA@meta.data$SampleID)[i], times = length(levels(IM_TMA@meta.data$PhenographClusteringAnnotated3))),table(IM_TMA@meta.data$PhenographClusteringAnnotated3[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]]) / length(IM_TMA@meta.data$PhenographClusteringAnnotated3[IM_TMA@meta.data$SampleID == unique(IM_TMA@meta.data$SampleID)[i]])))) 
  colnames(IM_Cell_Type_Abun) <- c("Cell.type","Patient","Abundance")
  IM_Cell_Type_Abun$Abundance <- as.numeric(IM_Cell_Type_Abun$Abundance)*100
  
  if(i == 1){
    IM_Cell_Type_Abun_all <- IM_Cell_Type_Abun
  } else{
    IM_Cell_Type_Abun_all <- rbind(IM_Cell_Type_Abun_all,IM_Cell_Type_Abun)
  }
  
}

ggbarplot(IM_Cell_Type_Abun_all, x="Patient", y="Abundance", 
          fill = "Cell.type", 
          color = "Cell.type", 
          palette = c("#FF0000","#00FF00", "#0000FF", "#FFFF00",  "#808000","#FF00FF", "#00FFFF", "#FFA500", "#800080", 
                      "#008000", "#C0C0C0", "#008080", "#800000", "#000080", "#808080"),
                      width = 0.85)
                       
# palette = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', 
#                      '#f032e6', '#bcf60c', '#fabebe','#800000', '#000075', '#808080','#ffffff', '#000000')

##########################################
# Stash Neighbourhood before masking
##########################################

neighbourhoodEnrichment(MISSILeObject = IM_TMA, neighbourhoodName = CellularNeighbourhood8)
IM_TMA@meta.data$GeneralNeighbourhood1 <- IM_TMA@neighbourhoods[["generalNeighbourhood"]][["CellularNeighbourhood8"]]@neighbourhoodID[rownames(IM_TMA@meta.data)]

##########################################
# Re-run neighbourhoods on select cells
##########################################

IM_TMA@current.identity[["ActiveIdents"]]
IM_TMA@current.identity[["ActiveIdents"]] <- "PhenographClusteringAnnotated3"

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


############################
# Expression on Clusters
############################
B_Genes <- c("CD20", "CD44", "CD21", "CD138", "CD38", "CD40", 
             "HLA.DR", "TIGIT", "CD45RO", "CD30", "LMP1", "IDO1","TIM3", 
             "PD1","TCF7","LAG3", "BCL2", "CD39", "CD69", "CD57","ICOS")

plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", 
                                    clusterID = "B cells", channels = B_Genes, fill = rainbow(52), yscale = c(0,100))

plotSingleClusterExpressionForCiara(MISSILeObject = IM_TMA, clusteringName = "Final_Annotations", 
                                    clusterID = "LMP1+ B cells", channels = B_Genes, fill = rainbow(52), yscale = c(0,120))
                                    
spatialPhenotype(MISSILeObject = IM_TMA, region = )

##########################
# B cells
##########################
mean(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "B cells"])
quantile(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "B cells"], probs = c(0.01,0.25,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
median(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "B cells"])

mean(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "LMP1+ B cells"])
quantile(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "LMP1+ B cells"], probs = c(0.01,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
median(imTonsil_TMA@Expression$MultiIHC@counts$LMP1[imTonsil_TMA@meta.data$PhenographClusteringAnnotated3 == "LMP1+ B cells"])

quantile(imTonsil_TMA@Expression$MultiIHC@counts$LMP1, probs = c(0.75,0.95))

B_Recluster <- splitMISSILe(MISSILeObject = IM_TMA, data = "PhenographClusteringAnnotated3", criteria = c("LMP1+ B cells", "B cells"))

B_cluster <- c("CD20","CD44","CD21","CD38","CD40","CD45RO","CD30","LMP1","BCL2","Ki67")  

B_Recluster <- MISSILe::clusterMISSILe(MISSILeObject = B_Recluster,
                                   markers = B_cluster,
                                   numNeighbours = 60,
                                   expSet = "counts",
                                   clusteringName = "PhenographClustering5")

MISSILe::BubblePlot(MISSILeObject = B_Recluster, markers = lineageMarkers,
                    threshold.percent = 0.05, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "PhenographClustering4",
                    transposePlot = TRUE) + grids(linetype = "dashed")

MISSILe::BubblePlot(MISSILeObject = B_Recluster, markers = all.genes.TMA,
                    threshold.percent = 0.05, ignoreClusters = NULL,
                    colour.scale = c(0,1.5), identities = "PhenographClustering4",
                    transposePlot = TRUE) + grids(linetype = "dashed")

spatialPhenotype(MISSILeObject = B_Recluster, region = 3, 
                 phenotypes = c("7"), clustering = "PhenographClustering4",
                 backgroundCol = "black", colours = c("green","red","blue"), pt.size = 0.25)


B_cell_annotations <- c("Epithelium","Undefined","Undefined","B cells","B cells","Epithelium", #5
                        "GC B cells","Undefined","B cells","Epithelium","M2 macrophages","Undefined", #11
                        "Undefined","B cells","Undefined","GC B cells","GC B cells","CD38+","GC B cells",
                        "GC B cells") #17 


##################################
#LMP1+ Thresholds
##################################
spatialPhenotype(MISSILeObject = IM_TMA, region = 1, 
                 phenotypes = "LMP1+ B cells", clustering = "Final_Annotations",
                 backgroundCol = "black", colours = c("green"), pt.size = 0.25)





































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