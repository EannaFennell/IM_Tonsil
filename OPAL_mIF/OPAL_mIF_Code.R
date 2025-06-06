library(grid)
library(ComplexHeatmap)
library(gridExtra)
library(ggthemes)
library(Seurat)
library(MISSILe)
library(MASS)
library(ggpubr)

Vectra_MISSILe <- readRDS(".../Vectra_MISSILe.rds")

pixelToDiameter <- function(pxArea, pixelSize = 0.5045)
{onePx <- pixelSize * pixelSize
cellArea <- pxArea/onePx
approximateDiameter <- sqrt((4*cellArea)/pi)
return(approximateDiameter)}

##Calculate the cell diameter for each cell
Vectra_MISSILe@meta.data$cellDiameter <- pixelToDiameter(pxArea = Vectra_MISSILe@meta.data$Size, pixelSize = 0.5008)

##Calculate the quantiles of the cell diameters
quantile(Vectra_MISSILe@meta.data$cellDiameter, probs = c(0.001,0.01,0.05,0.1,0.25,0.5,0.75,0.95,0.975,0.995,0.999))

##Exclude big/small cells
diaCut0ffMax <- 80
diaCut0ffMin <- 10
Vectra_MISSILe@Expression$MultiIHC@counts <- Vectra_MISSILe@Expression$MultiIHC@counts[Vectra_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Vectra_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Vectra_MISSILe@Expression$MultiIHC@spatial.data <- Vectra_MISSILe@Expression$MultiIHC@spatial.data[Vectra_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Vectra_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Vectra_MISSILe@meta.data <- Vectra_MISSILe@meta.data[Vectra_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Vectra_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]

##Find % of cells kept 
beforeCutOff <- nrow(cells_before_cutoff@meta.data)
afterCutOff <- sum(Vectra_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Vectra_MISSILe@meta.data$cellDiameter < diaCut0ffMax)
#afterCutOff <- nrow(Vectra_MISSILe@meta.data)
percentageKept <- (afterCutOff/beforeCutOff) * 100
print(percentageKept)

EBV_genes_Vectra <- colnames(Vectra_MISSILe@Expression$MultiIHC@counts)[2:6]

##Set to 0-255
Vectra_MISSILe@Expression$MultiIHC@counts[Vectra_MISSILe@Expression$MultiIHC@counts > 255] <- 255
Vectra_MISSILe@Expression$MultiIHC@counts[Vectra_MISSILe@Expression$MultiIHC@counts < 0] <- 0

##Look at expression
MISSILe::spatialExpression(MISSILeObject = Vectra_MISSILe, region = 1, marker = "EBNA1", pt.size = 0.5, colour = "black")

############################### 
# Make MISSILe from TRUE/FALSE mask ENKTL
############################### 
Mask_InOut <- importSeparateFiles(inputDirectory = "D:/PRINCE/IMTonsil/Vectra/Analysis/Masks/finalFCS/", multipleTMAs = FALSE, filetype = "csv")

for(i in 2:length(Mask_InOut)){colnames(Mask_InOut[[i]]) <- colnames(Mask_InOut[[1]])}

Mask_InOut_MISSILe <- MISSILe::createMISSILe(listFCS = Mask_InOut, markerChannels = c(15:21),
                                             spatialData = c(3,4), extraMetaDataColumns = c(1:2,5:14,22))

colnames(Mask_InOut_MISSILe@Expression$MultiIHC@counts)

length(Vectra_MISSILe@meta.data$Cell_ID) 
length(Mask_InOut_MISSILe@meta.data$Cell_ID)

table(Vectra_MISSILe@meta.data$allRegions)
table(Mask_InOut_MISSILe@meta.data$allRegions)
table(Mask_InOut_MISSILe@meta.data$Exclude)

Mask_InOut_MISSILe@meta.data$cellDiameter <- pixelToDiameter(pxArea = Mask_InOut_MISSILe@meta.data$Size, pixelSize = 0.5008)
quantile(Mask_InOut_MISSILe@meta.data$cellDiameter,na.rm = T,probs = c(0.01,0.25,0.5,0.75,0.95,0.995))
diaCut0ffMax <- 80
diaCut0ffMin <- 10
Mask_InOut_MISSILe@Expression$MultiIHC@counts <- Mask_InOut_MISSILe@Expression$MultiIHC@counts[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Mask_InOut_MISSILe@Expression$MultiIHC@spatial.data <- Mask_InOut_MISSILe@Expression$MultiIHC@spatial.data[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
Mask_InOut_MISSILe@meta.data <- Mask_InOut_MISSILe@meta.data[Mask_InOut_MISSILe@meta.data$cellDiameter > diaCut0ffMin & Mask_InOut_MISSILe@meta.data$cellDiameter < diaCut0ffMax,]
saveRDS(Mask_InOut_MISSILe, ".../Mask_InOut_MISSILe.rds")

length(Vectra_MISSILe@meta.data$Cell_ID) 
length(Mask_InOut_MISSILe@meta.data$Cell_ID)

Vectra_MISSILe@meta.data$Exclude <- Mask_InOut_MISSILe@meta.data$Exclude
table(Vectra_MISSILe@meta.data$Exclude)
saveRDS(Vectra_MISSILe, ".../Vectra_AfterMasks.rds")

Vectra_masked <- splitMISSILe(MISSILeObject = Vectra_MISSILe, data = "Exclude", criteria = "False")
Vectra_exclude <- splitMISSILe(MISSILeObject = Vectra_MISSILe, data = "Exclude", criteria = "True")

##Look at expression
MISSILe::spatialExpression(MISSILeObject = Vectra_exclude, region = 8, marker = "LMP1", pt.size = 0.5, colour = "black")
MISSILe::spatialExpression(MISSILeObject = Vectra_masked, region = 8, marker = "LMP1", pt.size = 0.5, colour = "black")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")
saveRDS(Vectra_exclude, ".../Vectra_Exclude.rds")

Vectra_masked <- readRDS(".../Vectra_Masked.rds")

Vectra_masked@Expression$MultiIHC@counts$LMP1[Vectra_masked@meta.data$allRegions == "6"] <- (Vectra_masked@Expression$MultiIHC@counts$LMP1[Vectra_masked@meta.data$allRegions == "6"]) - 10
Vectra_masked@Expression$MultiIHC@counts$LMP1[Vectra_masked@meta.data$allRegions == "6"] <- (Vectra_masked@Expression$MultiIHC@counts$LMP1[Vectra_masked@meta.data$allRegions == "6"]) - 10
table(Vectra_masked@meta.data$allRegions[Vectra_masked@Expression$MultiIHC@counts$LMP1 < 0])
Vectra_masked@Expression$MultiIHC@counts$LMP1[Vectra_masked@Expression$MultiIHC@counts$LMP1 < 0] <- 0

MISSILe::spatialExpression(MISSILeObject = Vectra_masked, region = 6, marker = "LMP1", pt.size = 0.5, colour = "black")
saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

Vectra_masked <- assignPositivity(MISSILeObject = Vectra_masked, markers = EBV_genes_Vectra, replicateTimes = 5)
table(Vectra_masked@meta.data$LMP1_Positive)
saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

### Manual LMP1
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "LMP1_Positive")
Vectra_masked@meta.data$LMP1_Manual <- Vectra_masked@meta.data$LMP1_Positive
Vectra_masked@meta.data$LMP1_Manual <-"0"
table(Vectra_masked@meta.data$LMP1_Manual)

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "1" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 25] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "2" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 60] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 2, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "3" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 25] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "4" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 30] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "5" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 40] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 5, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "6" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 70] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 6, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "7" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 55] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "1", clustering = "LMP1_Manual")

Vectra_masked@meta.data$LMP1_Manual[Vectra_masked@meta.data$allRegions == "8" & Vectra_masked@Expression$MultiIHC@counts$LMP1 > 15] <- "1"
table(Vectra_masked@meta.data$LMP1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 8, phenotypes = "1", clustering = "LMP1_Manual")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

### Manual LMP2
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "LMP2_Positive")
Vectra_masked@meta.data$LMP2_Manual <- Vectra_masked@meta.data$LMP2_Positive
Vectra_masked@meta.data$LMP2_Manual <-"0"
table(Vectra_masked@meta.data$LMP2_Manual)

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "1" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 40] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "2" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 45] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 2, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "3" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 45] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "4" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 60] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "5" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 60] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 5, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "6" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 40] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 6, phenotypes = "1", clustering = "LMP2_Manual")

#### Mask not perfect? ###
Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "7" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 70] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "1", clustering = "LMP2_Manual")

Vectra_masked@meta.data$LMP2_Manual[Vectra_masked@meta.data$allRegions == "8" & Vectra_masked@Expression$MultiIHC@counts$LMP2 > 70] <- "1"
table(Vectra_masked@meta.data$LMP2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 8, phenotypes = "1", clustering = "LMP2_Manual")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

#### Manual EBNA1
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "EBNA1_Positive")
Vectra_masked@meta.data$EBNA1_Manual <- Vectra_masked@meta.data$EBNA1_Positive
Vectra_masked@meta.data$EBNA1_Manual <-"0"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialExpression(MISSILeObject = Vectra_masked, region = 1, marker = "EBNA1")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "1" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 80] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "2" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 80] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 2, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "3" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 120] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "4" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 140] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "5" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 60] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 5, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "6" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 50] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 6, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "7" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 60] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "1", clustering = "EBNA1_Manual")

Vectra_masked@meta.data$EBNA1_Manual[Vectra_masked@meta.data$allRegions == "8" & Vectra_masked@Expression$MultiIHC@counts$EBNA1 > 80] <- "1"
table(Vectra_masked@meta.data$EBNA1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 8, phenotypes = "1", clustering = "EBNA1_Manual")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

#### Manual EBNA2
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "EBNA2_Positive")
Vectra_masked@meta.data$EBNA2_Manual <- Vectra_masked@meta.data$EBNA2_Positive
Vectra_masked@meta.data$EBNA2_Manual <-"0"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialExpression(MISSILeObject = Vectra_masked, region = 1, marker = "EBNA2")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "1" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 140] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "2" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 110] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 2, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "3" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 120] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "4" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 100] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "5" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 30] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 5, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "6" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 50] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 6, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "7" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 40] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "1", clustering = "EBNA2_Manual")

Vectra_masked@meta.data$EBNA2_Manual[Vectra_masked@meta.data$allRegions == "8" & Vectra_masked@Expression$MultiIHC@counts$EBNA2 > 50] <- "1"
table(Vectra_masked@meta.data$EBNA2_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 8, phenotypes = "1", clustering = "EBNA2_Manual")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

#### Manual BZLF1
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "BZLF1_Positive")
Vectra_masked@meta.data$BZLF1_Manual <- Vectra_masked@meta.data$BZLF1_Positive
Vectra_masked@meta.data$BZLF1_Manual <-"0"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialExpression(MISSILeObject = Vectra_masked, region = 1, marker = "BZLF1")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "1" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 40] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 1, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "2" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 20] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 2, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "3" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 20] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "4" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 40] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "5" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 25] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 5, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "6" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 40] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 6, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "7" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 80] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "1", clustering = "BZLF1_Manual")

Vectra_masked@meta.data$BZLF1_Manual[Vectra_masked@meta.data$allRegions == "8" & Vectra_masked@Expression$MultiIHC@counts$BZLF1 > 60] <- "1"
table(Vectra_masked@meta.data$BZLF1_Manual)
spatialPhenotype(MISSILeObject = Vectra_masked, region = 8, phenotypes = "1", clustering = "BZLF1_Manual")

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")
Vectra_masked <- readRDS(".../Vectra_Masked.rds")
saveRDS(Vectra_masked, ".../Vectra_MISSILe.rds")

####
# Assign Pos Cells
####

LMP1_positive<- Vectra_masked@meta.data[["LMP1_Manual"]] == "1"
EBNA2_positive <- Vectra_masked@meta.data[["EBNA2_Manual"]] == "1"
EBNA1_positive <- Vectra_masked@meta.data[["EBNA1_Manual"]] == "1"
LMP2_positive <- Vectra_masked@meta.data[["LMP2_Manual"]] == "1"
BZLF1_positive <- Vectra_masked@meta.data[["BZLF1_Manual"]] == "1"


#######
# All
#######

All_genes_positive <- EBNA1_positive & EBNA2_positive & LMP1_positive & LMP2_positive & BZLF1_positive
table(All_genes_positive)

counts_5_gene <- c(
  count_cells(All_genes_positive))

######
# FOUR GENES 
######

#Not BZLF1
Not_BZLF1 <- EBNA1_positive &
  EBNA2_positive &
  LMP1_positive &
  LMP2_positive &
  !BZLF1_positive

#Not LMP2
Not_LMP2 <- EBNA1_positive &
  EBNA2_positive &
  LMP1_positive &
  BZLF1_positive &
  !LMP2_positive

#Not EBNA1
Not_EBNA1 <- BZLF1_positive &
  EBNA2_positive &
  LMP1_positive &
  LMP2_positive &
  !EBNA1_positive

#Not EBNA2
Not_EBNA2 <- EBNA1_positive &
  LMP2_positive &
  LMP1_positive &
  BZLF1_positive &
  !EBNA2_positive

#Not LMP1
Not_LMP1 <- EBNA1_positive &
  LMP2_positive &
  EBNA2_positive &
  BZLF1_positive &
  !LMP1_positive

##########
# THREE GENES
#########

Not_BZLF1_Not_LMP2 <- EBNA1_positive &
  EBNA2_positive &
  LMP1_positive &
  !BZLF1_positive &
  !LMP2_positive

Not_BZLF1_Not_EBNA1 <- EBNA2_positive &
  LMP1_positive &
  LMP2_positive &
  !BZLF1_positive &
  !EBNA1_positive

Not_BZLF1_Not_EBNA2 <- EBNA1_positive &
  LMP1_positive &
  LMP2_positive &
  !BZLF1_positive &
  !EBNA2_positive

Not_BZLF1_Not_LMP1 <- EBNA1_positive &
  EBNA2_positive &
  LMP2_positive &
  !BZLF1_positive &
  !LMP1_positive

Not_LMP2_Not_EBNA1 <- EBNA2_positive &
  LMP1_positive &
  BZLF1_positive &
  !LMP2_positive &
  !EBNA1_positive

Not_LMP2_Not_EBNA2 <- EBNA1_positive &
  LMP1_positive &
  BZLF1_positive &
  !LMP2_positive &
  !EBNA2_positive

Not_LMP2_Not_LMP1 <- EBNA1_positive &
  EBNA2_positive &
  BZLF1_positive &
  !LMP2_positive &
  !LMP1_positive

Not_EBNA1_Not_EBNA2 <- LMP1_positive &
  LMP2_positive &
  BZLF1_positive &
  !EBNA1_positive &
  !EBNA2_positive

Not_EBNA1_Not_LMP1 <- EBNA2_positive &
  LMP2_positive &
  BZLF1_positive &
  !EBNA1_positive &
  !LMP1_positive

Not_EBNA2_Not_LMP1 <- EBNA1_positive &
  LMP2_positive &
  BZLF1_positive &
  !EBNA2_positive &
  !LMP1_positive

########
# TWO GENES
########

LMP1_LMP2 <- LMP1_positive &
  LMP2_positive &
  !BZLF1_positive &
  !EBNA1_positive &
  !EBNA2_positive

LMP1_BZLF1 <- LMP1_positive &
  BZLF1_positive &
  !LMP2_positive &
  !EBNA1_positive &
  !EBNA2_positive

LMP1_EBNA1 <- LMP1_positive &
  EBNA1_positive &
  !LMP2_positive &
  !BZLF1_positive &
  !EBNA2_positive

LMP1_EBNA2 <- LMP1_positive &
  EBNA2_positive &
  !LMP2_positive &
  !BZLF1_positive &
  !EBNA1_positive

LMP2_BZLF1 <- LMP2_positive &
  BZLF1_positive &
  !LMP1_positive &
  !EBNA1_positive &
  !EBNA2_positive

LMP2_EBNA1 <- LMP2_positive &
  EBNA1_positive &
  !LMP1_positive &
  !BZLF1_positive &
  !EBNA2_positive

LMP2_EBNA2 <- LMP2_positive &
  EBNA2_positive &
  !LMP1_positive &
  !BZLF1_positive &
  !EBNA1_positive

BZLF1_EBNA1 <- BZLF1_positive &
  EBNA1_positive &
  !LMP1_positive &
  !LMP2_positive &
  !EBNA2_positive

BZLF1_EBNA2 <- BZLF1_positive &
  EBNA2_positive &
  !LMP1_positive &
  !LMP2_positive &
  !EBNA1_positive

EBNA1_EBNA2 <- EBNA1_positive &
  EBNA2_positive &
  !LMP1_positive &
  !LMP2_positive &
  !BZLF1_positive

###
# Pos for one gene
###

LMP1 <- LMP1_positive &
  !LMP2_positive &
  !BZLF1_positive &
  !EBNA1_positive &
  !EBNA2_positive

LMP2 <- LMP2_positive &
  !LMP1_positive &
  !BZLF1_positive &
  !EBNA1_positive &
  !EBNA2_positive

BZLF1 <- BZLF1_positive &
  !LMP1_positive &
  !LMP2_positive &
  !EBNA1_positive &
  !EBNA2_positive

EBNA1 <- EBNA1_positive &
  !LMP1_positive &
  !LMP2_positive &
  !BZLF1_positive &
  !EBNA2_positive

EBNA2 <- EBNA2_positive &
  !LMP1_positive &
  !LMP2_positive &
  !BZLF1_positive &
  !EBNA1_positive

###################
# GRAPH 4 genes
###################

count_cells <- function(positive_combination) {
  sum(positive_combination)
}

counts_4_genes <- c(
  count_cells(Not_BZLF1),
  count_cells(Not_LMP2),
  count_cells(Not_EBNA1),
  count_cells(Not_EBNA2),
  count_cells(Not_LMP1))


barplot(counts_4_genes, 
        names.arg = c("Not_BZLF1", "Not_LMP2", "Not_EBNA1", "Not_EBNA2", "Not_LMP1"), 
        xlab = "Combinations", ylab = "Cell Counts", col = "skyblue",
        ylim = c(0, 120))  # Set the y-axis limits from 0 to 150


###################
# GRAPH 3 genes
###################

counts_3_genes <- c(
  count_cells(Not_BZLF1_Not_LMP2),
  count_cells(Not_BZLF1_Not_EBNA1),
  count_cells(Not_BZLF1_Not_EBNA2),
  count_cells(Not_BZLF1_Not_LMP1),
  count_cells(Not_LMP2_Not_EBNA1),
  count_cells(Not_LMP2_Not_EBNA2),
  count_cells(Not_LMP2_Not_LMP1),
  count_cells(Not_EBNA1_Not_EBNA2),
  count_cells(Not_EBNA1_Not_LMP1),
  count_cells(Not_EBNA2_Not_LMP1)
)

par(mar = c(12, 5, 5, 5)) #(bottom, left, top, right)

barplot(counts_3_genes, 
        names.arg = c("Not_BZLF1_Not_LMP2", "Not_BZLF1_Not_EBNA1", "Not_BZLF1_Not_EBNA2", 
                      "Not_BZLF1_Not_LMP1", "Not_LMP2_Not_EBNA1", "Not_LMP2_Not_EBNA2", 
                      "Not_LMP2_Not_LMP1", "Not_EBNA1_Not_EBNA2", "Not_EBNA1_Not_LMP1", 
                      "Not_EBNA2_Not_LMP1"), ylab = "Cell Counts", col = "skyblue",
        ylim = c(0, 400) , las = 2)  


mtext("Combinations", side = 1, line = 11)

##############
# Graph 2 genes
###############

# Create counts for each combination of two genes
counts_2_genes <- c(
  count_cells(LMP1_LMP2),
  count_cells(LMP1_BZLF1),
  count_cells(LMP1_EBNA1),
  count_cells(LMP1_EBNA2),
  count_cells(LMP2_BZLF1),
  count_cells(LMP2_EBNA1),
  count_cells(LMP2_EBNA2),
  count_cells(BZLF1_EBNA1),
  count_cells(BZLF1_EBNA2),
  count_cells(EBNA1_EBNA2)
)

barplot(counts_2_genes, 
        names.arg = c("LMP1_LMP2", "LMP1_BZLF1", "LMP1_EBNA1", "LMP1_EBNA2", 
                      "LMP2_BZLF1", "LMP2_EBNA1", "LMP2_EBNA2", 
                      "BZLF1_EBNA1", "BZLF1_EBNA2", "EBNA1_EBNA2"), ylab = "Cell Counts", col = "skyblue",
        ylim = c(0, 3700), las = 2)  

mtext("Combinations", side = 1, line = 8)

##############
# GRAPH ONE GENE ONLY
#############

counts_1_gene <- c(
  count_cells(LMP1),
  count_cells(LMP2),
  count_cells(BZLF1),
  count_cells(EBNA1),
  count_cells(EBNA2)
)


barplot(counts_1_gene, 
        names.arg = c("LMP1", "LMP2", "BZLF1", "EBNA1", "EBNA2"), ylab = "Cell Counts", col = "skyblue",
        ylim = c(0, 35000), las = 2)  
mtext("Combinations", side = 1, line = 5)

##################

Vectra_MISSILe@meta.data$phenotype_cells <- rep(0, times = nrow(Vectra_MISSILe@meta.data))

Vectra_masked@meta.data$phenotype_cells[All_genes_positive] <- 1
Vectra_masked@meta.data$phenotype_cells[Not_BZLF1] <- 2
Vectra_masked@meta.data$phenotype_cells[Not_LMP2] <- 3
Vectra_masked@meta.data$phenotype_cells[Not_EBNA1] <- 4
Vectra_masked@meta.data$phenotype_cells[Not_EBNA2] <- 5
Vectra_masked@meta.data$phenotype_cells[Not_LMP1] <- 6
Vectra_masked@meta.data$phenotype_cells[Not_BZLF1_Not_LMP2] <- 7
Vectra_masked@meta.data$phenotype_cells[Not_BZLF1_Not_EBNA1] <- 8
Vectra_masked@meta.data$phenotype_cells[Not_BZLF1_Not_EBNA2] <- 9
Vectra_masked@meta.data$phenotype_cells[Not_BZLF1_Not_LMP1] <- 10
Vectra_masked@meta.data$phenotype_cells[Not_LMP2_Not_EBNA1] <- 11
Vectra_masked@meta.data$phenotype_cells[Not_LMP2_Not_EBNA2] <- 12
Vectra_masked@meta.data$phenotype_cells[Not_LMP2_Not_LMP1] <- 13
Vectra_masked@meta.data$phenotype_cells[Not_EBNA1_Not_EBNA2] <- 14
Vectra_masked@meta.data$phenotype_cells[Not_EBNA1_Not_LMP1] <- 15
Vectra_masked@meta.data$phenotype_cells[Not_EBNA2_Not_LMP1] <- 16
Vectra_masked@meta.data$phenotype_cells[LMP1_LMP2] <- 17
Vectra_masked@meta.data$phenotype_cells[LMP1_BZLF1] <- 18
Vectra_masked@meta.data$phenotype_cells[LMP1_EBNA1] <- 19
Vectra_masked@meta.data$phenotype_cells[LMP1_EBNA2] <- 20
Vectra_masked@meta.data$phenotype_cells[LMP2_BZLF1] <- 21
Vectra_masked@meta.data$phenotype_cells[LMP2_EBNA1] <- 22
Vectra_masked@meta.data$phenotype_cells[LMP2_EBNA2] <- 23
Vectra_masked@meta.data$phenotype_cells[BZLF1_EBNA1] <- 24
Vectra_masked@meta.data$phenotype_cells[BZLF1_EBNA2] <- 25
Vectra_masked@meta.data$phenotype_cells[EBNA1_EBNA2] <- 26
Vectra_masked@meta.data$phenotype_cells[LMP1] <- 27
Vectra_masked@meta.data$phenotype_cells[LMP2] <- 28
Vectra_masked@meta.data$phenotype_cells[BZLF1] <- 29
Vectra_masked@meta.data$phenotype_cells[EBNA1] <- 30
Vectra_masked@meta.data$phenotype_cells[EBNA2] <- 31

table(Vectra_masked@meta.data$phenotype_cells)
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "1"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "2"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "3"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "4"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "5"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "6"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "7"])
table(Vectra_masked@meta.data$phenotype_cells[Vectra_masked@meta.data$allRegions == "8"])

saveRDS(Vectra_masked, ".../Vectra_Masked.rds")

############################

combination_table <- addmargins(table(Vectra_masked@meta.data$phenotype_cells, Vectra_masked@meta.data$allRegions))

unique_phenotypes <- unique(Vectra_masked@meta.data$phenotype_cells)


for (phenotype in unique_phenotypes) {
subset_data <- subset(Vectra_masked@meta.data, phenotype_cells == phenotype)
  
combination_table <- table(subset_data$phenotype_cells, subset_data$allRegions)
  
  
barplot(as.matrix(combination_table), beside = TRUE, col = rainbow(8), 
        xlab = "Region", ylab = "Frequency", 
        main = paste("Phenotype", phenotype),
        legend.text = paste("Region", 1:8),
        args.legend = list(x = "topright", bty = "n", inset = c(0.1, 0.1))
  )
}

table(Vectra_masked@meta.data$allRegions[Vectra_masked@meta.data$phenotype_cells == "3"])
spatialPhenotype(MISSILeObject = Vectra_masked, region = 7, phenotypes = "3", clustering = "phenotype_cells") 

table(Vectra_masked@meta.data$allRegions[Vectra_masked@meta.data$phenotype_cells == "2"])
spatialPhenotype(MISSILeObject = Vectra_masked, region = 4, phenotypes = "2", clustering = "phenotype_cells")

table(Vectra_masked@meta.data$allRegions[Vectra_masked@meta.data$phenotype_cells == "6"])
spatialPhenotype(MISSILeObject = Vectra_masked, region = 3, phenotypes = "6", clustering = "phenotype_cells")

#########################
# NAMED
########################

Vectra_MISSILe@meta.data$phenotype_names <- Vectra_MISSILe@meta.data$phenotype_cells 

Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 1] <- "LMP1+LMP2+EBNA1+EBNA2+BZF1+"

Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 2] <- "LMP1+LMP2+EBNA1+EBNA2+" 
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 3] <- "LMP1+EBNA1+EBNA2+BZLF1+" 
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 4] <- "LMP1+LMP2+EBNA2+BZLF1+"  
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 5] <- "LMP1+LMP2+EBNA1+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 6] <- "LMP2+EBNA1+EBNA2+BZLF1+"

Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 7] <- "LMP1+EBNA1+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 8] <- "LMP1+LMP2+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 9] <- "LMP1+LMP2+EBNA1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 10] <- "LMP2+EBNA1+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 11] <- "LMP1+EBNA2+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 12] <- "LMP1+EBNA1+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 13] <- "EBNA1+EBNA2+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 14] <- "LMP1+LMP2+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 15] <- "LMP2+EBNA2+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 16] <- "LMP2+EBNA1+BZLF1+"

Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 17] <- "LMP1+LMP2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 18] <- "LMP1+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 19] <- "LMP1+EBNA1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 20] <- "LMP1+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 21] <- "LMP2+BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 22] <- "LMP2+EBNA1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 23] <- "LMP2+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 24] <- "BZLF1+EBNA1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 25] <- "BZLF1+EBNA2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 26] <- "EBNA1+EBNA2+"

Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 27] <- "LMP1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 28] <- "LMP2+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 29] <- "BZLF1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 30] <- "EBNA1+"
Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$phenotype_names == 31] <- "EBNA2+"

table(Vectra_MISSILe@meta.data$phenotype_names)
table(Vectra_MISSILe@meta.data$phenotype_cells)

saveRDS(Vectra_MISSILe, ".../Vectra_MISSILe.rds")

##########

cluster_names <- rev(c(
  "LMP1+LMP2+EBNA1+EBNA2+BZLF1+",
  "LMP1+LMP2+EBNA1+EBNA2+",
  "LMP1+EBNA1+EBNA2+BZLF1+",
  "LMP1+LMP2+EBNA2+BZLF1+",
  "LMP1+LMP2+EBNA1+BZLF1+",
  "LMP2+EBNA1+EBNA2+BZLF1+",
  "LMP1+EBNA1+EBNA2+",
  "LMP1+LMP2+EBNA2+",
  "LMP1+LMP2+EBNA1+",
  "LMP2+EBNA1+EBNA2+",
  "LMP1+EBNA2+BZLF1+",
  "LMP1+EBNA1+BZLF1+",
  "EBNA1+EBNA2+BZLF1+",
  "LMP1+LMP2+BZLF1+",
  "LMP2+EBNA2+BZLF1+",
  "LMP2+EBNA1+BZLF1+",
  "LMP1+LMP2+",
  "LMP1+BZLF1+",
  "LMP1+EBNA1+",
  "LMP1+EBNA2+",
  "LMP2+BZLF1+",
  "LMP2+EBNA1+",
  "LMP2+EBNA2+",
  "BZLF1+EBNA1+",
  "BZLF1+EBNA2+",
  "EBNA1+EBNA2+",
  "LMP1+",
  "LMP2+",
  "BZLF1+",
  "EBNA1+",
  "EBNA2+"
))

duplicated_clusters <- cluster_names[duplicated(cluster_names)]
print(duplicated_clusters)

counts_df <- data.frame(
  Combination = c(cluster_names),
  Count = c(counts_5_gene, counts_4_genes, counts_3_genes, counts_2_genes, counts_1_gene)
)

print(counts_df)

counts_df$Log2_Counts <- log2(counts_df$Count)
counts_df$Log2_Counts[counts_df$Log2_Counts == "-Inf"] <- "0"
counts_df$Log2_Counts <- as.numeric(counts_df$Log2_Counts)

ggplot(counts_df, aes(x = reorder(Combination, -Log2_Counts), y = Log2_Counts)) +
geom_bar(stat = "identity", fill = "skyblue", alpha = 0.8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Combination", y = "Log2 Counts", title = "Log2 Counts for Combinations")


#ggbarplot(counts_df, x = "Combination", y = "Log2_Counts", sort.val = "desc")


################
#GC model map 
##############
Vectra_MISSILe <- readRDS(".../Vectra_MISSILe.rds")

Pheno_1 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "1" ]))
Pheno_2 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "2" ]))
Pheno_3 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "3" ]))
Pheno_4 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "4" ]))
Pheno_5 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "5" ]))
Pheno_6 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "6" ]))
Pheno_7 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "7" ]))
Pheno_8 <- names(table(Vectra_MISSILe@meta.data$phenotype_names[Vectra_MISSILe@meta.data$allRegions == "8" ]))

common_values <- Reduce(intersect, list(Pheno_1, Pheno_2, Pheno_3, Pheno_4, Pheno_5, Pheno_6, Pheno_7, Pheno_8))

#plot individual spatial

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Mapping/"
palette <- colorRampPalette(c("red", "white", "blue"))(32)

for (region in 1:8) {
  region_output_dir <- file.path(pdf_output, paste0("region", region))
  
  if (!dir.exists(region_output_dir)) {
    dir.create(region_output_dir, recursive = TRUE)
  }
  
  for (cluster_name in cluster_names) {
    pdf_filename <- file.path(region_output_dir, paste0(gsub("\\+|\\s", "", cluster_name), ".pdf"))
    
    pdf(pdf_filename)
    
    plot_result <- spatialPhenotype(
      MISSILeObject = Vectra_MISSILe,
      region = region,
      clustering = "phenotype_names",
      phenotypes = cluster_name, 
      colours = palette
    )
    
    print(plot_result)
    
    dev.off()
  }
}

spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 6,
  clustering = "phenotype_names",
  phenotypes = order, colours = colours, backgroundCol = "black")

####
# EBV+ cells per patient
####
unique(Vectra_MISSILe@meta.data$phenotype_cells)
table(Vectra_MISSILe@meta.data$phenotype_cells[Vectra_MISSILe@meta.data$allRegions == "1"])

results_df <- data.frame(Patient = numeric(), Sum_NonZero = numeric())

for (region in 1:8) {
  
  patient_cells <- Vectra_MISSILe@meta.data$phenotype_cells[Vectra_MISSILe@meta.data$allRegions == as.character(region)]
  
  sum_nonzero <- (sum(patient_cells != 0) / length(patient_cells)) * 100

  results_df <- rbind(results_df, data.frame(Patient = as.character(region), Sum_NonZero = sum_nonzero))
}

ggbarplot(results_df, x = "Patient", y = "Sum_NonZero", fill = "Patient", 
          palette = "jco", 
          position = position_dodge(0.8), 
          add.params = list(alpha = 0.7), order = c(1,2,4,8,5,6,7,3), title = "% of EBV+ cells per patient") +
  ylab("% of EBV+ Cells") +
  xlab("Patient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####
# Plot "trajectory"
####

library(ggpubr)

Vectra_MISSILe <- readRDS(".../Vectra_MISSILe.rds")

unique(Vectra_MISSILe@meta.data$phenotype_names)
Vectra_MISSILe@meta.data$phenotype_temporal <- as.factor(Vectra_MISSILe@meta.data$phenotype_names)

levels(Vectra_MISSILe@meta.data$phenotype_temporal)

levels(Vectra_MISSILe@meta.data$phenotype_temporal) <- c(0,8,8,1,7,3,0,
                                                         2,6,9,6,0,3,0,
                                                         3,0,5,0,5,3,3,4,
                                                         0,4,3,0,3,0)

Vectra_MISSILe@meta.data$disease <- rep("IM", times = nrow(Vectra_MISSILe@meta.data))

temporal_abundances <- MISSILe::plotCellAbundanceGrouped(MISSILeObject = Vectra_MISSILe, Idents = "phenotype_temporal", 
                                  cluster = levels(Vectra_MISSILe@meta.data$phenotype_temporal), 
                                  condition = "disease")

temporal_abundances <- temporal_abundances[!(temporal_abundances$Phenotype == "0"),]

for (i in 1:length(unique(temporal_abundances$Patient))) {
  temporal_abundances$Frequency[temporal_abundances$Patient == i] <- (temporal_abundances$Counts[temporal_abundances$Patient == i] / sum(temporal_abundances$Counts[temporal_abundances$Patient == i])) * 100
}

temporal_abundances$Phenotype <- as.numeric(temporal_abundances$Phenotype)

ggline(temporal_abundances, "Phenotype", "Frequency", color = "Patient", size = 1.1)

temporal_abundances$Phenotype <- as.factor(temporal_abundances$Phenotype)

ggbarplot(temporal_abundances, "Patient", "Frequency", fill = "Phenotype") + theme(legend.position = "none")

correlation.plot <- as.data.frame(cbind(temporal_abundances$Frequency[temporal_abundances$Phenotype == 2], temporal_abundances$Frequency[temporal_abundances$Phenotype == 6]))
colnames(correlation.plot) <- c("P2","P6")

ggscatter(correlation.plot, x = "P2", y = "P6", title = "Correlation between EBNA2+ and LMP1+ cells states",
          add = "reg.line",                                 
          conf.int = TRUE,                                  
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 3, label.y = 30) + xlab("A") + ylab("B")



####
# Plot "trajectory" of just LMP1 only and EBNA2 only 
####

unique(Vectra_MISSILe@meta.data$phenotype_names)
Vectra_MISSILe@meta.data$phenotype_temporal <- as.factor(Vectra_MISSILe@meta.data$phenotype_names)

levels(Vectra_MISSILe@meta.data$phenotype_temporal)

levels(Vectra_MISSILe@meta.data$phenotype_temporal) <- c(0,8,8,1,7,3,0,
                                                         2,6,9,0,0,3,0,
                                                         3,0,5,0,5,3,3,4,
                                                         0,4,3,0,3,0)

Vectra_MISSILe@meta.data$disease <- rep("IM", times = nrow(Vectra_MISSILe@meta.data))

temporal_abundances <- MISSILe::plotCellAbundanceGrouped(MISSILeObject = Vectra_MISSILe, Idents = "phenotype_temporal", 
                                                         cluster = levels(Vectra_MISSILe@meta.data$phenotype_temporal), 
                                                         condition = "disease")

temporal_abundances <- temporal_abundances[!(temporal_abundances$Phenotype == "0"),]

for (i in 1:length(unique(temporal_abundances$Patient))) {
  temporal_abundances$Frequency[temporal_abundances$Patient == i] <- (temporal_abundances$Counts[temporal_abundances$Patient == i] / sum(temporal_abundances$Counts[temporal_abundances$Patient == i])) * 100
}

temporal_abundances$Phenotype <- as.numeric(temporal_abundances$Phenotype)

ggline(temporal_abundances, "Phenotype", "Frequency", color = "Patient", size = 1.1)

correlation.plot <- as.data.frame(cbind(temporal_abundances$Frequency[temporal_abundances$Phenotype == 2], temporal_abundances$Frequency[temporal_abundances$Phenotype == 6]))
colnames(correlation.plot) <- c("P2","P6")

ggscatter(correlation.plot, x = "P2", y = "P6", title = "Correlation between EBNA2+ and LMP1+ cells",
          add = "reg.line",                                 
          conf.int = TRUE,                                  
          add.params = list(color = "blue",
                            fill = "lightgray")
)+
  stat_cor(method = "spearman", label.x = 3, label.y = 30) + xlab("A") + ylab("B")




order <- c("BZLF1+EBNA2+", 
           "EBNA2+",       
           "LMP1+LMP2+EBNA1+EBNA2+", "LMP1+EBNA1+EBNA2+", "LMP1+LMP2+EBNA2+", "EBNA1+EBNA2+", "LMP1+EBNA2+", "LMP2+EBNA2+", 
           "LMP2+EBNA1+", "LMP2+", 
           "LMP1+LMP2+EBNA1+", "LMP1+LMP2+",
           "LMP1+EBNA1+", 
           "EBNA1+",  
           "BZLF1+", "BZLF1+EBNA1+",
           "LMP1+BZLF1+")

          
colours <- c("green",
            "#FFA500",
            "#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000",
            "#FFA4A4","#FFA4A4",
            "#D5D5FF","#D5D5FF",
            "#4141FF",
            "cyan",
            "purple","purple",
            "#FFFF00")
            
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 6,
  clustering = "phenotype_names",
  phenotypes = order, colours = colours, backgroundCol = "black")

#####
# Test
####

#1#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = "BZLF1+EBNA2+", colours = "green", backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_1.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#2#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = "EBNA2+", colours = "#FFA500", backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_2.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#3#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = c("LMP1+LMP2+EBNA1+EBNA2+", "LMP1+EBNA1+EBNA2+", 
                 "LMP1+LMP2+EBNA2+", "EBNA1+EBNA2+", 
                 "LMP1+EBNA2+", "LMP2+EBNA2+"), colours = c("#FF0000","#FF0000","#FF0000","#FF0000","#FF0000","#FF0000"), backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_3.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#4#  
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = c("LMP2+EBNA1+", "LMP2+"), colours = c("#FFA4A4","#FFA4A4"), backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_4.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#5#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = c("LMP1+LMP2+EBNA1+", "LMP1+LMP2+"), colours = c("#D5D5FF","#D5D5FF"), backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_5.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#6#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = "LMP1+EBNA1+", colours = "#4141FF", backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_6.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#7#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = "EBNA1+", colours = "cyan", backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_7.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#8#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = c("BZLF1+", "BZLF1+EBNA1+"), colours = c("purple","purple"), backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_8.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()

#9#
spatialPhenotype(
  MISSILeObject = Vectra_MISSILe,
  region = 8,
  clustering = "phenotype_names",
  phenotypes = "LMP1+BZLF1+", colours = "#FFFF00", backgroundCol = "black", pt.size = 0.5)

pdf_output <- "D:/PRINCE/IMTonsil/Vectra/Analysis/Co-expression/Co-mapping/"
pdf_filename <- file.path(pdf_output, "plot_9.pdf")
dev.copy2pdf(file = pdf_filename)
dev.off()











# colours = "#FF0000","#FF1010","#FF2020","#FF3131","#FF4141","#FF5252","#FF6262","#FF7373","#FF8383","#FF9494","#FFA4A4","#FFB4B4",
# "#FFC5C5","#FFD5D5","#FFE6E6","#FFF6F6","#F6F6FF","#E6E6FF","#C5C5FF","#B4B4FF","#A4A4FF","#9494FF","#8383FF",
# "#7373FF","#6262FF","#5252FF","#4141FF","#3131FF","#2020FF","#1010FF","#0000FF")
# 
############################

#test_1b <- read.csv("D:/PRINCE/IMTonsil/Vectra/Analysis_1B/1B_test.csv", header = TRUE)
#test_1b <- read.csv("D:/PRINCE/IMTonsil/Vectra/Analysis_1A/CellSeg/fcsFiles/1A.csv", header = TRUE)
test_1a <- read.csv("D:/PRINCE/IMTonsil/Vectra/Analysis_1A/fcsFiles/1A.csv", header = TRUE)

myPalette <- colorRampPalette(c("grey99", "#1BB59A"))(n = 149)
sc <- scale_colour_gradientn(colours = myPalette, limits = c(min(test_1b$EBNA1), 
                                                             max(test_1b$EBNA1)))

ggplot() + geom_point(data = test_1a, aes(x = x, y = y, 
                                             colour = EBNA1), size = 1) + theme_bw() + 
  theme(axis.line = element_blank(), axis.text.x = element_blank(), 
        axis.text.y = element_blank(), axis.ticks = element_blank(), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank()) + ylab("") + xlab("") + 
  sc + coord_fixed()






