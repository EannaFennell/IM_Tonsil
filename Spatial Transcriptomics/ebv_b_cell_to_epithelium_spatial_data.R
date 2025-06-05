library(Seurat)
library(RANN)
library(ggsci)

im_cosmx_ebv_b_cells_2 <- readRDS(".../base_B_cell.rds")

im_cosmx_only_ebv_b <- readRDS(".../only_EBV_B_cell_2.rds")

im_cosmx <- readRDS(".../base_cleaned_IM_clustered.rds")

################################
# Put the EBV clusters back into the original seurat object

im_cosmx@meta.data$clustering.annotations.ebv.b.cells <- as.character(im_cosmx@meta.data$clustering.annotations)

im_cosmx@meta.data$clustering.annotations.ebv.b.cells[im_cosmx@meta.data$cell_id %in% im_cosmx_only_ebv_b@meta.data[["cell_id"]]] <- as.character(im_cosmx_only_ebv_b@meta.data$clustering.base)

phenotypes.of.interest <- c("Epithelial cells","0","1","2","3","4","5","6","7","8")

ebv.metadata <- im_cosmx@meta.data[im_cosmx@meta.data$clustering.annotations.ebv.b.cells %in% phenotypes.of.interest, c("x_slide_mm","y_slide_mm","fov","clustering.annotations.ebv.b.cells")]

ebv.expression <- as.data.frame(t(as.matrix(im_cosmx@assays[["SCT"]]@scale.data[ebv.features, im_cosmx@meta.data$clustering.annotations.ebv.b.cells %in% phenotypes.of.interest])))
ebv.expression <- as.data.frame(t(as.matrix(im_cosmx@assays[["SCT"]]@data[ebv.features, im_cosmx@meta.data$clustering.annotations.ebv.b.cells %in% phenotypes.of.interest])))
ebv.expression <- as.data.frame(t(as.matrix(im_cosmx@assays[["SCT"]]@counts[ebv.features, im_cosmx@meta.data$clustering.annotations.ebv.b.cells %in% phenotypes.of.interest])))

quantile(ebv.expression$EBNA2, probs = c(0, 0.25, 0.5, 0.75, 0.8,0.85,0.9, 0.95, 0.999))

quantile(ebv.expression$BZLF1, probs = c(0, 0.25, 0.5, 0.75, 0.8,0.85,0.9, 0.95, 0.999))

quantile(ebv.expression$LMP1, probs = c(0, 0.25, 0.5, 0.75, 0.8,0.85,0.9, 0.95, 0.999))



colnames(ebv.metadata)[4] <- "Phenotype"


ebv.metadata$easy.phenotype[ebv.expression$LMP1 > 1 & ebv.expression$LMP2A > 1 & ebv.expression$EBNA1 > 1 & ebv.expression$EBNA2 > 1] <- "Lat 3"


ebv.metadata$easy.phenotype <- rep("Other", times = nrow(ebv.metadata))

ebv.metadata$easy.phenotype[ebv.expression$LMP1 > 4] <- "LMP1+"
ebv.metadata$easy.phenotype[ebv.expression$EBNA2 > 3] <- "EBNA2+"
ebv.metadata$easy.phenotype[ebv.expression$EBNA2 > 3 & ebv.expression$BZLF1 > 3] <- "EBNA2+BZLF1+"
ebv.metadata$easy.phenotype[ebv.metadata$Phenotype == "Epithelial cells"] <- "Epithelial cells"

#ebv.metadata <- ebv.metadata[ebv.metadata$easy.phenotype %in% c("Epithelial cells","LMP1+","EBNA2+","EBNA2+BZLF1+"),]
phenotypes.of.interest <- c("Epithelial cells","LMP1+","EBNA2+","EBNA2+BZLF1+","Other")

ebv.metadata$Phenotype <- ebv.metadata$easy.phenotype
ebv.metadata$Phenotype <- as.factor(ebv.metadata$Phenotype)


## Plotting

#colnames(ebv.metadata)[4] <- "Phenotype_Orig"
#colnames(ebv.metadata)[6] <- "Phenotype"

pal_igv("default")(9)


palette <- c("#F0E685FF","#749B58FF","black","#CE3D32FF","#5DB1DDFF") #

c("#BA6338FF","#466983FF","#F0E685FF","#CE3D32FF","#6BD76BFF")

levels(ebv.metadata$Phenotype) <- c("EBV+","BHRF1+","RPMS1+","EBNA-LP+","LMP1+","BARF1+","BNLF2b+","EBNA3A+","BZLF1+","Epithelial")

ggpubr::ggscatter(data = ebv.metadata[ebv.metadata$fov == 17,], x = "x_slide_mm", y = "y_slide_mm", color = "Phenotype",
                  palette = palette) + theme_void() + theme(legend.position = "none") # 5 x 5 / 6 x 6

## Calculations

# loop over regions, get nn2, then filter out into bins, and plot
# Can also do just total distances per EBV phenotype

fovs <- unique(ebv.metadata$fov)

colnames(ebv.metadata)[4] <- "Phenotype"

phenotype.ids <- c(2:5)

a <- 0 

for(i in 1:length(fovs)){
  
  print(i)
  
  if(nrow(ebv.metadata[ebv.metadata$fov == fovs[i] & ebv.metadata$Phenotype == phenotypes.of.interest[1],]) != 0){
  
    a <- a + 1
    
  dists <- RANN::nn2(ebv.metadata[ebv.metadata$fov == fovs[i] & ebv.metadata$Phenotype == phenotypes.of.interest[1],c("x_slide_mm","y_slide_mm")], query = ebv.metadata[ebv.metadata$fov == fovs[i] & ebv.metadata$Phenotype %in% phenotypes.of.interest[phenotype.ids],c("x_slide_mm","y_slide_mm")])
  
  if(a == 1){
    results <- cbind(cbind(ebv.metadata[ebv.metadata$fov == fovs[i] & ebv.metadata$Phenotype %in% phenotypes.of.interest[phenotype.ids],c("Phenotype")],dists$nn.dists[,1]), rep(fovs[i], times = nrow(dists$nn.dists)))
  } else {
    resultsTemp <- cbind(cbind(ebv.metadata[ebv.metadata$fov == fovs[i] & ebv.metadata$Phenotype %in% phenotypes.of.interest[phenotype.ids],c("Phenotype")],dists$nn.dists[,1]), rep(fovs[i], times = nrow(dists$nn.dists)))
    results <- rbind(results,resultsTemp)
  }
  
  }
  
}

results <- as.data.frame(results)
colnames(results) <- c("Phenotype","dist","FOV")
results$dist <- as.numeric(results$dist)

results$Phenotype <- as.factor(results$Phenotype)

ggpubr::ggboxplot(results, x = "Phenotype", y = "dist", facet.by = "FOV")


levels(results$Phenotype) <- c("EBV+","BHRF1+","RPMS1+","EBNA-LP+","LMP1+","BARF1+","BNLF2b+","EBNA3A+","BZLF1+")

results$dist <- results$dist * 1000

ggpubr::ggboxplot(results, x = "Phenotype", y = "dist", 
                  order = c("EBV+","EBNA3A+","BHRF1+","BNLF2b+","BARF1+","LMP1+","EBNA-LP+","RPMS1+","BZLF1+"), 
                  palette = "igv", fill = "Phenotype") + theme(legend.position = "none") + xlab("") + ylab("Distance (um)")



#phenotypes.of.interest <- phenotypes.of.interest[c(2:4)]

#break_amounts <- c(0, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15,0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4) 

break_amounts <- c(0, 0.05, 0.1,0.15,0.2, 0.25, 0.3,0.35, 0.4) 

#break_amounts <- c(0, 0.1, 0.2, 0.3, 0.4) 

#results.binned <- results %>% mutate(new_bin = cut(dist, breaks=break_amounts))
results.binned <- results %>% mutate(new_bin = ntile(dist, n=length(break_amounts)))

for(i in 1:length(break_amounts)){
  if(i == 1){
    bin.results <- data.frame(counts = as.numeric(unname(table(results.binned$Phenotype[results.binned$new_bin == i]))), phenotype = phenotypes.of.interest[phenotype.ids], bin = rep(i, times = 4)) 
  } else{
    bin.results <- rbind(bin.results, data.frame(counts = as.numeric(unname(table(results.binned$Phenotype[results.binned$new_bin == i]))), phenotype = phenotypes.of.interest[phenotype.ids], bin = rep(i, times = 4))) 
  }
}

bin.results$abun <- bin.results$counts


for(i in 2:length(phenotypes.of.interest)){
  bin.results$abun[bin.results$phenotype == phenotypes.of.interest[i]] <- bin.results$counts[bin.results$phenotype == phenotypes.of.interest[i]] / sum(bin.results$counts[bin.results$phenotype == phenotypes.of.interest[i]])
}

bin.results$phenotype <- as.factor(bin.results$phenotype)
levels(bin.results$phenotype) # <- levels(ebv.metadata$Phenotype) 

bin.results$abun <- bin.results$abun * 100

#bin.results <- bin.results[bin.results$phenotype == "EBNA2+",]

palette2 <- c("#F0E685FF","#749B58FF","#CE3D32FF","#5DB1DDFF")

ggplot2::ggplot(data = bin.results, mapping = aes(x = bin, y = abun, group = phenotype, color = phenotype)) + 
  geom_line(size = 1.1) + 
  theme_classic() + ylab("Abundance per phenotype [%]") + theme(text = element_text(size = 16),
                                                                legend.position = "none") + xlab("Distance Bin") + scale_colour_manual(values = palette2)




