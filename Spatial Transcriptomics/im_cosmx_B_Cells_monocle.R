# Monocle

library(monocle3)
library(Seurat)
library(ggplot2)
source(".../asCellDataSet.R")

im_cosmx_ebv_b_cells_2 <- readRDS(".../base_B_cell.rds")

im_cosmx_ebv_b_cells_2_cds <- as.cell_data_set(im_cosmx_ebv_b_cells_2)

im_cosmx_ebv_b_cells_2_cds <- estimate_size_factors(im_cosmx_ebv_b_cells_2_cds) # necessary for running plot_gene_in_pseudotime

fData(im_cosmx_ebv_b_cells_2_cds)$gene_short_name <- rownames(fData(im_cosmx_ebv_b_cells_2_cds))
head(fData(im_cosmx_ebv_b_cells_2_cds))

head(counts(im_cosmx_ebv_b_cells_2_cds))

# Retrieve clustering information from Seurat object

recreate.partitions <- c(rep(1, length(im_cosmx_ebv_b_cells_2_cds@colData@rownames)))
names(recreate.partitions) <- im_cosmx_ebv_b_cells_2_cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions

im_cosmx_ebv_b_cells_2_cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- im_cosmx_ebv_b_cells_2@active.ident
im_cosmx_ebv_b_cells_2_cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

im_cosmx_ebv_b_cells_2_cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- im_cosmx_ebv_b_cells_2@reductions$umap@cell.embeddings

# Plot

cluster.before.traj <- plot_cells(im_cosmx_ebv_b_cells_2_cds, color_cells_by = "cluster", label_groups_by_cluster = F, 
                                  group_label_size = 5, cell_size = 2) + theme(legend.position = "right")

cluster.before.traj

# Learn trajectory - adjust partitions

im_cosmx_ebv_b_cells_2_cds@clusters@listData$UMAP$partitions <- im_cosmx_ebv_b_cells_2_cds@clusters@listData$UMAP$clusters
levels(im_cosmx_ebv_b_cells_2_cds@clusters@listData$UMAP$partitions) <- c("1","2","3","4","5","6","7","8","9") # Needs to be 1,2,3,etc. for whatever reason
levels(im_cosmx_ebv_b_cells_2_cds@clusters@listData$UMAP$partitions)

im_cosmx_ebv_b_cells_2_cds <- learn_graph(im_cosmx_ebv_b_cells_2_cds, use_partition = FALSE)

im_cosmx_ebv_b_cells_2_cds <- order_cells(im_cosmx_ebv_b_cells_2_cds, reduction_method = "UMAP")

plot_cells(im_cosmx_ebv_b_cells_2_cds, color_cells_by = "pseudotime", label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 2)

# - Subset out the EBV+ cells and Naive

ebv.features <- c("BALF2","BARF1","BcLF1","BCRF1","BGLF4",
                  "BHRF1","BLLF1","BNLF2a","BNLF2b","BNRF1",
                  "BRLF1","EBER1","EBER2","EBNA1","EBNA2","EBNA3A",
                  "EBNA3B","EBNA3C","EBNA-LP","RPMS1/A73","LMP1","LMP2A","LMP2A/B","BZLF1")

naive_lmp1_cds <- choose_graph_segments(im_cosmx_ebv_b_cells_2_cds, clear_cds=FALSE)
ebv.genes <- ebv.features[c(17:24)]
naive_lmp1_cds <- naive_lmp1_cds[rowData(naive_lmp1_cds)$gene_short_name %in% ebv.genes,
                                                         colData(naive_lmp1_cds)$ebv_positive %in% c("TRUE")]
naive_lmp1_cds <- order_cells(naive_lmp1_cds)
plot_genes_in_pseudotime(naive_lmp1_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.01, cell_size = 1)


gene_fits <- fit_models(naive_lmp1_cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs_terms <- fit_coefs %>% filter(term == "pseudotime")
fit_coefs_terms.reduced <- fit_coefs_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

write.csv(fit_coefs_terms.reduced,".../naive_to_lmp1_deg.csv", row.names = F)





naive_lytic_cds <- choose_graph_segments(im_cosmx_ebv_b_cells_2_cds, clear_cds=FALSE)
ebv.genes <- ebv.features[c(17:24)]
naive_lytic_cds <- naive_lytic_cds[rowData(naive_lytic_cds)$gene_short_name %in% ebv.genes,
                                 colData(naive_lytic_cds)$ebv_positive %in% c("TRUE")]
naive_lytic_cds <- order_cells(naive_lytic_cds)
plot_genes_in_pseudotime(naive_lytic_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.01, cell_size = 1)


gene_fits <- fit_models(naive_lytic_cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs_terms <- fit_coefs %>% filter(term == "pseudotime")
fit_coefs_terms.reduced <- fit_coefs_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

write.csv(fit_coefs_terms.reduced,".../naive_to_lytic_deg.csv", row.names = F)




lytic_lmp1_cds <- choose_graph_segments(im_cosmx_ebv_b_cells_2_cds, clear_cds=FALSE)

gene_fits <- fit_models(lytic_lmp1_cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs_terms <- fit_coefs %>% filter(term == "pseudotime")
fit_coefs_terms.reduced <- fit_coefs_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

write.csv(fit_coefs_terms.reduced,".../lytic_to_lmp1_deg.csv", row.names = F)







# Plot gene expression as a function of pseudotime

ebv.genes <- ebv.features[c(17:24)]
im_cosmx_ebv_b_cells_2_cds <- im_cosmx_ebv_b_cells_2_cds[rowData(im_cosmx_ebv_b_cells_2_cds)$gene_short_name %in% ebv.genes,
                                        colData(im_cosmx_ebv_b_cells_2_cds)$ebv_positive %in% c("TRUE")]
im_cosmx_ebv_b_cells_2_cds <- order_cells(im_cosmx_ebv_b_cells_2_cds)

plot_genes_in_pseudotime(im_cosmx_ebv_b_cells_2_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.01, cell_size = 1)

saveRDS(im_cosmx_ebv_b_cells_2_cds, ".../b_cell_monocle.rds")



