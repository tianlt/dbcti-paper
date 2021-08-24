library(monocle)
library(TSCAN)
library(SCORPIUS)
library(SLICER)
library(igraph)
library(tidyverse)
source('Data/cell_trajectory/rcode/package_function.R')
source('Data/cell_trajectory/rcode/plot_tools_function_upload.R')
path = '/home/tialan/Data/cell_trajectory/trajectory_inference/plot/package_comparison/'
#hesc(use) normalized#########################
hescmt<-readRDS('/home/tialan/Data/cell_trajectory/dataset/hescmt.rds')
hescmt<-hescmt[, 214:460]
hescsmsinfo_cycle_index<-readRDS('/home/tialan/Data/cell_trajectory/dataset/hescsmsinfo_cycle_index.rds')
hescsmsinfo_cycle_index<-hescsmsinfo_cycle_index[214:460]
cellcycle<-c('CCNB1','CDK4','GMNN',
             'CCNB2',
             'TOP2A','MKI67',
             'CDKN1A','CDKN1B','CDKN1C',
             'CDK1','UBE2C')

#1. monocle###############
hesc_gene <- as.data.frame(hescmt[, 1])
rownames(hesc_gene) <- rownames(hescmt)
colnames(hesc_gene) <- 'gene_short_name'
hesc_cell <- as.data.frame(data.frame(as.numeric(hescsmsinfo_cycle_index)))
rownames(hesc_cell) <- colnames(hescmt)

hesc_monocle <- infer_monocle(data = hescmt,
                   gene_annotation =  hesc_gene,
                   cell_annotation = hesc_cell,
                   gene_to_use = cellcycle)


plot_cell_trajectory(hesc_monocle)
plot_cell_trajectory(hesc_monocle, color_by = "Pseudotime")
#eval
hesc_monocle@phenoData@data$index = factor(t(hescsmsinfo_cycle_index))
plot_cell_trajectory(hesc_monocle, color_by = "index")

saveRDS(hesc_monocle, 'Data/cell_trajectory/dataset/hesc_monocle.rds')

#2. tscan###############
hesc_tscan_p <- preprocess(hescmt)
hesc_tscan_c <- exprmclust(hesc_tscan_p)
hesc_tscan_o <- TSCANorder(hesc_tscan_c)

plotmclust(hesc_tscan_c)
#eval
plot_index_tscan(hesc_tscan_c, index = factor(t(hescsmsinfo_cycle_index)))

saveRDS(hesc_tscan_c, 'Data/cell_trajectory/dataset/hesc_tscan_c.rds')
saveRDS(hesc_tscan_o, 'Data/cell_trajectory/dataset/hesc_tscan_o.rds')

#3. SCORPIUS###############
hesc_group = hescsmsinfo_cycle_index
hesc_group = as.factor(rep(1, length(hesc_group)))
hesc_sco_s <- reduce_dimensionality(t(as.matrix(hescmt)), "spearman")
hesc_sco_t <- infer_trajectory(hesc_sco_s)

#group same
draw_trajectory_plot(hesc_sco_s, hesc_group, hesc_sco_t$path, contour = TRUE)
#eval
draw_trajectory_plot(hesc_sco_s, factor(t(hescsmsinfo_cycle_index)), hesc_sco_t$path, contour = TRUE)

saveRDS(hesc_sco_s, 'Data/cell_trajectory/dataset/hesc_sco_s.rds')
saveRDS(hesc_sco_t, 'Data/cell_trajectory/dataset/hesc_sco_t.rds')

#young (use)#############
stem_mouse_C57BL6_data_young<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young.rds')
stem_mouse_C57BL6_index_young<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_index_young.rds')

#1. monocle###############
stem_mouse_C57BL6_data_young_gene <- as.data.frame(stem_mouse_C57BL6_data_young[, 1])
rownames(stem_mouse_C57BL6_data_young_gene) <- rownames(stem_mouse_C57BL6_data_young)
colnames(stem_mouse_C57BL6_data_young_gene) <- 'gene_short_name'
stem_mouse_C57BL6_data_young_cell <- data.frame(stem_mouse_C57BL6_index_young)
rownames(stem_mouse_C57BL6_data_young_cell) <- colnames(stem_mouse_C57BL6_data_young)

stem_mouse_C57BL6_data_young_monocle <- infer_monocle(data = stem_mouse_C57BL6_data_young, gene_annotation =  stem_mouse_C57BL6_data_young_gene, cell_annotation = stem_mouse_C57BL6_data_young_cell, gene_to_use = rownames(stem_mouse_C57BL6_data_young))
plot_cell_trajectory(stem_mouse_C57BL6_young_monocle)
plot_cell_trajectory(stem_mouse_C57BL6_young_monocle, color_by = "Pseudotime")
#eval
stem_mouse_C57BL6_data_young_monocle@phenoData@data$index = factor(t(stem_mouse_C57BL6_index_young))
plot_cell_trajectory(stem_mouse_C57BL6_data_young_monocle, color_by = "index")

saveRDS(stem_mouse_C57BL6_data_young_monocle, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young_monocle.rds')

#2. tscan###############
stem_mouse_C57BL6_data_young_tscan_p <- preprocess(stem_mouse_C57BL6_data_young)
stem_mouse_C57BL6_data_young_tscan_c <- exprmclust(stem_mouse_C57BL6_data_young_tscan_p)
stem_mouse_C57BL6_data_young_tscan_o <- TSCANorder(stem_mouse_C57BL6_data_young_tscan_c)

plotmclust(stem_mouse_C57BL6_data_young_tscan_c)
#eval
plot_index_tscan(stem_mouse_C57BL6_data_young_tscan_c, index = factor(t(stem_mouse_C57BL6_index_young)))

saveRDS(stem_mouse_C57BL6_data_young_tscan_c, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young_tscan_c.rds')
saveRDS(stem_mouse_C57BL6_data_young_tscan_o, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young_tscan_o.rds')

#3. SCORPIUS###############
stem_mouse_C57BL6_data_young_group = as.factor(rep(1, dim(stem_mouse_C57BL6_data_young)[2]))

stem_mouse_C57BL6_data_young_s <- reduce_dimensionality(t(as.matrix(stem_mouse_C57BL6_data_young)), "spearman")
stem_mouse_C57BL6_data_young_t <- infer_trajectory(stem_mouse_C57BL6_data_young_s)

#group same
draw_trajectory_plot(stem_mouse_C57BL6_data_young_s, stem_mouse_C57BL6_data_young_group, stem_mouse_C57BL6_data_young_t$path, contour = TRUE)

#eval
draw_trajectory_plot(stem_mouse_C57BL6_data_young_s, factor(t(stem_mouse_C57BL6_index_young)), stem_mouse_C57BL6_data_young_t$path, contour = TRUE)

saveRDS(stem_mouse_C57BL6_data_young_s, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young_s.rds')
saveRDS(stem_mouse_C57BL6_data_young_t, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young_t.rds')

#old (use)###########
stem_mouse_C57BL6_data_old<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old.rds')
stem_mouse_C57BL6_index_old<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_index_old.rds')

#1. monocle###############
stem_mouse_C57BL6_data_old_gene <- as.data.frame(stem_mouse_C57BL6_data_old[, 1])
rownames(stem_mouse_C57BL6_data_old_gene) <- rownames(stem_mouse_C57BL6_data_old)
colnames(stem_mouse_C57BL6_data_old_gene) <- 'gene_short_name'
stem_mouse_C57BL6_data_old_cell <- data.frame(stem_mouse_C57BL6_index_old)
rownames(stem_mouse_C57BL6_data_old_cell) <- colnames(stem_mouse_C57BL6_data_old)

stem_mouse_C57BL6_old_monocle <- infer_monocle(data = stem_mouse_C57BL6_data_old, gene_annotation =  stem_mouse_C57BL6_data_old_gene, cell_annotation = stem_mouse_C57BL6_data_old_cell, gene_to_use = rownames(stem_mouse_C57BL6_data_old))

#evaluate
plot_cell_trajectory(stem_mouse_C57BL6_old_monocle)
plot_cell_trajectory(stem_mouse_C57BL6_old_monocle, color_by = "Pseudotime")


#eval
stem_mouse_C57BL6_old_monocle@phenoData@data$index = factor(t(stem_mouse_C57BL6_index_old))
plot_cell_trajectory(stem_mouse_C57BL6_old_monocle, color_by = "index")

saveRDS(stem_mouse_C57BL6_old_monocle, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_old_monocle.rds')

#2. tscan###############
stem_mouse_C57BL6_data_old_tscan_p <- preprocess(stem_mouse_C57BL6_data_old)
stem_mouse_C57BL6_data_old_tscan_c <- exprmclust(stem_mouse_C57BL6_data_old_tscan_p)
stem_mouse_C57BL6_data_old_tscan_o <- TSCANorder(stem_mouse_C57BL6_data_old_tscan_c)

plotmclust(stem_mouse_C57BL6_data_old_tscan_c)

#eval
plot_index_tscan(stem_mouse_C57BL6_data_old_tscan_c, index = factor(t(stem_mouse_C57BL6_index_old)))

saveRDS(stem_mouse_C57BL6_data_old_tscan_c, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old_tscan_c.rds')
saveRDS(stem_mouse_C57BL6_data_old_tscan_o, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old_tscan_o.rds')

#3. SCORPIUS###############
stem_mouse_C57BL6_data_old_group = as.factor(rep(1, dim(stem_mouse_C57BL6_data_old)[2]))

stem_mouse_C57BL6_data_old_s <- reduce_dimensionality(t(as.matrix(stem_mouse_C57BL6_data_old)), "spearman")
stem_mouse_C57BL6_data_old_t <- infer_trajectory(stem_mouse_C57BL6_data_old_s)

#group same
draw_trajectory_plot(stem_mouse_C57BL6_data_old_s, stem_mouse_C57BL6_data_old_group, stem_mouse_C57BL6_data_old_t$path, contour = TRUE)

#eval
draw_trajectory_plot(stem_mouse_C57BL6_data_old_s, factor(t(stem_mouse_C57BL6_index_old)), stem_mouse_C57BL6_data_old_t$path, contour = TRUE)

saveRDS(stem_mouse_C57BL6_data_old_s, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old_s.rds')
saveRDS(stem_mouse_C57BL6_data_old_t, 'Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old_t.rds')

#batch1(use)###############
camp1_data_batch1 <- readRDS('/home/tialan/Data/cell_trajectory/dataset/camp1_data_batch1.rds')
camp1_index_batch1 <- readRDS('/home/tialan/Data/cell_trajectory/dataset/camp1_index_batch1.rds')
camp1_batch1_selected_gene <- readRDS('Data/cell_trajectory/dataset/camp1_batch1_selected_gene.rds')
camp1_index_batch1 <- camp1_index_batch1[2,]


#1. monocle###############
camp1_data_batch1_gene <- as.data.frame(camp1_data_batch1[, 1])
rownames(camp1_data_batch1_gene) <- rownames(camp1_data_batch1_gene)
colnames(camp1_data_batch1_gene) <- 'gene_short_name'
camp1_data_batch1_cell <- data.frame(t(data.frame(camp1_index_batch1)))
rownames(camp1_data_batch1_cell) <- colnames(camp1_data_batch1)

camp1_data_batch1_monocle <- infer_monocle(data = camp1_data_batch1, gene_annotation =  camp1_data_batch1_gene, cell_annotation = camp1_data_batch1_cell, gene_to_use = rownames(camp1_data_batch1))

#evaluate
plot_cell_trajectory(camp1_data_batch1_monocle)
plot_cell_trajectory(camp1_data_batch1_monocle, color_by = "Pseudotime")

#eval
camp1_data_batch1_monocle@phenoData@data$index = factor(t(camp1_index_batch1))
plot_cell_trajectory(camp1_data_batch1_monocle, color_by = "index")

saveRDS(camp1_data_batch1_monocle, 'Data/cell_trajectory/dataset/camp1_data_batch1_monocle.rds')

#2. tscan###############
camp1_data_batch1_tscan_p <- preprocess(camp1_data_batch1)
camp1_data_batch1_tscan_c <- exprmclust(camp1_data_batch1_tscan_p)
camp1_data_batch1_tscan_o <- TSCANorder(camp1_data_batch1_tscan_c)

plotmclust(camp1_data_batch1_tscan_c)

#eval
plot_index_tscan(camp1_data_batch1_tscan_c, index = factor(t(camp1_index_batch1)))

saveRDS(camp1_data_batch1_tscan_c, 'Data/cell_trajectory/dataset/camp1_data_batch1_tscan_c.rds')
saveRDS(camp1_data_batch1_tscan_o, 'Data/cell_trajectory/dataset/camp1_data_batch1_tscan_o.rds')

#3. SCORPIUS###############
camp1_data_batch1_group = as.factor(rep(1, dim(camp1_data_batch1)[2]))

camp1_data_batch1_s <- reduce_dimensionality(t(as.matrix(camp1_data_batch1)), "spearman")
camp1_data_batch1_t <- infer_trajectory(camp1_data_batch1_s)

#group same
draw_trajectory_plot(camp1_data_batch1_s, camp1_data_batch1_group, camp1_data_batch1_t$path, contour = TRUE)

#eval
draw_trajectory_plot(camp1_data_batch1_s, factor(t(camp1_index_batch1)), camp1_data_batch1_t$path, contour = TRUE)

saveRDS(camp1_data_batch1_s, 'Data/cell_trajectory/dataset/camp1_data_batch1_s.rds')
saveRDS(camp1_data_batch1_t, 'Data/cell_trajectory/dataset/camp1_data_batch1_t.rds')


#nestorowa_data_batch1_type_filter no stem(use1)####################
nestorowa_data_batch1_type_filter_no_stem <- readRDS('/home/tialan/Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem.rds')
nestorowa_index_batch1_type_filter_no_stem <- readRDS('/home/tialan/Data/cell_trajectory/dataset/nestorowa_index_batch1_type_filter_no_stem.rds')
nestorowa_batch1_type_filter_no_stem_selected_gene <- readRDS('Data/cell_trajectory/dataset/nestorowa_batch1_type_filter_no_stem_selected_gene.rds')
nestorowa_index_batch1_type_filter_no_stem <- nestorowa_index_batch1_type_filter_no_stem[2,]
nestorowa_data_batch1_type_filter_no_stem <- nestorowa_data_batch1_type_filter_no_stem[nestorowa_batch1_type_filter_no_stem_selected_gene, ]

#1. monocle###############
nestorowa_data_batch1_type_filter_no_stem_gene <- as.data.frame(nestorowa_data_batch1_type_filter_no_stem[, 1])
rownames(nestorowa_data_batch1_type_filter_no_stem_gene) <- rownames(nestorowa_data_batch1_type_filter_no_stem_gene)
colnames(nestorowa_data_batch1_type_filter_no_stem_gene) <- 'gene_short_name'

nestorowa_data_batch1_type_filter_no_stem_cell <- data.frame(nestorowa_index_batch1_type_filter_no_stem)
rownames(nestorowa_data_batch1_type_filter_no_stem_cell) <- colnames(nestorowa_data_batch1_type_filter_no_stem)

nestorowa_data_batch1_type_filter_no_stem_monocle <- infer_monocle(data = nestorowa_data_batch1_type_filter_no_stem, gene_annotation =  nestorowa_data_batch1_type_filter_no_stem_gene, cell_annotation = nestorowa_data_batch1_type_filter_no_stem_cell, gene_to_use = rownames(nestorowa_data_batch1_type_filter_no_stem))

#evaluate
plot_cell_trajectory(nestorowa_data_batch1_type_filter_no_stem_monocle)
plot_cell_trajectory(nestorowa_data_batch1_type_filter_no_stem_monocle, color_by = "Pseudotime")

#eval
nestorowa_data_batch1_type_filter_no_stem_monocle@phenoData@data$index = factor(t(nestorowa_index_batch1_type_filter_no_stem))
plot_cell_trajectory(nestorowa_data_batch1_type_filter_no_stem_monocle, color_by = "index")

saveRDS(nestorowa_data_batch1_type_filter_no_stem_monocle, 'Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem_monocle.rds')

#2. tscan###############
nestorowa_data_batch1_type_filter_no_stem_p <- preprocess(nestorowa_data_batch1_type_filter_no_stem, minexpr_percent = 0.2)
nestorowa_data_batch1_type_filter_no_stem_c <- exprmclust(nestorowa_data_batch1_type_filter_no_stem_p)
nestorowa_data_batch1_type_filter_no_stem_o <- TSCANorder(nestorowa_data_batch1_type_filter_no_stem_c)

plotmclust(nestorowa_data_batch1_type_filter_no_stem_c)

#eval
plot_index_tscan(nestorowa_data_batch1_type_filter_no_stem_c, index = factor(t(nestorowa_index_batch1_type_filter_no_stem)))

saveRDS(nestorowa_data_batch1_type_filter_no_stem_c, 'Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem_c.rds')
saveRDS(nestorowa_data_batch1_type_filter_no_stem_o, 'Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem_o.rds')

#3. SCORPIUS###############
nestorowa_data_batch1_type_filter_no_stem_group = as.factor(rep(1, dim(nestorowa_data_batch1_type_filter_no_stem)[2]))

nestorowa_data_batch1_type_filter_no_stem_s <- reduce_dimensionality(t(as.matrix(nestorowa_data_batch1_type_filter_no_stem)), "spearman")
nestorowa_data_batch1_type_filter_no_stem_t <- infer_trajectory(nestorowa_data_batch1_type_filter_no_stem_s)

#group same
draw_trajectory_plot(nestorowa_data_batch1_type_filter_no_stem_s, nestorowa_data_batch1_type_filter_no_stem_group, nestorowa_data_batch1_type_filter_no_stem_t$path, contour = TRUE)

#eval
draw_trajectory_plot(nestorowa_data_batch1_type_filter_no_stem_s, factor(t(nestorowa_index_batch1_type_filter_no_stem)), nestorowa_data_batch1_type_filter_no_stem_t$path, contour = TRUE)

saveRDS(nestorowa_data_batch1_type_filter_no_stem_s, 'Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem_s.rds')
saveRDS(nestorowa_data_batch1_type_filter_no_stem_t, 'Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem_t.rds')

#yan human embryo(gse36552) (hemberg lab github) (normalized) (use)#########################
yan_data <- readRDS('/home/tialan/Data/cell_trajectory/dataset/yan_data.rds')
yan_index <- readRDS('/home/tialan/Data/cell_trajectory/dataset/yan_index.rds')
yan_selected_gene <- readRDS('Data/cell_trajectory/dataset/yan_selected_gene.rds')
yan_index <- yan_index[1, ]

#1. monocle###############
yan_data_gene <- as.data.frame(yan_data[, 1])
rownames(yan_data_gene) <- rownames(yan_data)
colnames(yan_data_gene) <- 'gene_short_name'
yan_data_cell <- data.frame(t(data.frame(yan_index)))
rownames(yan_data_cell) <- colnames(yan_data)

yan_data_monocle <- infer_monocle(data = yan_data, gene_annotation =  yan_data_gene, cell_annotation = yan_data_cell, gene_to_use = rownames(yan_data))

#evaluate
plot_cell_trajectory(yan_data_monocle)
plot_cell_trajectory(yan_data_monocle, color_by = "Pseudotime")

#eval
yan_data_monocle@phenoData@data$index = factor(t(yan_index))
plot_cell_trajectory(yan_data_monocle, color_by = "index")

saveRDS(yan_data_monocle, 'Data/cell_trajectory/dataset/yan_data_monocle.rds')
#2. tscan###############
yan_data_p <- preprocess(yan_data)
yan_data_c <- exprmclust(yan_data_p)
yan_data_o <- TSCANorder(yan_data_c)

plotmclust(yan_data_c)
#eval
plot_index_tscan(yan_data_c, index = factor(t(yan_index)))

saveRDS(yan_data_c, 'Data/cell_trajectory/dataset/yan_data_c.rds')
saveRDS(yan_data_o, 'Data/cell_trajectory/dataset/yan_data_o.rds')

#3. SCORPIUS###############
yan_data_group = as.factor(rep(1, dim(yan_data)[2]))

yan_data_s <- reduce_dimensionality(t(as.matrix(yan_data)), "spearman")
yan_data_t <- infer_trajectory(yan_data_s)

#group same
draw_trajectory_plot(yan_data_s, yan_data_group, yan_data_t$path, contour = TRUE)

#eval
draw_trajectory_plot(yan_data_s, factor(t(yan_index)), yan_data_t$path, contour = TRUE)

saveRDS(yan_data_s, 'Data/cell_trajectory/dataset/yan_data_c.rds')
saveRDS(yan_data_t, 'Data/cell_trajectory/dataset/yan_data_t.rds')


#plot and measurement########################
#p  hesc(use) normalized#########################
#p1. monocle###############
pdf(file = paste0(path, 'hesc_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(hesc_monocle)
dev.off()
pdf(file = paste0(path, 'hesc_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(hesc_monocle, color_by = "Pseudotime")
dev.off()
#eval
pdf(file = paste0(path, 'hesc_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(hesc_monocle, color_by = "index")
dev.off()

#p2. tscan###############
pdf(file = paste0(path, 'hesc_tscan.pdf'),width = 6,height = 6)
tscan_plot(hesc_tscan_c, show_cell_names = FALSE)
dev.off()
#eval
pdf(file = paste0(path, 'hesc_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(hesc_tscan_c, index = factor(t(hescsmsinfo_cycle_index)))
dev.off()
#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'hesc_sco.pdf'),width = 6,height = 6)
sco_plot(hesc_sco_s, hesc_group, hesc_sco_t$path, contour = TRUE)
dev.off()
#eval
pdf(file = paste0(path, 'hesc_sco_index.pdf'),width = 6,height = 6)
sco_plot(hesc_sco_s, factor(t(hescsmsinfo_cycle_index)), hesc_sco_t$path, contour = TRUE)
dev.off()

#p  young (use)#############
#p1. monocle###############
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_young_monocle)
dev.off()
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_young_monocle, color_by = "Pseudotime")
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_data_young_monocle, color_by = "index")
dev.off()
#p2. tscan###############
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_tscan.pdf'),width = 6,height = 6)
tscan_plot(stem_mouse_C57BL6_data_young_tscan_c, show_cell_names = FALSE)
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(stem_mouse_C57BL6_data_young_tscan_c, index = factor(t(stem_mouse_C57BL6_index_young)))
dev.off()
#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_sco.pdf'),width = 6,height = 6)
sco_plot(stem_mouse_C57BL6_data_young_s, stem_mouse_C57BL6_data_young_group, stem_mouse_C57BL6_data_young_t$path, contour = TRUE)
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_sco_index.pdf'),width = 6,height = 6)
sco_plot(stem_mouse_C57BL6_data_young_s, factor(t(stem_mouse_C57BL6_index_young)), stem_mouse_C57BL6_data_young_t$path, contour = TRUE)
dev.off()
#p  old (use)###########

#p1. monocle###############

#evaluate
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_monocle.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_old_monocle)
dev.off()
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_old_monocle, color_by = "Pseudotime")
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(stem_mouse_C57BL6_old_monocle, color_by = "index")
dev.off()
#p2. tscan###############
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_tscan.pdf'),width = 6,height = 6)
tscan_plot(stem_mouse_C57BL6_data_old_tscan_c, show_cell_names = FALSE)
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(stem_mouse_C57BL6_data_old_tscan_c, index = factor(t(stem_mouse_C57BL6_index_old)))
dev.off()
#p3. SCORPIUS###############

#group same
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_sco.pdf'),width = 6,height = 6)
sco_plot(stem_mouse_C57BL6_data_old_s, stem_mouse_C57BL6_data_old_group, stem_mouse_C57BL6_data_old_t$path, contour = TRUE)
dev.off()
#eval
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_sco_index.pdf'),width = 6,height = 6)
sco_plot(stem_mouse_C57BL6_data_old_s, factor(t(stem_mouse_C57BL6_index_old)), stem_mouse_C57BL6_data_old_t$path, contour = TRUE)
dev.off()
#p  batch1(use)###############
#p1. monocle###############

#evaluate
pdf(file = paste0(path, 'camp1_data_batch1_monocle.pdf'),width = 6,height = 6)
monocle_plot(camp1_data_batch1_monocle)
dev.off()

pdf(file = paste0(path, 'camp1_data_batch1_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(camp1_data_batch1_monocle, color_by = "Pseudotime")
dev.off()

#eval
pdf(file = paste0(path, 'camp1_data_batch1_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(camp1_data_batch1_monocle, color_by = "index")
dev.off()

#p2. tscan###############
pdf(file = paste0(path, 'camp1_data_batch1_tscan.pdf'),width = 6,height = 6)
tscan_plot(camp1_data_batch1_tscan_c, show_cell_names = FALSE)
dev.off()

#eval
pdf(file = paste0(path, 'camp1_data_batch1_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(camp1_data_batch1_tscan_c, index = factor(t(camp1_index_batch1)))
dev.off()

#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'camp1_data_batch1_sco.pdf'),width = 6,height = 6)
sco_plot(camp1_data_batch1_s, camp1_data_batch1_group, camp1_data_batch1_t$path, contour = TRUE)
dev.off()

#eval
pdf(file = paste0(path, 'camp1_data_batch1_sco_index.pdf'),width = 6,height = 6)
sco_plot(camp1_data_batch1_s, factor(t(camp1_index_batch1)), camp1_data_batch1_t$path, contour = TRUE)
dev.off()

#p  nestorowa_data_batch1_type_filter(use2)####################
#p1. monocle###############
#evaluate
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_monocle.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_monocle)
dev.off()

pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_monocle, color_by = "Pseudotime")
dev.off()

#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_monocle, color_by = "index")
dev.off()
#p2. tscan###############
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_tscan.pdf'),width = 6,height = 6)
tscan_plot(nestorowa_data_batch1_type_filter_c, show_cell_names = FALSE)
dev.off()
#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(nestorowa_data_batch1_type_filter_c, index = factor(t(nestorowa_index_batch1_type_filter)))
dev.off()

#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_sco.pdf'),width = 6,height = 6)
sco_plot(nestorowa_data_batch1_type_filter_s, nestorowa_data_batch1_type_filter_group, nestorowa_data_batch1_type_filter_t$path, contour = TRUE)
dev.off()

#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_sco_index.pdf'),width = 6,height = 6)
sco_plot(nestorowa_data_batch1_type_filter_s, factor(t(nestorowa_index_batch1_type_filter)), nestorowa_data_batch1_type_filter_t$path, contour = TRUE)
dev.off()

#p  nestorowa_data_batch1_type_filter no stem(use1)####################
#p1. monocle###############
#evaluate
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_monocle.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_no_stem_monocle)
dev.off()

pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_no_stem_monocle, color_by = "Pseudotime")
dev.off()
#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(nestorowa_data_batch1_type_filter_no_stem_monocle, color_by = "index")
dev.off()

#p2. tscan###############
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_tscan.pdf'),width = 6,height = 6)
tscan_plot(nestorowa_data_batch1_type_filter_no_stem_c, show_cell_names = FALSE)
dev.off()

#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(nestorowa_data_batch1_type_filter_no_stem_c, index = factor(t(nestorowa_index_batch1_type_filter_no_stem)))
dev.off()

#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_sco.pdf'),width = 6,height = 6)
sco_plot(nestorowa_data_batch1_type_filter_no_stem_s, nestorowa_data_batch1_type_filter_no_stem_group, nestorowa_data_batch1_type_filter_no_stem_t$path, contour = TRUE)
dev.off()

#eval
pdf(file = paste0(path, 'nestorowa_data_batch1_type_filter_no_stem_sco_index.pdf'),width = 6,height = 6)
sco_plot(nestorowa_data_batch1_type_filter_no_stem_s, factor(t(nestorowa_index_batch1_type_filter_no_stem)), nestorowa_data_batch1_type_filter_no_stem_t$path, contour = TRUE)
dev.off()


#p  yan human embryo(gse36552) (hemberg lab github) (normalized) (use)#########################
#p1. monocle###############
#evaluate
pdf(file = paste0(path, 'yan_data_monocle.pdf'),width = 6,height = 6)
monocle_plot(yan_data_monocle)
dev.off()

pdf(file = paste0(path, 'yan_data_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(yan_data_monocle, color_by = "Pseudotime")
dev.off()
#eval
pdf(file = paste0(path, 'yan_data_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(yan_data_monocle, color_by = "index")
dev.off()

#p2. tscan###############
pdf(file = paste0(path, 'yan_data_tscan.pdf'),width = 6,height = 6)
tscan_plot(yan_data_c, show_cell_names = FALSE)
dev.off()
#eval
pdf(file = paste0(path, 'yan_data_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(yan_data_c, index = factor(t(yan_index)))
dev.off()
        
#p3. SCORPIUS###############
#group same
pdf(file = paste0(path, 'yan_data_sco.pdf'),width = 6,height = 6)
sco_plot(yan_data_s, yan_data_group, yan_data_t$path, contour = TRUE)
dev.off()

#eval
pdf(file = paste0(path, 'yan_data_sco_index.pdf'),width = 6,height = 6)
sco_plot(yan_data_s, factor(t(yan_index)), yan_data_t$path, contour = TRUE)
dev.off()