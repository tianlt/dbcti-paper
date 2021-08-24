source(file = 'Data/cell_trajectory/rcode/package_function.R')
source('Data/cell_trajectory/rcode/plot_tools_function_upload.R')
library(SCORPIUS)
library(monocle)
library(TSCAN)
library(igraph)
library(tidyverse)

#bifurcating_cycle##########################
bifurcating_cycle_data <- readRDS(file = 'Data/cell_trajectory/dataset/bifurcating_cycle_data.rds')
bifurcating_cycle_index <- readRDS(file = 'Data/cell_trajectory/dataset/bifurcating_cycle_index.rds')
bifurcating_cycle_loc_dimred <- readRDS(file = 'Data/cell_trajectory/dataset/bifurcating_cycle_loc_dimred.rds')

#1. monocle###############
bifurcating_cycle_gene <- as.data.frame(bifurcating_cycle_data[, 1])
rownames(bifurcating_cycle_gene) <- rownames(bifurcating_cycle_data)
colnames(bifurcating_cycle_gene) <- 'gene_short_name'
bifurcating_cycle_cell <- as.data.frame(bifurcating_cycle_index)
rownames(bifurcating_cycle_cell) <- colnames(bifurcating_cycle_data)


bifurcating_cycle_monocle <- infer_monocle(data = bifurcating_cycle_data,
                              gene_annotation =  bifurcating_cycle_gene,
                              cell_annotation = bifurcating_cycle_cell,
                              gene_to_use = rownames(bifurcating_cycle_data))

bifurcating_cycle_monocle@phenoData@data$index = factor(t(bifurcating_cycle_index[, 2]))

plot_cell_trajectory(bifurcating_cycle_monocle)
plot_cell_trajectory(bifurcating_cycle_monocle, color_by = "Pseudotime")

saveRDS(bifurcating_cycle_monocle, 'Data/cell_trajectory/dataset/bifurcating_cycle_monocle.rds')


#2. tscan###############
bifurcating_cycle_tscan_p <- preprocess(bifurcating_cycle_data, minexpr_percent = 0.01)
bifurcating_cycle_tscan_c <- exprmclust(bifurcating_cycle_tscan_p)
bifurcating_cycle_tscan_o <- TSCANorder(bifurcating_cycle_tscan_c)

#eval
plotmclust(bifurcating_cycle_tscan_c)

saveRDS(bifurcating_cycle_tscan_c, 'Data/cell_trajectory/dataset/bifurcating_cycle_tscan_c.rds')
saveRDS(bifurcating_cycle_tscan_o, 'Data/cell_trajectory/dataset/bifurcating_cycle_tscan_o.rds')

#3. SCORPIUS###############
bifurcating_cycle_group = bifurcating_cycle_index[, 2]
bifurcating_cycle_group = as.factor(rep(1, dim(bifurcating_cycle_group)[1]))
bifurcating_cycle_sco_s <- reduce_dimensionality(t(as.matrix(bifurcating_cycle_data)), "spearman")
bifurcating_cycle_sco_t <- SCORPIUS::infer_trajectory(bifurcating_cycle_sco_s)

#group same eval
draw_trajectory_plot(bifurcating_cycle_sco_s, bifurcating_cycle_group, bifurcating_cycle_sco_t$path, contour = TRUE)
draw_trajectory_plot(bifurcating_cycle_sco_s, factor(t(bifurcating_cycle_index[, 2])), bifurcating_cycle_sco_t$path, contour = TRUE)

saveRDS(bifurcating_cycle_sco_s, 'Data/cell_trajectory/dataset/bifurcating_cycle_sco_s.rds')
saveRDS(bifurcating_cycle_sco_t, 'Data/cell_trajectory/dataset/bifurcating_cycle_sco_t.rds')


#5. tool use###########################
bifurcating_cycle_data <- as.data.frame(bifurcating_cycle_data)
test_bifurcating_cycle <- create_object(bifurcating_cycle_data, normalized = TRUE) 
test_bifurcating_cycle <- filter_data(test_bifurcating_cycle, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_bifurcating_cycle <- select_var_feature(test_bifurcating_cycle, use_normalized_data = TRUE, n = 1000)
test_bifurcating_cycle <- tsneplot(test_bifurcating_cycle, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_bifurcating_cycle <- contour_plot(test_bifurcating_cycle)
test_bifurcating_cycle <- distribution_estimation(test_bifurcating_cycle, ndraw = 1000, expansion = 2.5, ... = 1,2,3,4,5) 
test_bifurcating_cycle <- point_possibility(test_bifurcating_cycle, r = 5)
test_bifurcating_cycle <- connect_cluster(test_bifurcating_cycle, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.05)
test_bifurcating_cycle <- infer_trajectory(test_bifurcating_cycle, iter_n =15)
test_bifurcating_cycle <- calculate_pseudotime(test_bifurcating_cycle, start_state_name = c('1'))
test_bifurcating_cycle <- plot_trajectory(test_bifurcating_cycle)
plot(test_bifurcating_cycle@trajectory_plot$plot)

#evaluation
eval_bifurcating_cycle <- plot_index(index = factor(t(bifurcating_cycle_index[, 2])), trajectory = test_bifurcating_cycle@trajectory, connection_matrix = test_bifurcating_cycle@connect_cluster$cluster_connection)
plot(eval_bifurcating_cycle$plot)
  
saveRDS(test_bifurcating_cycle, 'Data/cell_trajectory/dataset/test_bifurcating_cycle.rds')


#cycle##########################
cycle_data <- readRDS(file = 'Data/cell_trajectory/dataset/cycle_data.rds')
cycle_index <- readRDS(file = 'Data/cell_trajectory/dataset/cycle_index.rds')
cycle_loc_dimred <- readRDS(file = 'Data/cell_trajectory/dataset/cycle_loc_dimred.rds')

#1. monocle###############
cycle_gene <- as.data.frame(cycle_data[, 1])
rownames(cycle_gene) <- rownames(cycle_data)
colnames(cycle_gene) <- 'gene_short_name'
cycle_cell <- as.data.frame(cycle_index)
rownames(cycle_cell) <- colnames(cycle_data)


cycle_monocle <- infer_monocle(data = cycle_data,
                                           gene_annotation =  cycle_gene,
                                           cell_annotation = cycle_cell,
                                           gene_to_use = rownames(cycle_data))

cycle_monocle@phenoData@data$index = factor(t(cycle_index[, 2]))

plot_cell_trajectory(cycle_monocle)
plot_cell_trajectory(cycle_monocle, color_by = "Pseudotime")

saveRDS(cycle_monocle, 'Data/cell_trajectory/dataset/cycle_monocle.rds')

#2. tscan###############
cycle_tscan_p <- preprocess(cycle_data, minexpr_percent = 0.01)
cycle_tscan_c <- exprmclust(cycle_tscan_p)
cycle_tscan_o <- TSCANorder(cycle_tscan_c)

#eval
plotmclust(cycle_tscan_c)

saveRDS(cycle_tscan_c, 'Data/cell_trajectory/dataset/cycle_tscan_c.rds')
saveRDS(cycle_tscan_o, 'Data/cell_trajectory/dataset/cycle_tscan_o.rds')

#3. SCORPIUS###############
cycle_group = cycle_index[, 2]
cycle_group = as.factor(rep(1, dim(cycle_group)[1]))
cycle_sco_s <- reduce_dimensionality(t(as.matrix(cycle_data)), "spearman")
cycle_sco_t <- SCORPIUS::infer_trajectory(cycle_sco_s)

#group same eval
draw_trajectory_plot(cycle_sco_s, cycle_group, cycle_sco_t$path, contour = TRUE)
draw_trajectory_plot(cycle_sco_s, factor(t(cycle_index[, 2])), cycle_sco_t$path, contour = TRUE)

saveRDS(cycle_sco_s, 'Data/cell_trajectory/dataset/cycle_sco_s.rds')
saveRDS(cycle_sco_t, 'Data/cell_trajectory/dataset/cycle_sco_t.rds')

#5. tool (use)###########################
cycle_data <- as.data.frame(cycle_data)
test_cycle <- create_object(cycle_data, normalized = TRUE) 
test_cycle <- filter_data(test_cycle, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_cycle <- select_var_feature(test_cycle, use_normalized_data = TRUE, n = 1000)
test_cycle <- tsneplot(test_cycle, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_cycle <- contour_plot(test_cycle)
test_cycle <- distribution_estimation(test_cycle, ndraw = 1000, expansion = 2.5, ... = c(1,5), 4, 2, 3) 
test_cycle <- point_possibility(test_cycle, r = 5)
test_cycle <- connect_cluster(test_cycle, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_cycle <- infer_trajectory(test_cycle, iter_n =5)
test_cycle <- calculate_pseudotime(test_cycle, start_state_name = c('1'))
test_cycle <- plot_trajectory(test_cycle)
plot(test_cycle@trajectory_plot$plot)

#evaluation
eval_cycle <- plot_index(index = factor(t(cycle_index[, 2])), trajectory = test_cycle@trajectory, connection_matrix = test_cycle@connect_cluster$cluster_connection)
plot(eval_cycle$plot)

saveRDS(test_cycle, 'Data/cell_trajectory/dataset/test_cycle.rds')

#disconnected##########################
disconnected_data <- readRDS(file = 'Data/cell_trajectory/dataset/disconnected_data.rds')
disconnected_index <- readRDS(file = 'Data/cell_trajectory/dataset/disconnected_index.rds')
disconnected_loc_dimred <- readRDS(file = 'Data/cell_trajectory/dataset/disconnected_loc_dimred.rds')

#1. monocle###############
disconnected_gene <- as.data.frame(disconnected_data[, 1])
rownames(disconnected_gene) <- rownames(disconnected_data)
colnames(disconnected_gene) <- 'gene_short_name'
disconnected_cell <- as.data.frame(disconnected_index)
rownames(disconnected_cell) <- colnames(disconnected_data)


disconnected_monocle <- infer_monocle(data = disconnected_data,
                               gene_annotation =  disconnected_gene,
                               cell_annotation = disconnected_cell,
                               gene_to_use = rownames(disconnected_data))

disconnected_monocle@phenoData@data$index = factor(t(disconnected_index[, 2]))

plot_cell_trajectory(disconnected_monocle)
plot_cell_trajectory(disconnected_monocle, color_by = "Pseudotime")

saveRDS(disconnected_monocle, 'Data/cell_trajectory/dataset/disconnected_monocle.rds')

#2. tscan###############
disconnected_tscan_p <- preprocess(disconnected_data, minexpr_percent = 0.01)
disconnected_tscan_c <- exprmclust(disconnected_tscan_p)
disconnected_tscan_o <- TSCANorder(disconnected_tscan_c)

#eval
plotmclust(disconnected_tscan_c)

saveRDS(disconnected_tscan_c, 'Data/cell_trajectory/dataset/disconnected_tscan_c.rds')
saveRDS(disconnected_tscan_o, 'Data/cell_trajectory/dataset/disconnected_tscan_o.rds')

#3. SCORPIUS###############
disconnected_group = disconnected_index[, 2]
disconnected_group = as.factor(rep(1, dim(disconnected_group)[1]))
disconnected_sco_s <- reduce_dimensionality(t(as.matrix(disconnected_data)), "spearman")
disconnected_sco_t <- SCORPIUS::infer_trajectory(disconnected_sco_s)

#group same eval
draw_trajectory_plot(disconnected_sco_s, disconnected_group, disconnected_sco_t$path, contour = TRUE)
draw_trajectory_plot(disconnected_sco_s, factor(t(disconnected_index[, 2])), disconnected_sco_t$path, contour = TRUE)

saveRDS(disconnected_sco_s, 'Data/cell_trajectory/dataset/disconnected_sco_s.rds')
saveRDS(disconnected_sco_t, 'Data/cell_trajectory/dataset/disconnected_sco_t.rds')

#5. tool use###########################
disconnected_data <- as.data.frame(disconnected_data)
test_disconnected <- create_object(disconnected_data, normalized = TRUE) 
test_disconnected <- filter_data(test_disconnected, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_disconnected <- select_var_feature(test_disconnected, use_normalized_data = TRUE, n = 991)
test_disconnected <- tsneplot(test_disconnected, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_disconnected <- contour_plot(test_disconnected)
test_disconnected <- distribution_estimation(test_disconnected, ndraw = 1000, expansion = 2.5, ... = c(3,4), 6, 5, 1, 2) 
test_disconnected <- point_possibility(test_disconnected, r = 5)
test_disconnected <- connect_cluster(test_disconnected, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_disconnected <- infer_trajectory(test_disconnected, iter_n =5)
test_disconnected <- calculate_pseudotime(test_disconnected, start_state_name = c('3','4'))
test_disconnected <- plot_trajectory(test_disconnected)
plot(test_disconnected@trajectory_plot$plot)

#evaluation
eval_disconnected <- plot_index(index = factor(t(disconnected_index[, 2])), trajectory = test_disconnected@trajectory, connection_matrix = test_disconnected@connect_cluster$cluster_connection)
plot(eval_disconnected$plot)

saveRDS(test_disconnected, 'Data/cell_trajectory/dataset/test_disconnected.rds')

#binary_tree##########################  
binary_tree_data <- readRDS(file = 'Data/cell_trajectory/dataset/binary_tree_data.rds')
binary_tree_index <- readRDS(file = 'Data/cell_trajectory/dataset/binary_tree_index.rds')
binary_tree_loc_dimred <- readRDS(file = 'Data/cell_trajectory/dataset/binary_tree_loc_dimred.rds')

#1. monocle###############
binary_tree_gene <- as.data.frame(binary_tree_data[, 1])
rownames(binary_tree_gene) <- rownames(binary_tree_data)
colnames(binary_tree_gene) <- 'gene_short_name'
binary_tree_cell <- as.data.frame(binary_tree_index)
rownames(binary_tree_cell) <- colnames(binary_tree_data)


binary_tree_monocle <- infer_monocle(data = binary_tree_data,
                                      gene_annotation =  binary_tree_gene,
                                      cell_annotation = binary_tree_cell,
                                      gene_to_use = rownames(binary_tree_data))

binary_tree_monocle@phenoData@data$index = factor(t(binary_tree_index[, 2]))


plot_cell_trajectory(binary_tree_monocle)
plot_cell_trajectory(binary_tree_monocle, color_by = "Pseudotime")

saveRDS(binary_tree_monocle, 'Data/cell_trajectory/dataset/binary_tree_monocle.rds')

#2. tscan###############
binary_tree_tscan_p <- preprocess(binary_tree_data, minexpr_percent = 0.01)
binary_tree_tscan_c <- exprmclust(binary_tree_tscan_p)
binary_tree_tscan_o <- TSCANorder(binary_tree_tscan_c)

#eval
plotmclust(binary_tree_tscan_c)

saveRDS(binary_tree_tscan_c, 'Data/cell_trajectory/dataset/binary_tree_tscan_c.rds')
saveRDS(binary_tree_tscan_o, 'Data/cell_trajectory/dataset/binary_tree_tscan_o.rds')

#3. SCORPIUS###############
binary_tree_group = binary_tree_index[, 2]
binary_tree_group = as.factor(rep(1, dim(binary_tree_group)[1]))
binary_tree_sco_s <- reduce_dimensionality(t(as.matrix(binary_tree_data)), "spearman")
binary_tree_sco_t <- SCORPIUS::infer_trajectory(binary_tree_sco_s)

#group same eval
draw_trajectory_plot(binary_tree_sco_s, binary_tree_group, binary_tree_sco_t$path, contour = TRUE)
draw_trajectory_plot(binary_tree_sco_s, factor(t(binary_tree_index[, 2])), binary_tree_sco_t$path, contour = TRUE)

saveRDS(binary_tree_sco_s, 'Data/cell_trajectory/dataset/binary_tree_sco_s.rds')
saveRDS(binary_tree_sco_t, 'Data/cell_trajectory/dataset/binary_tree_sco_t.rds')

#5. tool use###########################
binary_tree_data <- as.data.frame(binary_tree_data)
test_binary_tree <- create_object(binary_tree_data, normalized = TRUE) 
test_binary_tree <- filter_data(test_binary_tree, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_binary_tree <- select_var_feature(test_binary_tree, use_normalized_data = TRUE, n = 991)
test_binary_tree <- tsneplot(test_binary_tree, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_binary_tree <- contour_plot(test_binary_tree)
test_binary_tree <- distribution_estimation(test_binary_tree, ndraw = 1000, expansion = 2.5, ... = 1, 6, 3, 4, 2, 5) 
test_binary_tree <- point_possibility(test_binary_tree, r = 5)
test_binary_tree <- connect_cluster(test_binary_tree, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_binary_tree <- infer_trajectory(test_binary_tree, iter_n =5)
test_binary_tree <- calculate_pseudotime(test_binary_tree, start_state_name = c('1'))
test_binary_tree <- plot_trajectory(test_binary_tree)
plot(test_binary_tree@trajectory_plot$plot)

#evaluation
eval_binary_tree <- plot_index(index = factor(t(binary_tree_index[, 2])), trajectory = test_binary_tree@trajectory, connection_matrix = test_binary_tree@connect_cluster$cluster_connection)
plot(eval_binary_tree$plot)

saveRDS(test_binary_tree, 'Data/cell_trajectory/dataset/test_binary_tree.rds')

#linear##########################  
linear_data <- readRDS(file = 'Data/cell_trajectory/dataset/linear_data.rds')
linear_index <- readRDS(file = 'Data/cell_trajectory/dataset/linear_index.rds')
linear_loc_dimred <- readRDS(file = 'Data/cell_trajectory/dataset/linear_loc_dimred.rds')

#1. monocle###############
linear_gene <- as.data.frame(linear_data[, 1])
rownames(linear_gene) <- rownames(linear_data)
colnames(linear_gene) <- 'gene_short_name'
linear_cell <- as.data.frame(linear_index)
rownames(linear_cell) <- colnames(linear_data)


linear_monocle <- infer_monocle(data = linear_data,
                                          gene_annotation =  linear_gene,
                                          cell_annotation = linear_cell,
                                          gene_to_use = rownames(linear_data))
linear_monocle@phenoData@data$index = factor(t(linear_index[, 2]))

plot_cell_trajectory(linear_monocle)
plot_cell_trajectory(linear_monocle, color_by = "Pseudotime")

saveRDS(linear_monocle, 'Data/cell_trajectory/dataset/test_binary_tree.rds')

#2. tscan###############
linear_tscan_p <- preprocess(linear_data, minexpr_percent = 0.01)
linear_tscan_c <- exprmclust(linear_tscan_p)
linear_tscan_o <- TSCANorder(linear_tscan_c)

#eval
plotmclust(linear_tscan_c)

saveRDS(linear_tscan_c, 'Data/cell_trajectory/dataset/linear_tscan_c.rds')
saveRDS(linear_tscan_o, 'Data/cell_trajectory/dataset/linear_tscan_o.rds')

#3. SCORPIUS###############
linear_group = linear_index[, 2]
linear_group = as.factor(rep(1, dim(linear_group)[1]))
linear_sco_s <- reduce_dimensionality(t(as.matrix(linear_data)), "spearman")
linear_sco_t <- SCORPIUS::infer_trajectory(linear_sco_s)

#group same eval
draw_trajectory_plot(linear_sco_s, linear_group, linear_sco_t$path, contour = TRUE)
draw_trajectory_plot(linear_sco_s, factor(t(linear_index[, 2])), linear_sco_t$path, contour = TRUE)

saveRDS(linear_sco_s, 'Data/cell_trajectory/dataset/linear_sco_s.rds')
saveRDS(linear_sco_t, 'Data/cell_trajectory/dataset/linear_sco_t.rds')

#5. tool use###########################
linear_data <- as.data.frame(linear_data)
test_linear <- create_object(linear_data, normalized = TRUE) 
test_linear <- filter_data(test_linear, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_linear <- select_var_feature(test_linear, use_normalized_data = TRUE, n = 1000)
test_linear <- tsneplot(test_linear, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 30)
test_linear <- contour_plot(test_linear)
test_linear <- distribution_estimation(test_linear, ndraw = 1000, expansion = 2.5, ... = c(9, 4, 3, 8), 7, c(1,6), c(2,5)) 
test_linear <- point_possibility(test_linear, r = 5)
test_linear <- connect_cluster(test_linear, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_linear <- infer_trajectory(test_linear, iter_n = 10)
test_linear <- calculate_pseudotime(test_linear, start_state_name = c('1'))
test_linear <- plot_trajectory(test_linear)
plot(test_linear@trajectory_plot$plot)

#evaluation
eval_linear <- plot_index(index = factor(t(linear_index[, 2])), trajectory = test_linear@trajectory, connection_matrix = test_linear@connect_cluster$cluster_connection)
plot(eval_linear$plot)

saveRDS(test_linear, 'Data/cell_trajectory/dataset/test_linear.rds')

#plot and measurement##########################################################
path = '/home/tialan/Data/cell_trajectory/trajectory_inference/plot/simulated_dataset_comparison/'
#p  bifurcating_cycle##########################
#p1. monocle###########

pdf(file = paste0(path, 'bifurcating_cycle_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(bifurcating_cycle_monocle)
dev.off()

pdf(file = paste0(path, 'bifurcating_cycle_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(bifurcating_cycle_monocle, color_by = "Pseudotime")
dev.off()


pdf(file = paste0(path, 'bifurcating_cycle_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(bifurcating_cycle_monocle, color_by = "index")
dev.off()

#p2. tscan###########
pdf(file = paste0(path, 'bifurcating_cycle_tscan.pdf'),width = 6,height = 6)
tscan_plot(bifurcating_cycle_tscan_c, show_cell_names = FALSE)
dev.off()

pdf(file = paste0(path, 'bifurcating_cycle_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(bifurcating_cycle_tscan_c, index = factor(t(bifurcating_cycle_index[, 2])))
dev.off()


#p3. scorpius###########
#group same eval
pdf(file = paste0(path, 'bifurcating_cycle_sco.pdf'),width = 6,height = 6)
sco_plot(bifurcating_cycle_sco_s, bifurcating_cycle_group, bifurcating_cycle_sco_t$path, contour = TRUE)
dev.off()

pdf(file = paste0(path, 'bifurcating_cycle_sco_index.pdf'),width = 6,height = 6)
sco_plot(bifurcating_cycle_sco_s, factor(t(bifurcating_cycle_index[, 2])), bifurcating_cycle_sco_t$path, contour = TRUE)
dev.off()
#p5. tool###########
pdf(file = paste0(path, 'bifurcating_cycle_tool.pdf'),width = 6,height = 6)
plot(test_bifurcating_cycle@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'bifurcating_cycle_tool_index.pdf'),width = 6,height = 6)
plot(eval_bifurcating_cycle$plot)
dev.off()

#p  cycle##########################
#p1. monocle###########
pdf(file = paste0(path, 'cycle_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(cycle_monocle)
dev.off()

pdf(file = paste0(path, 'cycle_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(cycle_monocle, color_by = "Pseudotime")
dev.off()


pdf(file = paste0(path, 'cycle_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(cycle_monocle, color_by = "index")
dev.off()
#p2. tscan###########
pdf(file = paste0(path, 'cycle_tscan.pdf'),width = 6,height = 6)
tscan_plot(cycle_tscan_c, show_cell_names = FALSE)
dev.off()

pdf(file = paste0(path, 'cycle_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(cycle_tscan_c, index = factor(t(cycle_index[, 2])))
dev.off()


#p3. scorpius###########
#group same eval
pdf(file = paste0(path, 'cycle_sco.pdf'),width = 6,height = 6)
sco_plot(cycle_sco_s, cycle_group, cycle_sco_t$path, contour = TRUE)
dev.off()

pdf(file = paste0(path, 'cycle_sco_index.pdf'),width = 6,height = 6)
sco_plot(cycle_sco_s, factor(t(cycle_index[, 2])), cycle_sco_t$path, contour = TRUE)
dev.off()

#p5. tool###########
pdf(file = paste0(path, 'cycle_tool.pdf'),width = 6,height = 6)
plot(test_cycle@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'cycle_tool_index.pdf'),width = 6,height = 6)
plot(eval_cycle$plot)
dev.off()

#p  disconnected##########################
#p1. monocle###########

pdf(file = paste0(path, 'disconnected_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(disconnected_monocle)
dev.off()

pdf(file = paste0(path, 'disconnected_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(disconnected_monocle, color_by = "Pseudotime")
dev.off()


pdf(file = paste0(path, 'disconnected_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(disconnected_monocle, color_by = "index")
dev.off()
#p2. tscan###########
pdf(file = paste0(path, 'disconnected_tscan.pdf'),width = 6,height = 6)
tscan_plot(disconnected_tscan_c, show_cell_names = FALSE)
dev.off()

pdf(file = paste0(path, 'disconnected_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(disconnected_tscan_c, index = factor(t(disconnected_index[, 2])))
dev.off()

#p3. scorpius###########
#group same eval
pdf(file = paste0(path, 'disconnected_sco.pdf'),width = 6,height = 6)
sco_plot(disconnected_sco_s, disconnected_group, disconnected_sco_t$path, contour = TRUE)
dev.off()

pdf(file = paste0(path, 'disconnected_sco_index.pdf'),width = 6,height = 6)
sco_plot(disconnected_sco_s, factor(t(disconnected_index[, 2])), disconnected_sco_t$path, contour = TRUE)
dev.off()

#p5. tool###########
pdf(file = paste0(path, 'disconnected_tool.pdf'),width = 6,height = 6)
plot(test_disconnected@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'disconnected_tool_index.pdf'),width = 6,height = 6)
plot(eval_disconnected$plot)
dev.off()
#p  binary_tree#################
#p1. monocle###########

pdf(file = paste0(path, 'binary_tree_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(binary_tree_monocle)
dev.off()

pdf(file = paste0(path, 'binary_tree_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(binary_tree_monocle, color_by = "Pseudotime")
dev.off()


pdf(file = paste0(path, 'binary_tree_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(binary_tree_monocle, color_by = "index")
dev.off()
#p2. tscan###########
pdf(file = paste0(path, 'binary_tree_tscan.pdf'),width = 6,height = 6)
tscan_plot(binary_tree_tscan_c, show_cell_names = FALSE)
dev.off()

pdf(file = paste0(path, 'binary_tree_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(binary_tree_tscan_c, index = factor(t(binary_tree_index[, 2])))
dev.off()


#p3. scorpius###########
#group same eval
pdf(file = paste0(path, 'binary_tree_sco.pdf'),width = 6,height = 6)
sco_plot(binary_tree_sco_s, binary_tree_group, binary_tree_sco_t$path, contour = TRUE)
dev.off()

pdf(file = paste0(path, 'binary_tree_sco_index.pdf'),width = 6,height = 6)
sco_plot(binary_tree_sco_s, factor(t(binary_tree_index[, 2])), binary_tree_sco_t$path, contour = TRUE)
dev.off()

#p5. tool###########
pdf(file = paste0(path, 'binary_tree_tool.pdf'),width = 6,height = 6)
plot(test_binary_tree@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'binary_tree_tool_index.pdf'),width = 6,height = 6)
plot(eval_binary_tree$plot)
dev.off()

#p  linear##########################  
#p1. monocle###########

pdf(file = paste0(path, 'linear_monocle_state.pdf'),width = 6,height = 6)
monocle_plot(linear_monocle)
dev.off()

pdf(file = paste0(path, 'linear_monocle_Pseudotime.pdf'),width = 6,height = 6)
monocle_plot(linear_monocle, color_by = "Pseudotime")
dev.off()


pdf(file = paste0(path, 'linear_monocle_index.pdf'),width = 6,height = 6)
monocle_plot(linear_monocle, color_by = "index")
dev.off()

#p2. tscan###########
pdf(file = paste0(path, 'linear_tscan.pdf'),width = 6,height = 6)
tscan_plot(linear_tscan_c, show_cell_names = FALSE)
dev.off()

pdf(file = paste0(path, 'linear_tscan_index.pdf'),width = 6,height = 6)
plot_index_tscan(linear_tscan_c, index = factor(t(linear_index[, 2])))
dev.off()


#p3. scorpius###########
#group same eval
pdf(file = paste0(path, 'linear_sco.pdf'),width = 6,height = 6)
sco_plot(linear_sco_s, linear_group, linear_sco_t$path, contour = TRUE)
dev.off()

pdf(file = paste0(path, 'linear_sco_index.pdf'),width = 6,height = 6)
sco_plot(linear_sco_s, factor(t(linear_index[, 2])), linear_sco_t$path, contour = TRUE)
dev.off()

#p5. tool###########
pdf(file = paste0(path, 'linear_tool.pdf'),width = 6,height = 6)
plot(test_linear@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'linear_tool_index.pdf'),width = 6,height = 6)
plot(eval_linear$plot)
dev.off()
