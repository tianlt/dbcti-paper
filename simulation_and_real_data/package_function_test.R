source('/home/tialan/Data/cell_trajectory/rcode/package_function.R')
#hesc(use)#########################
hescmt<-readRDS('/home/tialan/Data/cell_trajectory/dataset/hescmt.rds')
hescmt<-hescmt[, 214:460]
hescsmsinfo_cycle_index<-readRDS('/home/tialan/Data/cell_trajectory/dataset/hescsmsinfo_cycle_index.rds')
hescsmsinfo_cycle_index<-hescsmsinfo_cycle_index[214:460]
cellcycle<-c('CCNB1','CDK4','GMNN',
             'CCNB2',
             'TOP2A','MKI67',
             'CDKN1A','CDKN1B','CDKN1C',
             'CDK1','UBE2C','H4C3')

test_hescmt <- create_object(hescmt, normalized = TRUE) 
test_hescmt <- data_nmf(test_hescmt, use_normalized_data = TRUE, thread = 30)
test_hescmt <- feature_selection_knn(test_hescmt,feature = cellcycle,k = 50)
test_hescmt <- tsneplot(test_hescmt, use_normalized_data = TRUE)
test_hescmt <- contour_plot(test_hescmt)
test_hescmt <- distribution_estimation(test_hescmt, ndraw = 1000, expansion = 1.5, ... = 1,2,3)
test_hescmt <- point_possibility(test_hescmt, r = 2)
test_hescmt <- connect_cluster(test_hescmt, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)
test_hescmt <- infer_trajectory(test_hescmt, iter_n = 50)
test_hescmt <- calculate_pseudotime(test_hescmt, start_state_name = c('1'))
test_hescmt <- plot_trajectory(test_hescmt)
plot(test_hescmt@trajectory_plot$plot)

#evaluate
eval_hescmt_plot <- plot_index(index = hescsmsinfo_cycle_index, trajectory = test_hescmt@trajectory, connection_matrix = test_hescmt@connect_cluster$cluster_connection)
plot(eval_hescmt_plot$plot)

saveRDS(test_hescmt, 'Data/cell_trajectory/dataset/test_hescmt.rds')


#young (use)#############
stem_mouse_C57BL6_data_young<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_young.rds')
stem_mouse_C57BL6_index_young<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_index_young.rds')


stem_mouse_C57BL6_data_young <- as.data.frame(stem_mouse_C57BL6_data_young)
test_stem_mouse_C57BL6_young <- create_object(stem_mouse_C57BL6_data_young, normalized = TRUE) 
test_stem_mouse_C57BL6_young <- filter_data(test_stem_mouse_C57BL6_young, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_stem_mouse_C57BL6_young@specified_gene <- rownames(test_stem_mouse_C57BL6_young@normalized_data)
test_stem_mouse_C57BL6_young <- tsneplot(test_stem_mouse_C57BL6_young, use_normalized_data = TRUE, specified_gene = TRUE, pca = TRUE)
test_stem_mouse_C57BL6_young <- contour_plot(test_stem_mouse_C57BL6_young)
test_stem_mouse_C57BL6_young <- distribution_estimation(test_stem_mouse_C57BL6_young, ndraw = 1000, expansion = 1, ... = c(4,2), 1, c(3,5))
test_stem_mouse_C57BL6_young <- point_possibility(test_stem_mouse_C57BL6_young, r = 2)
test_stem_mouse_C57BL6_young <- connect_cluster(test_stem_mouse_C57BL6_young, sum_cri = 0.9, diff_cri = 0.3, vague_cri = 0.01)
test_stem_mouse_C57BL6_young <- infer_trajectory(test_stem_mouse_C57BL6_young, iter_n = 100)
test_stem_mouse_C57BL6_young <- calculate_pseudotime(test_stem_mouse_C57BL6_young, start_state_name = c('1'))
test_stem_mouse_C57BL6_young <- plot_trajectory(test_stem_mouse_C57BL6_young)
plot(test_stem_mouse_C57BL6_young@trajectory_plot$plot)
  
#evaluation
eval_stem_mouse_C57BL6_young_plot <- plot_index(index = stem_mouse_C57BL6_index_young, trajectory = test_stem_mouse_C57BL6_young@trajectory, connection_matrix = test_stem_mouse_C57BL6_young@connect_cluster$cluster_connection)
plot(eval_stem_mouse_C57BL6_young_plot$plot)

saveRDS(test_stem_mouse_C57BL6_young, 'Data/cell_trajectory/dataset/test_stem_mouse_C57BL6_young.rds')

#old (use)###########
stem_mouse_C57BL6_data_old<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_data_old.rds')
stem_mouse_C57BL6_index_old<-readRDS('/home/tialan/Data/cell_trajectory/dataset/stem_mouse_C57BL6_index_old.rds')

  
stem_mouse_C57BL6_data_old <- as.data.frame(stem_mouse_C57BL6_data_old)
test_stem_mouse_C57BL6_old <- create_object(stem_mouse_C57BL6_data_old, normalized = TRUE) 
test_stem_mouse_C57BL6_old <- filter_data(test_stem_mouse_C57BL6_old, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_stem_mouse_C57BL6_old@specified_gene <- rownames(test_stem_mouse_C57BL6_old@normalized_data)
test_stem_mouse_C57BL6_old <- tsneplot(test_stem_mouse_C57BL6_old, use_normalized_data = TRUE, specified_gene = TRUE, pca = TRUE)
test_stem_mouse_C57BL6_old <- contour_plot(test_stem_mouse_C57BL6_old)
test_stem_mouse_C57BL6_old <- distribution_estimation(test_stem_mouse_C57BL6_old, ndraw = 1000, expansion = 1.5, ... = 1, 2, 3)
test_stem_mouse_C57BL6_old <- point_possibility(test_stem_mouse_C57BL6_old, r = 2)
test_stem_mouse_C57BL6_old <- connect_cluster(test_stem_mouse_C57BL6_old, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)
test_stem_mouse_C57BL6_old <- infer_trajectory(test_stem_mouse_C57BL6_old, iter_n = 1000)
test_stem_mouse_C57BL6_old <- calculate_pseudotime(test_stem_mouse_C57BL6_old, start_state_name = c('1'))
test_stem_mouse_C57BL6_old <- plot_trajectory(test_stem_mouse_C57BL6_old, width = 0.1, height = 0.1)
plot(test_stem_mouse_C57BL6_old@trajectory_plot$plot)


#evaluation
eval_stem_mouse_C57BL6_old_plot <- plot_index(index = stem_mouse_C57BL6_index_old, trajectory = test_stem_mouse_C57BL6_old@trajectory, connection_matrix = test_stem_mouse_C57BL6_old@connect_cluster$cluster_connection, jitter = 0.1)
plot(eval_stem_mouse_C57BL6_old_plot$plot)

saveRDS(test_stem_mouse_C57BL6_old, 'Data/cell_trajectory/dataset/test_stem_mouse_C57BL6_old.rds')


#batch1(use) 2000gene###############
camp1_data_batch1 <- readRDS('/home/tialan/Data/cell_trajectory/dataset/camp1_data_batch1.rds')
camp1_index_batch1 <- readRDS('/home/tialan/Data/cell_trajectory/dataset/camp1_index_batch1.rds')

camp1_data_batch1 <- as.data.frame(camp1_data_batch1)
test_camp1_batch1 <- create_object(camp1_data_batch1, normalized = TRUE) 
test_camp1_batch1 <- filter_data(test_camp1_batch1, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_camp1_batch1 <- select_var_feature(test_camp1_batch1, use_normalized_data = TRUE, n = 2000)
test_camp1_batch1 <- tsneplot(test_camp1_batch1, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE)
test_camp1_batch1 <- contour_plot(test_camp1_batch1)
test_camp1_batch1 <- distribution_estimation(test_camp1_batch1, ndraw = 1000, expansion = 2, ... = 5, 4, c(1, 2), 3)
test_camp1_batch1 <- point_possibility(test_camp1_batch1, r = 3)
test_camp1_batch1 <- connect_cluster(test_camp1_batch1, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_camp1_batch1 <- infer_trajectory(test_camp1_batch1, iter_n =500)
test_camp1_batch1 <- calculate_pseudotime(test_camp1_batch1, start_state_name = c('1', '2'))
test_camp1_batch1 <- plot_trajectory(test_camp1_batch1)
plot(test_camp1_batch1@trajectory_plot$plot)

#evaluation
eval_camp1_batch1_plot <- plot_index(index = factor(camp1_index_batch1[2,]), trajectory = test_camp1_batch1@trajectory, connection_matrix = test_camp1_batch1@connect_cluster$cluster_connection)
plot(eval_camp1_batch1_plot$plot)

saveRDS(test_camp1_batch1@selected_feature$selected_gene, file = 'Data/cell_trajectory/dataset/camp1_batch1_selected_gene.rds')
saveRDS(test_camp1_batch1, 'Data/cell_trajectory/dataset/test_camp1_batch1.rds')




#nestorowa_data_batch1_type_filter no stem(use1) 2000gene####################
nestorowa_data_batch1_type_filter_no_stem <- readRDS('/home/tialan/Data/cell_trajectory/dataset/nestorowa_data_batch1_type_filter_no_stem.rds')
nestorowa_index_batch1_type_filter_no_stem <- readRDS('/home/tialan/Data/cell_trajectory/dataset/nestorowa_index_batch1_type_filter_no_stem.rds')

nestorowa_data_batch1_type_filter_no_stem <- as.data.frame(nestorowa_data_batch1_type_filter_no_stem)
test_nestorowa_batch1_type_filter_no_stem <- create_object(nestorowa_data_batch1_type_filter_no_stem, normalized = TRUE) 
test_nestorowa_batch1_type_filter_no_stem <- filter_data(test_nestorowa_batch1_type_filter_no_stem, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_nestorowa_batch1_type_filter_no_stem <- select_var_feature(test_nestorowa_batch1_type_filter_no_stem, use_normalized_data = TRUE, n = 2000)
test_nestorowa_batch1_type_filter_no_stem <- tsneplot(test_nestorowa_batch1_type_filter_no_stem, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE)
test_nestorowa_batch1_type_filter_no_stem <- contour_plot(test_nestorowa_batch1_type_filter_no_stem)
test_nestorowa_batch1_type_filter_no_stem <- distribution_estimation(test_nestorowa_batch1_type_filter_no_stem, ndraw = 1000, expansion = 1.5, ... = 1, 2, 3, c(4,5)) 
test_nestorowa_batch1_type_filter_no_stem <- point_possibility(test_nestorowa_batch1_type_filter_no_stem, r = 3)
test_nestorowa_batch1_type_filter_no_stem <- connect_cluster(test_nestorowa_batch1_type_filter_no_stem, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_nestorowa_batch1_type_filter_no_stem <- infer_trajectory(test_nestorowa_batch1_type_filter_no_stem, iter_n =200)
test_nestorowa_batch1_type_filter_no_stem <- calculate_pseudotime(test_nestorowa_batch1_type_filter_no_stem, start_state_name = c('1'))
test_nestorowa_batch1_type_filter_no_stem <- plot_trajectory(test_nestorowa_batch1_type_filter_no_stem)
plot(test_nestorowa_batch1_type_filter_no_stem@trajectory_plot$plot)

#evaluation
eval_nestorowa_batch1_type_filter_no_stem <- plot_index(index = factor(nestorowa_index_batch1_type_filter_no_stem[2,]), trajectory = test_nestorowa_batch1_type_filter_no_stem@trajectory, connection_matrix = test_nestorowa_batch1_type_filter_no_stem@connect_cluster$cluster_connection)
plot(eval_nestorowa_batch1_type_filter_no_stem$plot)

saveRDS(test_nestorowa_batch1_type_filter_no_stem@selected_feature$selected_gene, file = 'Data/cell_trajectory/dataset/nestorowa_batch1_type_filter_no_stem_selected_gene.rds')
saveRDS(test_nestorowa_batch1_type_filter_no_stem, 'Data/cell_trajectory/dataset/test_nestorowa_batch1_type_filter_no_stem.rds')

#yan human embryo(gse36552) (hemberg lab github) (normalized) (use) 2000gene#########################
yan_data <- readRDS('/home/tialan/Data/cell_trajectory/dataset/yan_data.rds')
yan_index <- readRDS('/home/tialan/Data/cell_trajectory/dataset/yan_index.rds')

yan_data <- as.data.frame(yan_data)
test_yan <- create_object(yan_data, normalized = TRUE) 
test_yan <- filter_data(test_yan, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
test_yan <- select_var_feature(test_yan, use_normalized_data = TRUE, n = 2000)
test_yan <- tsneplot(test_yan, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 15)
test_yan <- contour_plot(test_yan)
test_yan <- distribution_estimation(test_yan, ndraw = 1000, expansion = 2.5, ... = c(1,2), 3, c(4,5), c(6,7,8)) 
test_yan <- point_possibility(test_yan, r = 5)
test_yan <- connect_cluster(test_yan, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.01)
test_yan <- infer_trajectory(test_yan, iter_n =200)
test_yan <- calculate_pseudotime(test_yan, start_state_name = c('1'))
test_yan <- plot_trajectory(test_yan)
plot(test_yan@trajectory_plot$plot)

#evaluation
eval_yan <- plot_index(index = factor(yan_index[1,]), trajectory = test_yan@trajectory, connection_matrix = test_yan@connect_cluster$cluster_connection)
plot(eval_yan$plot)

saveRDS(test_yan@selected_feature$selected_gene, file = 'Data/cell_trajectory/dataset/yan_selected_gene.rds')
saveRDS(test_yan, 'Data/cell_trajectory/dataset/test_yan.rds')


#plot and measurement#####################
path = '/home/tialan/Data/cell_trajectory/trajectory_inference/plot/package_comparison/'
#p  hesc(use)#########################
pdf(file = paste0(path, 'hesc_tool.pdf'),width = 6,height = 6)
plot(test_hescmt@trajectory_plot$plot)
dev.off()

pdf(file = paste0(path, 'hesc_tool_index.pdf'),width = 6,height = 6)
plot(eval_hescmt_plot$plot)
dev.off()
#p  young (use)#############
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_tool.pdf'),width = 6,height = 6)
plot(test_stem_mouse_C57BL6_young@trajectory_plot$plot)
dev.off()
pdf(file = paste0(path, 'stem_mouse_C57BL6_young_tool_index.pdf'),width = 6,height = 6)
plot(eval_stem_mouse_C57BL6_young_plot$plot)
dev.off()
#p  old (use)###########
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_tool.pdf'),width = 6,height = 6)
plot(test_stem_mouse_C57BL6_old@trajectory_plot$plot)
dev.off()
pdf(file = paste0(path, 'stem_mouse_C57BL6_old_tool_index.pdf'),width = 6,height = 6)
plot(eval_stem_mouse_C57BL6_old_plot$plot)
dev.off()
#p  batch1(use) 2000gene###############
pdf(file = paste0(path, 'camp1_batch1_tool.pdf'),width = 6,height = 6)
plot(test_camp1_batch1@trajectory_plot$plot)
dev.off()
pdf(file = paste0(path, 'camp1_batch1_tool_index.pdf'),width = 6,height = 6)
plot(eval_camp1_batch1_plot$plot)
dev.off()

#p  nestorowa_data_batch1_type_filter no stem(use1) 2000gene####################
pdf(file = paste0(path, 'nestorowa_batch1_type_filter_no_stem_tool.pdf'),width = 6,height = 6)
plot(test_nestorowa_batch1_type_filter_no_stem@trajectory_plot$plot)
dev.off()
pdf(file = paste0(path, 'nestorowa_batch1_type_filter_no_stem_tool_index.pdf'),width = 6,height = 6)
plot(eval_nestorowa_batch1_type_filter_no_stem$plot)
dev.off()


#p  yan human embryo(gse36552) (hemberg lab github) (normalized) (use) 2000gene#########################
pdf(file = paste0(path, 'yan_tool.pdf'),width = 6,height = 6)
plot(test_yan@trajectory_plot$plot)
dev.off()
pdf(file = paste0(path, 'yan_tool_index.pdf'),width = 6,height = 6)
plot(eval_yan$plot)
dev.off()