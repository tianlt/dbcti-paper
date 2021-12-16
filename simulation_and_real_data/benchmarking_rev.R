library(dynplot)
library(dbcti)
library(monocle)
library(tidyverse)
library(SCORPIUS)
library(TSCAN)
library(dyno)
library(dyngen)
library(purrr)
#backbone_bifurcating_cycle###############
bifurcating_cycle <- function(x){
  set.seed(x)
  backbone_bifurcating_cycle <- backbone_bifurcating_cycle()
  config_bifurcating_cycle <- 
    initialise_model(
      backbone = backbone_bifurcating_cycle,
      num_tfs = nrow(backbone_bifurcating_cycle$module_info),
      num_targets = 500,
      num_hks = 500,
      verbose = FALSE, 
      num_cores = 8
    )
  
  out_bifurcating_cycle <- generate_dataset(
    config_bifurcating_cycle,
    format = "dyno",
    make_plots = TRUE
  )
  return(out_bifurcating_cycle)
}

#backbone_cycle###############
cycle <- function(x){
  set.seed(x)
  backbone_cycle <- backbone_cycle()
  config_cycle <- 
    initialise_model(
      backbone = backbone_cycle,
      num_tfs = nrow(backbone_cycle$module_info),
      num_targets = 500,
      num_hks = 500,
      verbose = FALSE
    )
  
  out_cycle <- generate_dataset(
    config_cycle,
    format = "dyno",
    make_plots = TRUE
  )
  return(out_cycle)
  
  
}

#backbone_disconnected###############
disconnected <- function(x){
  set.seed(x)
  backbone_disconnected <- backbone_disconnected()
  config_disconnected <- 
    initialise_model(
      backbone = backbone_disconnected,
      num_tfs = nrow(backbone_disconnected$module_info),
      num_targets = 500,
      num_hks = 500,
      verbose = FALSE
    )
  
  out_disconnected <- generate_dataset(
    config_disconnected,
    format = "dyno",
    make_plots = TRUE
  )
  return(out_disconnected)
}

#backbone_binary_tree###############
binary_tree <- function(x){
  set.seed(x)
  backbone_binary_tree <- backbone_binary_tree(
    num_modifications = 2
  )
  
  config_binary_tree <- 
    initialise_model(
      backbone = backbone_binary_tree,
      num_tfs = nrow(backbone_binary_tree$module_info),
      num_targets = 500,
      num_hks = 500,
      verbose = FALSE,
      num_cores = 8
    )
  
  out_binary_tree <- generate_dataset(
    config_binary_tree,
    format = "dyno",
    make_plots = TRUE
  )
  return(out_binary_tree)
}

#backbone_linear###############
linear <- function(x){
  set.seed(x)
  backbone_linear <- backbone_linear()
  
  config_linear <- 
    initialise_model(
      backbone = backbone_linear,
      num_tfs = nrow(backbone_linear$module_info),
      num_targets = 500,
      num_hks = 500,
      verbose = FALSE,
      num_cores = 8
    )
  
  out_linear <- generate_dataset(config_linear,
                                 format = 'dyno', 
                                 make_plots = TRUE)
  return(out_linear)
}
#data import#############################
bc_0 <- readRDS('out_bifurcating_cycle.rds')
bc_2 <- readRDS('bc_2.rds')
bc_3 <- readRDS('bc_3.rds')

bt_0 <- readRDS('out_binary_tree.rds')
bt_2 <- readRDS('bt_2.rds')
bt_3 <- readRDS('bt_3.rds')

l_0 <- readRDS('out_linear.rds')
l_2 <- readRDS('l_2.rds')
l_3<- readRDS('l_3.rds')

c_0 <- readRDS('out_cycle.rds')
c_2 <- readRDS('c_2.rds')
c_3 <- readRDS('c_3.rds')


d_0 <- readRDS('out_disconnected.rds')
d_2 <- readRDS('d_2.rds')
d_3 <- readRDS('d_3.rds')

leng_data <- readRDS('hescmt_for_new.rds')
leng_index <- readRDS('hescsmsinfo_cycle_index_for_new.rds')

kowalczyk_y_data <- readRDS('stem_mouse_C57BL6_data_young_for_new.rds')
kowalczyk_y_index <- readRDS('stem_mouse_C57BL6_index_young_for_new.rds')

kowalczyk_o_data <- readRDS('stem_mouse_C57BL6_data_old_for_new.rds')
kowalczyk_o_index <- readRDS('stem_mouse_C57BL6_index_old_for_new.rds')

camp_data <- readRDS('camp1_data_batch1_for_new.rds')
camp_index <- readRDS('camp1_index_batch1_for_new.rds')

nestorowa_data <- readRDS('nestorowa_data_batch1_type_filter_no_stem_for_new.rds')
nestorowa_index <- readRDS('nestorowa_index_batch1_type_filter_no_stem_for_new.rds')

yan_data <- readRDS('yan_data_for_new.rds')
yan_index <- readRDS('yan_index_for_new.rds')


ref_name <- c('bc_0','bc_2','bc_3','bt_0','bt_2','bt_3','l_0','l_2','l_3','c_0','c_2','c_3','d_0','d_2','d_3')
for (i in ref_name) {
  pdf(file = paste0('diagram/', i, '.pdf'), width = 6,height = 6)
  print(plot_dimred(eval(parse(text = paste0(i,'$dataset')))))
  dev.off()
}







dataset_bc_2_df = as.data.frame(t(as.matrix(bc_2$dataset[["expression"]]))) #dataset matrix
bc_2_data = dataset_bc_2_df[apply(dataset_bc_2_df, 1, sum) != 0, ]
dataset_bc_3_df = as.data.frame(t(as.matrix(bc_3$dataset[["expression"]]))) #dataset matrix
bc_3_data = dataset_bc_3_df[apply(dataset_bc_3_df, 1, sum) != 0, ]


plot_dimred(c_1$dataset)
plot_dimred(c_2$dataset)




dataset_c_2_df = as.data.frame(t(as.matrix(c_2$dataset[["expression"]]))) #dataset matrix
c_2_data = dataset_c_2_df[apply(dataset_c_2_df, 1, sum) != 0, ]
dataset_c_3_df = as.data.frame(t(as.matrix(c_3$dataset[["expression"]]))) #dataset matrix
c_3_data = dataset_c_3_df[apply(dataset_c_3_df, 1, sum) != 0, ]



dataset_l_2_df = as.data.frame(t(as.matrix(l_2$dataset[["expression"]]))) #dataset matrix
l_2_data = dataset_l_2_df[apply(dataset_l_2_df, 1, sum) != 0, ]
dataset_l_3_df = as.data.frame(t(as.matrix(l_3$dataset[["expression"]]))) #dataset matrix
l_3_data = dataset_l_3_df[apply(dataset_l_3_df, 1, sum) != 0, ]


dataset_d_2_df = as.data.frame(t(as.matrix(d_2$dataset[["expression"]]))) #dataset matrix
d_2_data = dataset_d_2_df[apply(dataset_d_2_df, 1, sum) != 0, ]
dataset_d_3_df = as.data.frame(t(as.matrix(d_3$dataset[["expression"]]))) #dataset matrix
d_3_data = dataset_d_3_df[apply(dataset_d_3_df, 1, sum) != 0, ]



dataset_bt_2_df = as.data.frame(t(as.matrix(bt_2$dataset[["expression"]]))) #dataset matrix
bt_2_data = dataset_bt_2_df[apply(dataset_bt_2_df, 1, sum) != 0, ]
dataset_bt_3_df = as.data.frame(t(as.matrix(bt_3$dataset[["expression"]]))) #dataset matrix
bt_3_data = dataset_bt_3_df[apply(dataset_bt_3_df, 1, sum) != 0, ]


dataset_bt_0_df = as.data.frame(t(as.matrix(bt_0$dataset[["expression"]]))) #dataset matrix
bt_0_data = dataset_bt_0_df[apply(dataset_bt_0_df, 1, sum) != 0, ]
dataset_bc_0_df = as.data.frame(t(as.matrix(bc_0$dataset[["expression"]]))) #dataset matrix
bc_0_data = dataset_bc_0_df[apply(dataset_bc_0_df, 1, sum) != 0, ]
dataset_l_0_df = as.data.frame(t(as.matrix(l_0$dataset[["expression"]]))) #dataset matrix
l_0_data = dataset_l_0_df[apply(dataset_l_0_df, 1, sum) != 0, ]
dataset_c_0_df = as.data.frame(t(as.matrix(c_0$dataset[["expression"]]))) #dataset matrix
c_0_data = dataset_c_0_df[apply(dataset_c_0_df, 1, sum) != 0, ]
dataset_d_0_df = as.data.frame(t(as.matrix(d_0$dataset[["expression"]]))) #dataset matrix
d_0_data = dataset_d_0_df[apply(dataset_d_0_df, 1, sum) != 0, ]


bc_2_data <- as.data.frame(bc_2_data)
bc_2_dbcti <- create_object(bc_2_data, normalized = TRUE) 
bc_2_dbcti <- filter_data(bc_2_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
bc_2_dbcti <- select_var_feature(bc_2_dbcti, use_normalized_data = TRUE, n = 1000)
bc_2_dbcti <- tsneplot(bc_2_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
bc_2_dbcti <- contour_plot(bc_2_dbcti)


bc_2_dbcti <- distribution_estimation(bc_2_dbcti, ndraw = 1000, expansion = 2.5, ... = c(2,6),c(1,8),c(5,7),c(3,4)) 
bc_2_dbcti <- point_possibility(bc_2_dbcti, r = 5)
bc_2_dbcti <- connect_cluster(bc_2_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
bc_2_dbcti <- infer_trajectory(bc_2_dbcti, iter_n =3)
bc_2_dbcti <- calculate_pseudotime(bc_2_dbcti, start_state_name = c('1'))
bc_2_dbcti <- plot_trajectory(bc_2_dbcti)
plot(bc_2_dbcti@trajectory_plot$plot)



bc_3_data <- as.data.frame(bc_3_data)
bc_3_dbcti <- create_object(bc_3_data, normalized = TRUE) 
bc_3_dbcti <- filter_data(bc_3_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
bc_3_dbcti <- select_var_feature(bc_3_dbcti, use_normalized_data = TRUE, n = 1000)
bc_3_dbcti <- tsneplot(bc_3_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
bc_3_dbcti <- contour_plot(bc_3_dbcti)


bc_3_dbcti <- distribution_estimation(bc_3_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,2,3,4,5,6,7) 
bc_3_dbcti <- point_possibility(bc_3_dbcti, r = 5)
bc_3_dbcti <- connect_cluster(bc_3_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
bc_3_dbcti <- infer_trajectory(bc_3_dbcti, iter_n =5)
bc_3_dbcti <- calculate_pseudotime(bc_3_dbcti, start_state_name = c('1'))
bc_3_dbcti <- plot_trajectory(bc_3_dbcti)
plot(bc_3_dbcti@trajectory_plot$plot)


l_2_data <- as.data.frame(l_2_data)
l_2_dbcti <- create_object(l_2_data, normalized = TRUE) 
l_2_dbcti <- filter_data(l_2_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
l_2_dbcti <- select_var_feature(l_2_dbcti, use_normalized_data = TRUE, n = 500)
l_2_dbcti <- tsneplot(l_2_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
l_2_dbcti <- contour_plot(l_2_dbcti)


l_2_dbcti <- distribution_estimation(l_2_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,2,3) 
l_2_dbcti <- point_possibility(l_2_dbcti, r = 5)
l_2_dbcti <- connect_cluster(l_2_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
l_2_dbcti <- infer_trajectory(l_2_dbcti, iter_n =15)
l_2_dbcti <- calculate_pseudotime(l_2_dbcti, start_state_name = c('1'))
l_2_dbcti <- plot_trajectory(l_2_dbcti)
plot(l_2_dbcti@trajectory_plot$plot)


l_3_data <- as.data.frame(l_3_data)
l_3_dbcti <- create_object(l_3_data, normalized = TRUE) 
l_3_dbcti <- filter_data(l_3_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
l_3_dbcti <- select_var_feature(l_3_dbcti, use_normalized_data = TRUE, n = 500)
l_3_dbcti <- tsneplot(l_3_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
l_3_dbcti <- contour_plot(l_3_dbcti)


l_3_dbcti <- distribution_estimation(l_3_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,2,3) 
l_3_dbcti <- point_possibility(l_3_dbcti, r = 5)
l_3_dbcti <- connect_cluster(l_3_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
l_3_dbcti <- infer_trajectory(l_3_dbcti, iter_n =25)
l_3_dbcti <- calculate_pseudotime(l_3_dbcti, start_state_name = c('2'))
l_3_dbcti <- plot_trajectory(l_3_dbcti)
plot(l_3_dbcti@trajectory_plot$plot)



#c#####################################
#c_2###########################################
c_2_data <- as.data.frame(c_2_data)
c_2_dbcti <- create_object(c_2_data, normalized = TRUE) 
c_2_dbcti <- filter_data(c_2_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
c_2_dbcti <- select_var_feature(c_2_dbcti, use_normalized_data = TRUE, n = 500)
c_2_dbcti <- tsneplot(c_2_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
c_2_dbcti <- contour_plot(c_2_dbcti)


c_2_dbcti <- distribution_estimation(c_2_dbcti, ndraw = 1000, expansion = 2.5, ... = c(1,3),2,4)
c_2_dbcti <- point_possibility(c_2_dbcti, r = 5)
c_2_dbcti <- connect_cluster(c_2_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
c_2_dbcti <- infer_trajectory(c_2_dbcti, iter_n =15)
c_2_dbcti <- calculate_pseudotime(c_2_dbcti, start_state_name = c('1'))
c_2_dbcti <- plot_trajectory(c_2_dbcti)
plot(c_2_dbcti@trajectory_plot$plot)

#c_3#################################################
c_3_data <- as.data.frame(c_3_data)
c_3_dbcti <- create_object(c_3_data, normalized = TRUE) 
c_3_dbcti <- filter_data(c_3_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
c_3_dbcti <- select_var_feature(c_3_dbcti, use_normalized_data = TRUE, n = 500)
c_3_dbcti <- tsneplot(c_3_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
c_3_dbcti <- contour_plot(c_3_dbcti)


c_3_dbcti <- distribution_estimation(c_3_dbcti, ndraw = 1000, expansion = 2.5, ... = c(1,7),2,3,c(4,5,6)) 
c_3_dbcti <- point_possibility(c_3_dbcti, r = 5)
c_3_dbcti <- connect_cluster(c_3_dbcti, sum_cri = 0.7, diff_cri = 0.4, vague_cri = 0.01)
c_3_dbcti <- infer_trajectory(c_3_dbcti, iter_n =15)
c_3_dbcti <- calculate_pseudotime(c_3_dbcti, start_state_name = c('1'))
c_3_dbcti <- plot_trajectory(c_3_dbcti)
plot(c_3_dbcti@trajectory_plot$plot)

#d#####################################
#d_2#################################################
d_2_data <- as.data.frame(d_2_data)
d_2_dbcti <- create_object(d_2_data, normalized = TRUE) 
d_2_dbcti <- filter_data(d_2_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
d_2_dbcti <- select_var_feature(d_2_dbcti, use_normalized_data = TRUE, n = 500)
d_2_dbcti <- tsneplot(d_2_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
d_2_dbcti <- contour_plot(d_2_dbcti)


d_2_dbcti <- distribution_estimation(d_2_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,2) 
d_2_dbcti <- point_possibility(d_2_dbcti, r = 5)
d_2_dbcti <- connect_cluster(d_2_dbcti, sum_cri = 0.7, diff_cri = 0.3, vague_cri = 0.01)
d_2_dbcti <- infer_trajectory(d_2_dbcti, iter_n =15)
d_2_dbcti <- calculate_pseudotime(d_2_dbcti, start_state_name = c('1','2'))
d_2_dbcti <- plot_trajectory(d_2_dbcti)
plot(d_2_dbcti@trajectory_plot$plot)





#d_3#################################################
d_3_data <- as.data.frame(d_3_data)
d_3_dbcti <- create_object(d_3_data, normalized = TRUE) 
d_3_dbcti <- filter_data(d_3_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
d_3_dbcti <- select_var_feature(d_3_dbcti, use_normalized_data = TRUE, n = 500)
d_3_dbcti <- tsneplot(d_3_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
d_3_dbcti <- contour_plot(d_3_dbcti)


d_3_dbcti <- distribution_estimation(d_3_dbcti, ndraw = 1000, expansion = 2.5, ... = c(2,3),7,c(1,8),c(4,5,6)) 
d_3_dbcti <- point_possibility(d_3_dbcti, r = 5)
d_3_dbcti <- connect_cluster(d_3_dbcti, sum_cri = 0.7, diff_cri = 0.2, vague_cri = 0.01)
d_3_dbcti <- infer_trajectory(d_3_dbcti, iter_n =15)
d_3_dbcti <- calculate_pseudotime(d_3_dbcti, start_state_name = c('1','3'))
d_3_dbcti <- plot_trajectory(d_3_dbcti)
plot(d_3_dbcti@trajectory_plot$plot)



#bt#################################################
#bt_2#################################################
bt_2_data <- as.data.frame(bt_2_data)
bt_2_dbcti <- create_object(bt_2_data, normalized = TRUE) 
bt_2_dbcti <- filter_data(bt_2_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
bt_2_dbcti <- select_var_feature(bt_2_dbcti, use_normalized_data = TRUE, n = 1000)
bt_2_dbcti <- tsneplot(bt_2_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
bt_2_dbcti <- contour_plot(bt_2_dbcti)


bt_2_dbcti <- distribution_estimation(bt_2_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,c(3,4),5,c(2,6)) 
bt_2_dbcti <- point_possibility(bt_2_dbcti, r = 5)
bt_2_dbcti <- connect_cluster(bt_2_dbcti, sum_cri = 0.7, diff_cri = 0.2, vague_cri = 0.01)
bt_2_dbcti <- infer_trajectory(bt_2_dbcti, iter_n =15)
bt_2_dbcti <- calculate_pseudotime(bt_2_dbcti, start_state_name = c('1'))
bt_2_dbcti <- plot_trajectory(bt_2_dbcti)
plot(bt_2_dbcti@trajectory_plot$plot)


#bt_3#################################################
bt_3_data <- as.data.frame(bt_3_data)
bt_3_dbcti <- create_object(bt_3_data, normalized = TRUE) 
bt_3_dbcti <- filter_data(bt_3_dbcti, gene_cri = 1, cell_cri = 1, use_normalized_data = TRUE) 
bt_3_dbcti <- select_var_feature(bt_3_dbcti, use_normalized_data = TRUE, n = 1000)
bt_3_dbcti <- tsneplot(bt_3_dbcti, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
bt_3_dbcti <- contour_plot(bt_3_dbcti)


bt_3_dbcti <- distribution_estimation(bt_3_dbcti, ndraw = 1000, expansion = 2.5, ... = 1,c(3,5),2,4,6) 
bt_3_dbcti <- point_possibility(bt_3_dbcti, r = 5)
bt_3_dbcti <- connect_cluster(bt_3_dbcti, sum_cri = 0.7, diff_cri = 0.3, vague_cri = 0.01)
bt_3_dbcti <- infer_trajectory(bt_3_dbcti, iter_n =15)
bt_3_dbcti <- calculate_pseudotime(bt_3_dbcti, start_state_name = c('1'))
bt_3_dbcti <- plot_trajectory(bt_3_dbcti)
plot(bt_3_dbcti@trajectory_plot$plot)


#save_plot
simulated_dbcti_plot_list <- vector(mode = 'list', length = 10)
names(simulated_dbcti_plot_list) <- c('bc_2','bc_3','bt_2','bt_3','l_2','l_3','c_2','c_3','d_2','d_3')


for (i in names(simulated_dbcti_plot_list)) {
  simulated_dbcti_plot_list[[i]] <- eval(parse(text = paste0(i,'_dbcti@trajectory_plot$plot')))
}

plot_data_dbcti <- function(data){
  for (i in names(data)){
    pdf(file = paste0('diagram/', i, '_', 'model_dbcti', '.pdf'),width = 6,height = 6)
    print(plot(data[[i]]))
    dev.off()
  }
}
plot_data_dbcti(simulated_dbcti_plot_list)

########################################################################################################################
c_2_monocle <- monocle_infer(c_2_data)
c_2_monocle <- monocle_infer(c_2_data)c_2_monocle <- monocle_infer(c_2_data)c_2_monocle <- monocle_infer(c_2_data)
#monocle
function(data){
  model_monocle <- monocle_infer(data)
  
}











#bc2
bc_2_gene <- as.data.frame(bc_2_data[, 1])
rownames(bc_2_gene) <- rownames(bc_2_data)
colnames(bc_2_gene) <- 'gene_short_name'
bc_2_cell <- as.data.frame(t(bc_2_data[1,]))
rownames(bc_2_cell) <- colnames(bc_2_data)


bc_2_monocle <- infer_monocle(data = bc_2_data,
                              gene_annotation =  bc_2_gene,
                              cell_annotation = bc_2_cell,
                              gene_to_use = rownames(bc_2_data))


monocle_plot(bc_2_monocle, color_by = "Pseudotime")

#bc9
bc_3_gene <- as.data.frame(bc_3_data[, 1])
rownames(bc_3_gene) <- rownames(bc_3_data)
colnames(bc_3_gene) <- 'gene_short_name'
bc_3_cell <- as.data.frame(t(bc_3_data[1,]))
rownames(bc_3_cell) <- colnames(bc_3_data)


bc_3_monocle <- infer_monocle(data = bc_3_data,
                              gene_annotation =  bc_3_gene,
                              cell_annotation = bc_3_cell,
                              gene_to_use = rownames(bc_3_data))



monocle_plot(bc_3_monocle, color_by = "Pseudotime")

#c
c_2_monocle <- monocle_infer(c_2_data)
monocle_plot(c_2_monocle, color_by = "Pseudotime")
c_3_monocle <- monocle_infer(c_3_data)
monocle_plot(c_3_monocle, color_by = "Pseudotime")

#d
d_2_monocle <- monocle_infer(d_2_data)
monocle_plot(d_2_monocle, color_by = "Pseudotime")
d_3_monocle <- monocle_infer(d_39_data)
monocle_plot(d_3_monocle, color_by = "Pseudotime")

#d
d_2_monocle <- monocle_infer(d_2_data)
monocle_plot(d_2_monocle, color_by = "Pseudotime")
d_3_monocle <- monocle_infer(d_3_data)
monocle_plot(d_3_monocle, color_by = "Pseudotime")

#l
l_2_monocle <- monocle_infer(l_2_data)
monocle_plot(l_2_monocle, color_by = "Pseudotime")
l_3_monocle <- monocle_infer(l_3_data)
monocle_plot(l_3_monocle, color_by = "Pseudotime")

#bt
bt_2_monocle <- monocle_infer(bt_2_data)
monocle_plot(bt_2_monocle, color_by = "Pseudotime")
bt_3_monocle <- monocle_infer(bt_3_data)
monocle_plot(bt_3_monocle, color_by = "Pseudotime")

#0
bt_0_monocle <- monocle_infer(bt_0_data)
bc_0_monocle <- monocle_infer(bc_0_data)
l_0_monocle <- monocle_infer(l_0_data)
c_0_monocle <- monocle_infer(c_0_data)
d_0_monocle <- monocle_infer(d_0_data)

simulated_monocle_plot_list <- vector(mode = 'list', length = 15)
names(simulated_monocle_plot_list) <- c('bc_0','bc_2','bc_3','bt_0','bt_2','bt_3','l_0','l_2','l_3','c_0','c_2','c_3','d_0','d_2','d_3')


for (i in names(simulated_monocle_plot_list)) {
  simulated_monocle_plot_list[[i]] <- eval(parse(text = paste0(i,'_monocle')))
}
saveRDS(simulated_monocle_plot_list, 'simulated_monocle_plot_list.rds')
plot_data_monocle <- function(data){
  for (i in names(data)){
    pdf(file = paste0('diagram/', i, '_', 'model_monocle', '.pdf'),width = 6,height = 6)
    print(monocle_plot(data[[i]], color_by = "Pseudotime"))
    dev.off()
  }
}

plot_data_monocle(simulated_monocle_plot_list)

#########save rds
saveRDS(bc_2_monocle, 'bc_2_monocle.rds')
saveRDS(bc_3_monocle, 'bc_3_monocle.rds')

saveRDS(c_2_monocle, 'c_2_monocle.rds')
saveRDS(c_3_monocle, 'c_3_monocle.rds')

saveRDS(bt_2_monocle, 'bt_2_monocle.rds')
saveRDS(bt_3_monocle, 'bt_3_monocle.rds')

saveRDS(d_2_monocle, 'd_2_monocle.rds')
saveRDS(d_3_monocle, 'd_3_monocle.rds')

saveRDS(l_2_monocle, 'l_2_monocle.rds')
saveRDS(l_3_monocle, 'l_3_monocle.rds')

saveRDS(bc_2_dbcti, 'bc_2_dbcti.rds')
saveRDS(bc_3_dbcti, 'bc_3_dbcti.rds')

saveRDS(c_2_dbcti, 'c_2_dbcti.rds')
saveRDS(c_3_dbcti, 'c_3_dbcti.rds')

saveRDS(bt_2_dbcti, 'bt_2_dbcti.rds')
saveRDS(bt_3_dbcti, 'bt_3_dbcti.rds')

saveRDS(d_2_dbcti, 'd_2_dbcti.rds')
saveRDS(d_3_dbcti, 'd_3_dbcti.rds')

saveRDS(l_2_dbcti, 'l_2_dbcti.rds')
saveRDS(l_3_dbcti, 'l_3_dbcti.rds')
###############################################################################################
monocle_infer <- function(data){
  gene <- as.data.frame(data[, 1])
  rownames(gene) <- rownames(data)
  colnames(gene) <- 'gene_short_name'
  cell <- as.data.frame(t(data[1,]))
  rownames(cell) <- colnames(data)
  
  
  monocle <- infer_monocle(data = data,
                           gene_annotation =  gene,
                           cell_annotation = cell,
                           gene_to_use = rownames(data))
  return(monocle)
}


#simulate and real
babelwhale::set_default_config(babelwhale::create_singularity_config())

##################################################################
simulated_data <- list(bc_0$dataset,bc_2$dataset,bc_3$dataset,bt_0$dataset,bt_2$dataset,bt_3$dataset,l_0$dataset,l_2$dataset,l_3$dataset,c_0$dataset,c_2$dataset,c_3$dataset,d_0$dataset,d_2$dataset,d_3$dataset)
names(simulated_data) <- c('bc_0','bc_2','bc_3','bt_0','bt_2','bt_3','l_0','l_2','l_3','c_0','c_2','c_3','d_0','d_2','d_3')

tool_name <- c('slingshot','paga', 'monocle2','tscan')

tool_name_list <- vector(mode = 'list', length = 5)
names(tool_name_list) <- tool_name

simulated_data_tool <- simulated_data


infer <- function(data){
  set.seed(1)
  model_slingshot <- infer_trajectory(dataset = data, ti_slingshot(), verbose = TRUE)
  set.seed(1)
  model_paga <- infer_trajectory(add_prior_information(dataset = data, start_id = 'cell1'), ti_paga(filter_features = FALSE), verbose = TRUE)
  set.seed(1)
  
  
  result <- list(model_slingshot,model_paga)
  names(result) <- c('model_slingshot','model_paga')
  return(result)
}

simulated_data_tool$bc_0 <- infer(simulated_data$bc_0)
simulated_data_tool$bc_2 <- infer(simulated_data$bc_2)
simulated_data_tool$bc_3 <- infer(simulated_data$bc_3)
simulated_data_tool$bt_0 <- infer(simulated_data$bt_0)
simulated_data_tool$bt_2 <- infer(simulated_data$bt_2)
simulated_data_tool$bt_3 <- infer(simulated_data$bt_3)
simulated_data_tool$c_0 <- infer(simulated_data$c_0)
simulated_data_tool$c_2 <- infer(simulated_data$c_2)
simulated_data_tool$c_3 <- infer(simulated_data$c_3)
simulated_data_tool$l_0 <- infer(simulated_data$l_0)
simulated_data_tool$l_2 <- infer(simulated_data$l_2)
simulated_data_tool$l_3 <- infer(simulated_data$l_3)
simulated_data_tool$d_0 <- infer(simulated_data$d_0)
simulated_data_tool$d_2 <- infer(simulated_data$d_2)
simulated_data_tool$d_3 <- infer(simulated_data$d_3)


##########################################################################

infer_1 <- function(data){
  
  set.seed(1)
  model_tscan <- infer_trajectory(dataset = data, ti_tscan(), verbose = TRUE)
  
  result <- list(model_tscan)
  names(result) <- c('model_tscan')
  return(result)
}

simulated_data_tool_2 <- simulated_data

simulated_data_tool_2$bc_0 <- infer_1(simulated_data$bc_0)
simulated_data_tool_2$bc_2 <- infer_1(simulated_data$bc_2)
simulated_data_tool_2$bc_3 <- infer_1(simulated_data$bc_3)
simulated_data_tool_2$bt_0 <- infer_1(simulated_data$bt_0)
simulated_data_tool_2$bt_2 <- infer_1(simulated_data$bt_2)
simulated_data_tool_2$bt_3 <- infer_1(simulated_data$bt_3)
simulated_data_tool_2$c_0 <- infer_1(simulated_data$c_0)
simulated_data_tool_2$c_2 <- infer_1(simulated_data$c_2)
simulated_data_tool_2$c_3 <- infer_1(simulated_data$c_3)
simulated_data_tool_2$l_0 <- infer_1(simulated_data$l_0)
simulated_data_tool_2$l_2 <- infer_1(simulated_data$l_2)
simulated_data_tool_2$l_3 <- infer_1(simulated_data$l_3)
simulated_data_tool_2$d_0 <- infer_1(simulated_data$d_0)
simulated_data_tool_2$d_2 <- infer_1(simulated_data$d_2)
simulated_data_tool_2$d_3 <- infer_1(simulated_data$d_3)







#######################################################################################################
#real
infer_real <- function(data){
  data <- wrap_expression(expression = t(data),
                          counts = t(data))
  set.seed(1)
  model_slingshot <- infer_trajectory(dataset = data, ti_slingshot(), verbose = TRUE)
  set.seed(1)
  model_paga <- infer_trajectory(add_prior_information(dataset = data, start_id = data$cell_ids[1]), ti_paga(filter_features = FALSE), verbose = TRUE)
  
  result <- list(model_slingshot,model_paga)
  names(result) <- c('model_slingshot','model_paga')
  return(result)
}

infer_1_real <- function(data){
  data <- wrap_expression(expression = t(data),
                          counts = t(data))
  
  set.seed(1)
  model_tscan <- infer_trajectory(dataset = data, ti_tscan(), verbose = TRUE)
  
  result <- list(model_tscan)
  names(result) <- c('model_tscan')
  return(result)
}





real_data <- list(as.matrix(leng_data),as.matrix(kowalczyk_y_data),as.matrix(kowalczyk_o_data),camp_data,nestorowa_data,yan_data)
names(real_data) <- c('leng_data','kowalczyk_y_data','kowalczyk_o_data','camp_data','nestorowa_data','yan_data')

real_data_tool <- real_data
real_data_tool_2 <- real_data

real_data_tool$leng_data <- infer_real(real_data$leng_data)
real_data_tool$kowalczyk_y_data <- infer_real(real_data$kowalczyk_y_data)
real_data_tool$kowalczyk_o_data <- infer_real(real_data$kowalczyk_o_data)
real_data_tool$camp_data <- infer_real(real_data$camp_data)
real_data_tool$nestorowa_data <- infer_real(real_data$nestorowa_data)
real_data_tool$yan_data <- infer_real(real_data$yan_data)



#
real_data_tool_2$leng_data <- infer_1_real(real_data$leng_data)
real_data_tool_2$kowalczyk_y_data <- infer_1_real(real_data$kowalczyk_y_data)
real_data_tool_2$kowalczyk_o_data <- infer_1_real(real_data$kowalczyk_o_data)
real_data_tool_2$camp_data <- infer_1_real(real_data$camp_data)
real_data_tool_2$nestorowa_data <- infer_1_real(real_data$nestorowa_data)
real_data_tool_2$yan_data <- infer_1_real(real_data$yan_data)

