library(dbcti)
library(tidyverse)
library(dyngen)
library(dyno)

# path = 'Data/cell_trajectory/trajectory_inference/plot/simulated_dataset_comparison/'
setwd('..')
path = ''
#backbone_bifurcating_cycle###############
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

dataset_bifurcating_cycle <- out_bifurcating_cycle$dataset
model_bifurcating_cycle <- out_bifurcating_cycle$model
print(out_bifurcating_cycle$plot)

pdf(file = paste0(path, 'dataset_bifurcating_cycle.pdf'),width = 6,height = 6)
plot_dimred(dataset_bifurcating_cycle, label_milestones = TRUE)
dev.off()

plot_graph(dataset_bifurcating_cycle)
############### rds
dataset_bifurcating_cycle_df = as.data.frame(t(as.matrix(dataset_bifurcating_cycle[["expression"]]))) #dataset matrix
bifurcating_cycle_data = dataset_bifurcating_cycle_df[apply(dataset_bifurcating_cycle_df, 1, sum) != 0, ]
bifurcating_cycle_index = dataset_bifurcating_cycle[["progressions"]] #cell info
bifurcating_cycle_loc_dimred = dataset_bifurcating_cycle[["dimred"]] #cell location

saveRDS(bifurcating_cycle_data, file = 'datasets/bifurcating_cycle_data.rds')
saveRDS(bifurcating_cycle_index, file = 'datasets/bifurcating_cycle_index.rds')
saveRDS(bifurcating_cycle_loc_dimred, file = 'datasets/bifurcating_cycle_loc_dimred.rds')





#backbone_cycle###############
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

dataset_cycle <- out_cycle$dataset
model_cycle <- out_cycle$model
print(out_cycle$plot)


pdf(file = paste0(path, 'dataset_cycle.pdf'),width = 6,height = 6)
plot_dimred(dataset_cycle, label_milestones = TRUE)
dev.off()
plot_graph(dataset_cycle)

############### rds
dataset_cycle_df = as.data.frame(t(as.matrix(dataset_cycle[["expression"]]))) #dataset matrix
cycle_data = dataset_cycle_df[apply(dataset_cycle_df, 1, sum) != 0, ]

cycle_index = dataset_cycle[["progressions"]] #cell info
cycle_loc_dimred = dataset_cycle[["dimred"]] #cell location

saveRDS(cycle_data, file = 'datasets/cycle_data.rds')
saveRDS(cycle_index, file = 'datasets/cycle_index.rds')
saveRDS(cycle_loc_dimred, file = 'datasets/cycle_loc_dimred.rds')

#backbone_disconnected###############
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

dataset_disconnected <- out_disconnected$dataset
model_disconnected <- out_disconnected$model
print(out_disconnected$plot)

pdf(file = paste0(path, 'dataset_disconnected.pdf'),width = 6,height = 6)
plot_dimred(dataset_disconnected, label_milestones = TRUE)
dev.off()
plot_graph(dataset_disconnected)

############### rds
dataset_disconnected_df = as.data.frame(t(as.matrix(dataset_disconnected[["expression"]]))) #dataset matrix
disconnected_data = dataset_disconnected_df[apply(dataset_disconnected_df, 1, sum) != 0, ]

disconnected_index = dataset_disconnected[["progressions"]] #cell info
disconnected_loc_dimred = dataset_disconnected[["dimred"]] #cell location

saveRDS(disconnected_data, file = 'dataset/disconnected_data.rds')
saveRDS(disconnected_index, file = 'dataset/disconnected_index.rds')
saveRDS(disconnected_loc_dimred, file = 'dataset/disconnected_loc_dimred.rds')

#backbone_binary_tree###############
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

dataset_binary_tree <- out_binary_tree$dataset
model_binary_tree <- out_binary_tree$model
print(out_binary_tree$plot)

pdf(file = paste0(path, 'dataset_binary_tree.pdf'),width = 6,height = 6)
plot_dimred(dataset_binary_tree, label_milestones = TRUE)
dev.off()
plot_graph(dataset_binary_tree)

############### rds
dataset_binary_tree_df = as.data.frame(t(as.matrix(dataset_binary_tree[["expression"]]))) #dataset matrix
binary_tree_data = dataset_binary_tree_df[apply(dataset_binary_tree_df, 1, sum) != 0, ]

binary_tree_index = dataset_binary_tree[["progressions"]] #cell info
binary_tree_loc_dimred = dataset_binary_tree[["dimred"]] #cell location

saveRDS(binary_tree_data, file = 'dataset/binary_tree_data.rds')
saveRDS(binary_tree_index, file = 'dataset/binary_tree_index.rds')
saveRDS(binary_tree_loc_dimred, file = 'dataset/binary_tree_loc_dimred.rds')




#backbone_linear###############
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
out_linear$plot


dataset_linear <- out_linear$dataset
model_linear <- out_linear$model
pdf(file = paste0(path, 'dataset_linear.pdf'),width = 6,height = 6)
plot_dimred(dataset_linear, label_milestones = TRUE)
dev.off()
plot_graph(dataset_linear)
############### rds
dataset_linear_df = as.data.frame(t(as.matrix(dataset_linear[["expression"]]))) #dataset matrix
linear_data = dataset_linear_df[apply(dataset_linear_df, 1, sum) != 0, ]
linear_index = dataset_linear[["progressions"]] #cell info
linear_loc_dimred = dataset_linear[["dimred"]] #cell location

saveRDS(linear_data, file = 'dataset/linear_data.rds')
saveRDS(linear_index, file = 'dataset/linear_index.rds')
saveRDS(linear_loc_dimred, file = 'dataset/linear_loc_dimred.rds')
