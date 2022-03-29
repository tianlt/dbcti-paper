# setwd('..')
library(monocle)
library(TSCAN)
library(SCORPIUS)
library(igraph)
library(ggplot2)
library(NbClust)
library(dbcti)
path = ''
#dataset######################
#real dataset 
test_hescmt <- readRDS('datasets/test_hescmt.rds')
test_stem_mouse_C57BL6_young<- readRDS('datasets/test_stem_mouse_C57BL6_young.rds')
test_stem_mouse_C57BL6_old<- readRDS('datasets/test_stem_mouse_C57BL6_old.rds')
test_camp1_batch1<- readRDS('datasets/test_camp1_batch1.rds')
test_nestorowa_batch1_type_filter_no_stem<- readRDS('datasets/test_nestorowa_batch1_type_filter_no_stem.rds')
test_yan<- readRDS('datasets/test_yan.rds')

hesc_monocle<- readRDS('datasets/hesc_monocle.rds')
hesc_tscan_c<- readRDS('datasets/hesc_tscan_c.rds')
hesc_tscan_o<- readRDS('datasets/hesc_tscan_o.rds')
hesc_sco_s<- readRDS('datasets/hesc_sco_s.rds')
hesc_sco_t<- readRDS('datasets/hesc_sco_t.rds')

stem_mouse_C57BL6_data_young_monocle<- readRDS('datasets/stem_mouse_C57BL6_data_young_monocle.rds')
stem_mouse_C57BL6_data_young_tscan_c<- readRDS('datasets/stem_mouse_C57BL6_data_young_tscan_c.rds')
stem_mouse_C57BL6_data_young_tscan_o<- readRDS('datasets/stem_mouse_C57BL6_data_young_tscan_o.rds')
stem_mouse_C57BL6_data_young_s<- readRDS('datasets/stem_mouse_C57BL6_data_young_s.rds')
stem_mouse_C57BL6_data_young_t<- readRDS('datasets/stem_mouse_C57BL6_data_young_t.rds')

stem_mouse_C57BL6_data_old_monocle<- readRDS('datasets/stem_mouse_C57BL6_old_monocle.rds')
stem_mouse_C57BL6_data_old_tscan_c<- readRDS('datasets/stem_mouse_C57BL6_data_old_tscan_c.rds')
stem_mouse_C57BL6_data_old_tscan_o<- readRDS('datasets/stem_mouse_C57BL6_data_old_tscan_o.rds')
stem_mouse_C57BL6_data_old_s<- readRDS('datasets/stem_mouse_C57BL6_data_old_s.rds')
stem_mouse_C57BL6_data_old_t<- readRDS('datasets/stem_mouse_C57BL6_data_old_t.rds')

camp1_data_batch1_monocle<- readRDS('datasets/camp1_data_batch1_monocle.rds')
camp1_data_batch1_tscan_c<- readRDS('datasets/camp1_data_batch1_tscan_c.rds')
camp1_data_batch1_tscan_o<- readRDS('datasets/camp1_data_batch1_tscan_o.rds')
camp1_data_batch1_s<- readRDS('datasets/camp1_data_batch1_s.rds')
camp1_data_batch1_t<- readRDS('datasets/camp1_data_batch1_t.rds')

nestorowa_data_batch1_type_filter_no_stem_monocle<- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_monocle.rds')
nestorowa_data_batch1_type_filter_no_stem_c<- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_c.rds')
nestorowa_data_batch1_type_filter_no_stem_o<- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_o.rds')
nestorowa_data_batch1_type_filter_no_stem_s<- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_s.rds')
nestorowa_data_batch1_type_filter_no_stem_t<- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_t.rds')

yan_data_monocle<- readRDS('datasets/yan_data_monocle.rds')
yan_data_c<- readRDS('datasets/yan_data_c.rds')
yan_data_o<- readRDS('datasets/yan_data_o.rds')
yan_data_s<- readRDS('datasets/yan_data_c.rds')
yan_data_t<- readRDS('datasets/yan_data_t.rds')

#simulated dataset
bifurcating_cycle_monocle<- readRDS('datasets/bifurcating_cycle_monocle.rds')
bifurcating_cycle_tscan_c<- readRDS('datasets/bifurcating_cycle_tscan_c.rds')
bifurcating_cycle_tscan_o<- readRDS('datasets/bifurcating_cycle_tscan_o.rds')
bifurcating_cycle_sco_s<- readRDS('datasets/bifurcating_cycle_sco_s.rds')
bifurcating_cycle_sco_t<- readRDS('datasets/bifurcating_cycle_sco_t.rds')
test_bifurcating_cycle<- readRDS('datasets/test_bifurcating_cycle.rds')

cycle_monocle<- readRDS('datasets/cycle_monocle.rds')
cycle_tscan_c<- readRDS('datasets/cycle_tscan_c.rds')
cycle_tscan_o<- readRDS('datasets/cycle_tscan_o.rds')
cycle_sco_s<- readRDS('datasets/cycle_sco_s.rds')
cycle_sco_t<- readRDS('datasets/cycle_sco_t.rds')
test_cycle<- readRDS('datasets/test_cycle.rds')

disconnected_monocle<- readRDS('datasets/disconnected_monocle.rds')
disconnected_tscan_c<- readRDS('datasets/disconnected_tscan_c.rds')
disconnected_tscan_o<- readRDS('datasets/disconnected_tscan_o.rds')
disconnected_sco_s<- readRDS('datasets/disconnected_sco_s.rds')
disconnected_sco_t<- readRDS('datasets/disconnected_sco_t.rds')
test_disconnected<- readRDS('datasets/test_disconnected.rds')

binary_tree_monocle<- readRDS('datasets/binary_tree_monocle.rds')
binary_tree_tscan_c<- readRDS('datasets/binary_tree_tscan_c.rds')
binary_tree_tscan_o<- readRDS('datasets/binary_tree_tscan_o.rds')
binary_tree_sco_s<- readRDS('datasets/binary_tree_sco_s.rds')
binary_tree_sco_t<- readRDS('datasets/binary_tree_sco_t.rds')
test_binary_tree<- readRDS('datasets/test_binary_tree.rds')

linear_monocle<- readRDS('datasets/test_binary_tree.rds')
linear_tscan_c<- readRDS('datasets/linear_tscan_c.rds')
linear_tscan_o<- readRDS('datasets/linear_tscan_o.rds')
linear_sco_s<- readRDS('datasets/linear_sco_s.rds')
linear_sco_t<- readRDS('datasets/linear_sco_t.rds')
test_linear<- readRDS('datasets/test_linear.rds')


#metric ################
#structure 1
eval_topo <- read.csv('datasets/performance_evaluation.csv', stringsAsFactors = FALSE)

simulated_relative_score <- ggplot(subset(eval_topo, dataset_index == 'simulated'),aes(x = dataset, y = score_1, fill = tool))+
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') +
  theme_classic() +
  xlab('dataset')+ylab('score')+
  ggtitle(label = 'Relative topological similarity score')+
  theme(axis.text.x = element_text(face = 'bold', size = 12, angle = -5), axis.title = element_text(face = 'bold', size = 15), title = element_text(face = 'bold'))

plot(simulated_relative_score)

real_relative_score <- ggplot(subset(eval_topo, dataset_index == 'real'),aes(x = dataset, y = score_1, fill = tool))+
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') +
  theme_classic() +
  xlab('dataset')+ylab('score')+
  ggtitle(label = 'Relative topological similarity score') +
  theme(axis.text.x = element_text(face = 'bold', size = 12, angle = -5), axis.title = element_text(face = 'bold', size = 15), title = element_text(face = 'bold'))

plot(real_relative_score)

#structure 2
eval_topo_2 <- read.csv('datasets/performance_evaluation_2.csv', stringsAsFactors = FALSE)
eval_topo_2 = eval_topo_2 [1:16, 1:3]
pattern_matching_score <- ggplot(eval_topo_2, aes(x = structure_reference, y = score, fill = tool))+
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black') +
  theme_classic() +
  xlab('dataset')+ylab('score')+
  ggtitle(label = 'Pattern matching score') +
  theme(axis.text.x = element_text(face = 'bold', size = 12, angle = -5), axis.title = element_text(face = 'bold', size = 15), title = element_text(face = 'bold'))

plot(pattern_matching_score)

#classification purity
#real
hesc_purity = calculate_purity(DBCTI = test_hescmt, monocle = hesc_monocle, tscan_c = hesc_tscan_c, sco_s = hesc_sco_s)
camp_purity = calculate_purity(DBCTI = test_camp1_batch1, monocle =camp1_data_batch1_monocle, tscan_c = camp1_data_batch1_tscan_c, sco_s = camp1_data_batch1_s)
yan_purity = calculate_purity(DBCTI = test_yan, monocle = yan_data_monocle, tscan_c = yan_data_c, sco_s = yan_data_s)
stem_mouse_young_purity = calculate_purity(DBCTI = test_stem_mouse_C57BL6_young, monocle = stem_mouse_C57BL6_data_young_monocle, tscan_c = stem_mouse_C57BL6_data_young_tscan_c, sco_s = stem_mouse_C57BL6_data_young_s)
stem_mouse_old_purity = calculate_purity(DBCTI = test_stem_mouse_C57BL6_old, monocle = stem_mouse_C57BL6_data_old_monocle, tscan_c = stem_mouse_C57BL6_data_old_tscan_c, sco_s = stem_mouse_C57BL6_data_old_s)
nestorowa_purity = calculate_purity(DBCTI = test_nestorowa_batch1_type_filter_no_stem, monocle = nestorowa_data_batch1_type_filter_no_stem_monocle, tscan_c = nestorowa_data_batch1_type_filter_no_stem_c, sco_s = nestorowa_data_batch1_type_filter_no_stem_s)

#to data frame
real_purity_score = cbind(hesc_purity, camp_purity, yan_purity, stem_mouse_young_purity, stem_mouse_old_purity, nestorowa_purity)
real_purity_score[2, ] = c(rep('Leng_data', 4), rep('Camp_data', 4), rep('Yan_data', 4), rep('Kowalczyk_young_data', 4), rep('Kowalczyk_old_data', 4), rep('Nestorowa_data', 4))
real_purity_score[3, ] = unlist(lapply(strsplit(colnames(real_purity_score), '_'), function(x) x[1]))
real_purity_score[4, ] = rep('real', 24)
rownames(real_purity_score) = c('score', 'dataset', 'tool', 'dataset_index')

#simulated
linear_purity = calculate_purity(DBCTI = test_linear, monocle = linear_monocle, tscan_c = linear_tscan_c, sco_s = linear_sco_s)
cycle_purity = calculate_purity(DBCTI = test_cycle, monocle = cycle_monocle, tscan_c = cycle_tscan_c, sco_s = cycle_sco_s)
binary_tree_purity = calculate_purity(DBCTI = test_binary_tree, monocle = binary_tree_monocle, tscan_c = binary_tree_tscan_c, sco_s = binary_tree_sco_s)
disconnected_purity = calculate_purity(DBCTI = test_disconnected, monocle = disconnected_monocle, tscan_c = disconnected_tscan_c, sco_s = disconnected_sco_s)
bifurcating_cycle_purity = calculate_purity(DBCTI = test_bifurcating_cycle, monocle = bifurcating_cycle_monocle, tscan_c = bifurcating_cycle_tscan_c, sco_s = bifurcating_cycle_sco_s)

#to data frame
simulated_purity_score = cbind(linear_purity, cycle_purity, binary_tree_purity, disconnected_purity, bifurcating_cycle_purity)
simulated_purity_score[2, ] = c(rep('linear_simulated', 4), rep('cycle_simulated', 4), rep('binary_tree_simulated', 4), rep('disconnected_simulated', 4), rep('bifurcating_cycle_simulated', 4))
simulated_purity_score[3, ] = unlist(lapply(strsplit(colnames(simulated_purity_score), '_'), function(x) x[1]))
simulated_purity_score[4, ] = rep('simulated', 20)
rownames(simulated_purity_score) = c('score', 'dataset', 'tool', 'dataset_index')

#plot purity score
purity_score = data.frame(t(cbind(simulated_purity_score, real_purity_score)), stringsAsFactors = FALSE)
purity_score[, 1] = round(as.numeric(purity_score[, 1]), 2)


simulated_purity_score_plot <- ggplot(subset(purity_score, dataset_index == 'simulated'),aes(x = dataset, y = score, fill = tool, color = tool))+
  geom_point(size=3) + geom_line(aes(group = tool), size = 1) +
  theme_classic() + 
  xlab('dataset') + ylab('score') +
  ggtitle(label = 'Purity score on simulated dataset line plot') +
  theme(axis.text.x = element_text(face = 'bold', size = 12, angle = -5), axis.title = element_text(face = 'bold', size = 15), title = element_text(face = 'bold'))

plot(simulated_purity_score_plot)

real_purity_score_plot <- ggplot(subset(purity_score, dataset_index == 'real'),aes(x = dataset, y = score, fill = tool, color = tool))+
  geom_point(size=3) + geom_line(aes(group = tool), size = 1) +
  theme_classic() +
  xlab('dataset') + ylab('score') +
  ggtitle(label = 'Purity score on real dataset line plot') +
  theme(axis.text.x = element_text(face = 'bold', size = 12, angle = -5), axis.title = element_text(face = 'bold', size = 15), title = element_text(face = 'bold'))

plot(real_purity_score_plot)

purity_score_plot <- ggplot(purity_score, aes(x = dataset, y = score, fill = tool, color = tool))+
  geom_point(size=3) + geom_line(aes(group = tool), size = 2) +
  theme_classic() +
  xlab('dataset') + ylab('score') +
  ggtitle(label = 'Purity score line plot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22, angle = 0, vjust = 0.5), axis.text.y = element_text(face = 'bold', size = 16), axis.title = element_text(face = 'bold', size = 22), title = element_text(face = 'bold', size = 22),
        , legend.title=element_text(size=25, face="bold"),legend.text=element_text(size=25)) + coord_flip()

plot(purity_score_plot)

purity_score_boxplot <- ggplot(purity_score, aes(x = tool, y = score, fill = tool), color = 'black')+
  geom_boxplot() + 
  theme_classic() +
  xlab('dataset') + ylab('score') +
  ggtitle(label = 'Purity score boxplot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22, angle = 0), axis.text.y = element_text(size = 22, angle = 0), axis.title = element_text(face = 'bold', size = 22), 
        title = element_text(face = 'bold', size = 22), legend.title=element_text(size=25, face="bold"),legend.text=element_text(size=25))

plot(purity_score_boxplot)


#
pdf(file = paste0(path, 'simulated_relative_score.pdf'),width = 16,height = 6)
plot(simulated_relative_score)
dev.off()

pdf(file = paste0(path, 'real_relative_score.pdf'),width = 16,height = 6)
plot(real_relative_score)
dev.off()

pdf(file = paste0(path, 'pattern_matching_score.pdf'),width = 16,height = 6)
plot(pattern_matching_score)
dev.off()

pdf(file = paste0(path, 'purity_score_plot.pdf'),width = 10,height = 6)
plot(purity_score_plot)
dev.off()

pdf(file = paste0(path, 'purity_score_boxplot.pdf'),width = 10,height = 6)
plot(purity_score_boxplot)
dev.off()



#function
calculate_purity <- function(DBCTI, monocle, tscan_c, sco_s){
  DBCTI_purity <- NMF::purity(as.factor(DBCTI@distribution_estimation[["cluster_index"]]), monocle@phenoData@data[["index"]])
  
  monocle_purity <- NMF::purity(monocle@phenoData@data[["State"]], monocle@phenoData@data[["index"]])
  
  tscan_purity <- NMF::purity(as.factor(tscan_c[["clusterid"]]), monocle@phenoData@data[["index"]])
  
  sco_pos = sco_s[, c(1, 2)]
  sco_kmeans = NbClust::NbClust(sco_pos, method = 'kmeans')$Best.partition
  sco_purity <- NMF::purity(as.factor(sco_kmeans), monocle@phenoData@data[["index"]])
  res = data.frame(monocle_purity, tscan_purity, sco_purity, DBCTI_purity)
  
  return(res)
}

