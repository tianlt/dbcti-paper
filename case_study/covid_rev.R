library(data.table)
library(Seurat)
library(tidyverse)
library(dbcti)
setwd('..')


covid_combined.nc <- readRDS(url("https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds"))
View(covid_combined.nc)
covid_combined.nc@meta.data[["Donor"]]
covid_combined.nc@meta.data[["cell.type"]]
unique(covid_combined.nc@meta.data[["Donor"]])
unique(covid_combined.nc@meta.data[["cell.type"]])
unique(covid_combined.nc@meta.data[["Status"]])
unique(covid_combined.nc@meta.data[["cell.type.fine"]])
unique(covid_combined.nc@meta.data[["cell.type.coarse"]])

IgG_PB_index = covid_combined.nc@meta.data[["cell.type"]] == 'IgG PB' 
IgA_PB_index = covid_combined.nc@meta.data[["cell.type"]] == 'IgA PB' 
Neutrophil_index = covid_combined.nc@meta.data[["cell.type"]] == 'Neutrophil' 
donor = covid_combined.nc@meta.data[["Donor"]]

data_C1 <- SubsetData(covid_combined.nc, cells = (IgG_PB_index | Neutrophil_index) & donor == 'C1')
data_raw_C1 <- as.data.frame(data_C1@assays[["RNA"]]@counts)

data_C2 <- SubsetData(covid_combined.nc, cells = (IgG_PB_index | Neutrophil_index) & donor == 'C2')
data_raw_C2 <- as.data.frame(data_C2@assays[["RNA"]]@counts)


set.seed(122)
data_raw_C1 <- as.data.frame(data_raw_C1)
test_C1 <- create_object(data_raw_C1, normalized = FALSE) 
test_C1 <- normalize(test_C1, gene_cri = 1, cell_cri = 1)
test_C1 <- select_var_feature(test_C1, use_normalized_data = TRUE, n = 2000)
test_C1 <- tsneplot(test_C1, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_C1 <- contour_plot(test_C1)

test_C1 <- distribution_estimation(test_C1, ndraw = 1000, expansion = 2.5, ... =c(4,5),c(2,3),6,1) 
test_C1 <- point_possibility(test_C1, r = 5)
test_C1 <- connect_cluster(test_C1, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.05)
test_C1 <- infer_trajectory(test_C1, iter_n =150)
test_C1 <- calculate_pseudotime(test_C1, start_state_name = c('1','2'))
test_C1 <- plot_trajectory(test_C1)
plot(test_C1@trajectory_plot$plot)

eval_C1_ori <- plot_index(index = factor(t(data_C1@meta.data[["cell.type"]])), trajectory = test_C1@trajectory, connection_matrix = test_C1@connect_cluster$cluster_connection)
plot(eval_C1_ori$plot)
eval_C1_clu <- plot_index(index = factor((test_C1@distribution_estimation[["cluster_index"]])), trajectory = test_C1@trajectory, connection_matrix = test_C1@connect_cluster$cluster_connection)
plot(eval_C1_clu$plot)

Idents(data_C1) <- factor((test_C1@distribution_estimation[["cluster_index"]]))
de_C1_14<-FindMarkers(data_C1, ident.1 = 1, ident.2 = 4, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_41<-FindMarkers(data_C1, ident.1 = 4, ident.2 = 1, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_14
de_C1_41

de_C1_21<-FindMarkers(data_C1, ident.1 = 2, ident.2 = 1, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_12<-FindMarkers(data_C1, ident.1 = 1, ident.2 = 2, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_21
de_C1_12

de_C1_24<-FindMarkers(data_C1, ident.1 = 2, ident.2 = 4, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_42<-FindMarkers(data_C1, ident.1 = 4, ident.2 = 2, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_24
de_C1_42

de_C1_2_14<-FindMarkers(data_C1, ident.1 = 2, ident.2 = c(1,4), logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_14_2<-FindMarkers(data_C1, ident.1 = c(1,4), ident.2 = 2, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C1_2_14
de_C1_14_2


set.seed(122)
data_raw_C2 <- as.data.frame(data_raw_C2)
test_C2 <- create_object(data_raw_C2, normalized = FALSE) 
test_C2 <- normalize(test_C2, gene_cri = 1, cell_cri = 1)
test_C2 <- select_var_feature(test_C2, use_normalized_data = TRUE, n = 2000)
test_C2 <- tsneplot(test_C2, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_C2 <- contour_plot(test_C2)

test_C2 <- distribution_estimation(test_C2, ndraw = 1000, expansion = 2.5, ... =c(2,3),c(1,4,5),6) 
test_C2 <- point_possibility(test_C2, r = 5)
test_C2 <- connect_cluster(test_C2, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.05)
test_C2 <- infer_trajectory(test_C2, iter_n =15)
test_C2 <- calculate_pseudotime(test_C2, start_state_name = c('1','2'))
test_C2 <- plot_trajectory(test_C2)
plot(test_C2@trajectory_plot$plot)

eval_C2_ori <- plot_index(index = factor(t(data_C2@meta.data[["cell.type"]])), trajectory = test_C2@trajectory, connection_matrix = test_C2@connect_cluster$cluster_connection)
plot(eval_C2_ori$plot)
eval_C2_clu <- plot_index(index = factor((test_C2@distribution_estimation[["cluster_index"]])), trajectory = test_C2@trajectory, connection_matrix = test_C2@connect_cluster$cluster_connection)
plot(eval_C2_clu$plot)

Idents(data_C2) <- factor((test_C2@distribution_estimation[["cluster_index"]]))
de_C2_12<-FindMarkers(data_C2, ident.1 = 1, ident.2 = 2, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C2_21<-FindMarkers(data_C2, ident.1 = 2, ident.2 = 1, logfc.threshold = 0.25, only.pos = TRUE) %>% filter(p_val_adj <= 0.05)
de_C2_12
de_C2_21


#gene expression
eval_C1_IGLC2 <- plot_expression(gene_exp = scale_exp(as.numeric(test_C1@normalized_data['IGLC2', ][1,])), trajectory = test_C1@trajectory, connection_matrix = test_C1@connect_cluster$cluster_connection)
plot(eval_C1_IGLC2$plot)
eval_C1_IGLC3 <- plot_expression(gene_exp = scale_exp(as.numeric(test_C1@normalized_data['IGLC3', ][1,])), trajectory = test_C1@trajectory, connection_matrix = test_C1@connect_cluster$cluster_connection)
plot(eval_C1_IGLC3$plot)
eval_C1_IGKC <- plot_expression(gene_exp = scale_exp(as.numeric(test_C1@normalized_data['IGKC', ][1,])), trajectory = test_C1@trajectory, connection_matrix = test_C1@connect_cluster$cluster_connection)
plot(eval_C1_IGKC$plot)

eval_C2_IGLC2 <- plot_expression(gene_exp = scale_exp(as.numeric(test_C2@normalized_data['IGLC2', ][1,])), trajectory = test_C2@trajectory, connection_matrix = test_C2@connect_cluster$cluster_connection)
plot(eval_C2_IGLC2$plot)
eval_C2_IGLC3 <- plot_expression(gene_exp = scale_exp(as.numeric(test_C2@normalized_data['IGLC3', ][1,])), trajectory = test_C2@trajectory, connection_matrix = test_C2@connect_cluster$cluster_connection)
plot(eval_C2_IGLC3$plot)
eval_C2_IGKC <- plot_expression(gene_exp = scale_exp(as.numeric(test_C2@normalized_data['IGKC', ][1,])), trajectory = test_C2@trajectory, connection_matrix = test_C2@connect_cluster$cluster_connection)
plot(eval_C2_IGKC$plot)


gene_diff <- c('IGLC2','IGLC3', 'IGKC')

# de_C1_24[gene_diff,]#the cluster not connect 'IGLC2','IGLC3'
# de_C1_42[gene_diff,]#'IGKC'
# 
# de_C1_21[gene_diff,]#the cluster not connect 'IGLC2','IGLC3'
# de_C1_12[gene_diff,]#'IGKC'
#
de_C1_2_14[gene_diff,]#the cluster not connect 'IGLC2','IGLC3'
de_C1_14_2[gene_diff,]#'IGKC'
#
de_C2_12[gene_diff,]#the cluster not connect 'IGLC2','IGLC3'
de_C2_21[gene_diff,]#'IGKC'


#difference in liner plsmablast
de_C1_14[gene_diff,]
de_C1_41[gene_diff,]

#plot############################
plot(eval_C1_ori$plot)
plot(eval_C2_ori$plot)
plot(eval_C1_clu$plot)
plot(eval_C2_clu$plot)

pdf(file = 'eval_C1_ori.pdf',width = 6,height = 6)
plot(eval_C1_ori$plot)
dev.off()

pdf(file = 'eval_C2_ori.pdf',width = 6,height = 6)
plot(eval_C2_ori$plot)
dev.off()

pdf(file = 'eval_C1_clu.pdf',width = 6,height = 6)
plot(eval_C1_clu$plot)
dev.off()

pdf(file = 'eval_C2_clu.pdf',width = 6,height = 6)
plot(eval_C2_clu$plot)
dev.off()



pdf(file = 'eval_C1_IGLC2.pdf',width = 6,height = 6)
plot(eval_C1_IGLC2$plot)
dev.off()

pdf(file = 'eval_C1_IGLC3.pdf',width = 6,height = 6)
plot(eval_C1_IGLC3$plot)
dev.off()

pdf(file = 'eval_C1_IGKC.pdf',width = 6,height = 6)
plot(eval_C1_IGKC$plot)
dev.off()

pdf(file = 'eval_C2_IGLC2.pdf',width = 6,height = 6)
plot(eval_C2_IGLC2$plot)
dev.off()

pdf(file = 'eval_C2_IGLC3.pdf',width = 6,height = 6)
plot(eval_C2_IGLC3$plot)
dev.off()

pdf(file = 'eval_C2_IGKC.pdf',width = 6,height = 6)
plot(eval_C2_IGKC$plot)
dev.off()
##############################
plot_index <- function(index, trajectory, connection_matrix, alpha = 1, jitter = 0.2, plot_title = ''){
  
  ordered_trajectory<-trajectory$result_list[[1]]
  for (i in trajectory$result_list[-1]) ordered_trajectory<-rbind(as.matrix(i), as.matrix(ordered_trajectory))
  
  ordered_trajectory<-ordered_trajectory[order(as.numeric(rownames(ordered_trajectory))),]
  
  gg_data<-as.data.frame(ordered_trajectory)
  g<-ggplot()
  
  for (i in trajectory[["lines_list"]]) {
    line <- as.data.frame(i$s[i$ord,])
    colnames(line) <- c('x', 'y')
    g <- g + geom_path(data = line, aes(x,y), size = 1.5)
  }
  
  g <- g + geom_point(data = gg_data, aes(x,y, fill = index), color = 'black', size=15, pch =21, alpha = alpha, stroke = 0.5, position = position_jitter(width = jitter, height = jitter)) +
    scale_fill_brewer(name = 'Index', palette = 'Dark2') + theme_classic()+xlab('x')+ylab('y') + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                             axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                                             axis.title.y=element_blank(), legend.title=element_text(size=25, face="bold"),legend.text=element_text(size=20))
  
  if (length(plot_title) >= 1) g <- g + ggtitle(label = plot_title) 
  result = list(ordered_trajectory = ordered_trajectory, plot = g)
  return(result)
}

plot_expression <- function(gene_exp, trajectory, connection_matrix, alpha = 1, jitter = 0.2, plot_title = ''){
  
  ordered_trajectory<-trajectory$result_list[[1]]
  for (i in trajectory$result_list[-1]) ordered_trajectory<-rbind(as.matrix(i), as.matrix(ordered_trajectory))
  
  ordered_trajectory<-ordered_trajectory[order(as.numeric(rownames(ordered_trajectory))),]
  
  gg_data<-as.data.frame(ordered_trajectory)
  g<-ggplot()
  
  for (i in trajectory[["lines_list"]]) {
    line <- as.data.frame(i$s[i$ord,])
    colnames(line) <- c('x', 'y')
    g <- g + geom_path(data = line, aes(x,y), size = 1.5)
  }
  exp <- gene_exp
  g <- g + geom_point(data = gg_data, aes(x,y, fill = exp), size=15, pch =21, alpha = alpha, stroke = 0.5, position = position_jitter(width = jitter, height = jitter)) +
    scale_fill_gradient(low = "white", high = "black") + theme_classic()+xlab('x')+ylab('y') + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                           axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                                           axis.title.y=element_blank(), legend.title=element_text(size=25, face="bold"),legend.text=element_text(size=20)) + labs('exp')
  
  if (length(plot_title) >= 1) g <- g + ggtitle(label = plot_title) 
  result = list(ordered_trajectory = ordered_trajectory, plot = g)
  return(result)
}



scale_exp <- function(x){(x-min(x))/(max(x)-min(x))}

