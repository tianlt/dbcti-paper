library(dbcti)
source('simulation_and_real_data/plot_tools_function_upload.R')
library(ggplot2)
library(tidyverse)
path = ''
set.seed(122)

#data
darmanisnormalizeddata<-read.csv('GBM_normalized_gene_counts.csv',sep = '')
darmanisrawdata<-read.csv('GBM_raw_gene_counts.csv',sep = '')
darmanismetadata<-read.csv('GBM_metadata.csv',sep = '')

cellloc<-darmanismetadata$Location
patient<-darmanismetadata$Sample.name
celltype<-darmanismetadata$Selection

bts2data<-darmanisrawdata[,cellloc=='Tumor'&patient=='BT_S2'] #723


#bts2data##########################
set.seed(122)
bts2data <- as.data.frame(bts2data)
test_bts2data <- create_object(bts2data, normalized = FALSE) 
test_bts2data <- normalize(test_bts2data, gene_cri = 1, cell_cri = 1)
test_bts2data <- select_var_feature(test_bts2data, use_normalized_data = TRUE, n = 2000)
test_bts2data <- tsneplot(test_bts2data, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)
test_bts2data <- contour_plot(test_bts2data)
test_bts2data <- distribution_estimation(test_bts2data, ndraw = 1000, expansion = 2.5, ... = c(6,9),c(4,2,7),3,5,8,1) 
test_bts2data <- point_possibility(test_bts2data, r = 5)
test_bts2data <- connect_cluster_vague_name(test_bts2data, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.05)
test_bts2data <- infer_trajectory(test_bts2data, iter_n =150)
test_bts2data <- calculate_pseudotime(test_bts2data, start_state_name = c('1','4','5'))
test_bts2data <- plot_trajectory_vague_name(test_bts2data)
plot(test_bts2data@trajectory_plot$plot)

#evaluation
eval_bts2data_clu <- plot_index(index = factor(t(test_bts2data@distribution_estimation[["cluster_index"]])), trajectory = test_bts2data@trajectory, connection_matrix = test_bts2data@connect_cluster$cluster_connection)
plot(eval_bts2data_clu$plot)


#de
library(Seurat)
bts2_seurat <- CreateSeuratObject(counts = bts2data)
bts2_seurat <- NormalizeData(bts2_seurat, normalization.method = "LogNormalize")
bts2_seurat <- ScaleData(bts2_seurat)
bts2_seurat <- FindVariableFeatures(bts2_seurat, selection.method = "vst", nfeatures = 2000)
bts2_seurat <- RunPCA(bts2_seurat, do.print=FALSE)
Idents(bts2_seurat) <- factor((test_bts2data@distribution_estimation[["cluster_index"]]))

de_bts2_all<-FindAllMarkers(bts2_seurat, logfc.threshold = 0.1)

de_bts2_1<-FindMarkers(bts2_seurat, ident.1 = c(1),logfc.threshold = 0.1, only.pos = TRUE)

de_bts2_all_tb_pos <- de_bts2_all %>% as_tibble() %>% group_by(cluster) %>% filter(avg_logFC > 0.2 & p_val_adj < 0.05) 
de_bts2_all_tb_neg <- de_bts2_all %>% as_tibble() %>% group_by(cluster) %>% filter(avg_logFC < -0.2 & p_val_adj < 0.05) 
de_bts2_bk_gene <- rownames(bts2data)

de_bts5_2<-FindMarkers(bts2_seurat, ident.1 = 5, ident.2 = 2, logfc.threshold = 0.1, only.pos = TRUE)
de_bts2_6<-FindMarkers(bts2_seurat, ident.1 = 2, ident.2 = 6, logfc.threshold = 0.1, only.pos = TRUE)
de_bts6_3<-FindMarkers(bts2_seurat, ident.1 = 6, ident.2 = 3, logfc.threshold = 0.1, only.pos = TRUE)
de_bts5_6<-FindMarkers(bts2_seurat, ident.1 = 5, ident.2 = 6, logfc.threshold = 0.1, only.pos = TRUE)

de_bts5_2_tb <- de_bts5_2 %>% as_tibble(rownames = 'gene')  %>% filter(avg_logFC > 0.2 & p_val_adj < 0.05) 
de_bts2_6_tb <- de_bts2_6 %>% as_tibble(rownames = 'gene')  %>% filter(avg_logFC > 0.2 & p_val_adj < 0.05) 
de_bts6_3_tb <- de_bts6_3 %>% as_tibble(rownames = 'gene')  %>% filter(avg_logFC > 0.2 & p_val_adj < 0.05) 
de_bts5_6_tb <- de_bts5_6 %>% as_tibble(rownames = 'gene')  %>% filter(avg_logFC > 0.2 & p_val_adj < 0.05) 

write.csv(de_bts5_2_tb , file = 'de_bts5_2.csv')
write.csv(de_bts2_6_tb , file = 'de_bts2_6.csv')
write.csv(de_bts6_3_tb , file = 'de_bts6_3.csv')
write.csv(de_bts5_6_tb , file = 'de_bts5_6.csv')

#annotation
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(AnnotationDbi)
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                             keys = de_bts2_bk_gene,  # data to use for retrieval
                                             columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                             keytype = "SYMBOL") # type of data given in 'keys' argument



annotations_orgDb <- annotations_orgDb[which(duplicated(annotations_orgDb$SYMBOL) == FALSE), ]



annotation_ahb <- read.csv('datasets/annotations_ahb.txt')
de_bts2_all_tb_pos_ahb <- inner_join(de_bts2_all_tb_pos, annotation_ahb, by = c('gene'='gene_name'))
de_bts2_bk_gene <- tibble(gene = de_bts2_bk_gene)
de_bts2_bk_gene_ahb <- inner_join(de_bts2_bk_gene, annotation_ahb, by = c('gene'='gene_name'))

sig_bts2_pos_genes <- group_map(de_bts2_all_tb_pos_ahb, function(x, ...) as.character(x$gene_id))
all_bts2_genes <-  as.character(de_bts2_bk_gene_ahb$gene_id)




enrichGO_map <- function(x){
  ego<-enrichGO(gene = sig_bts2_pos_genes[[x]],
              universe = all_bts2_genes,
              keyType = 'ENSEMBL',
              OrgDb = org.Hs.eg.db,
              ont = 'ALL',
              qvalueCutoff = 0.05,
              readable = TRUE)
  ego
}
enrich_bts2_all_pos <- map(c(1:6), enrichGO_map)

enrichGO_cc_map <- function(x){
  ego<-enrichGO(gene = sig_bts2_pos_genes[[x]],
                universe = all_bts2_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'CC',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrich_bts2_all_pos_cc <- map(c(1:6), enrichGO_cc_map)


summary_enrich_bts2_all_pos <- lapply(enrich_bts2_all_pos, data.frame)
dotplot(enrich_bts2_all_pos[[1]], showCategory=20)
dotplot(enrich_bts2_all_pos[[5]], showCategory=20)

get_dotplot <- function(x, ...) {
  dotplot(enrich_bts2_all_pos[[x]], showCategory=20)
}


map(1:6,  get_dotplot, i = 20)



dotplot(enrich_bts2_all_pos[[2]], showCategory=20)
dotplot(enrich_bts2_all_pos[[3]], showCategory=20)
dotplot(enrich_bts2_all_pos[[5]], showCategory=20)
dotplot(enrich_bts2_all_pos[[6]], showCategory=20)


dotplot(enrich_bts2_all_pos_cc[[2]], showCategory=20)
dotplot(enrich_bts2_all_pos_cc[[3]], showCategory=20)
dotplot(enrich_bts2_all_pos_cc[[5]], showCategory=20)
dotplot(enrich_bts2_all_pos_cc[[6]], showCategory=20)
tibble(cluster1='microglia', cluster4='neuron')


#cell marker
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>% 
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
 
  
  
   
de_bts2_all_tb_pos_entrez <- inner_join(de_bts2_all_tb_pos, annotation_ahb, by = c('gene'='gene_name'))
de_bts2_all_tb_pos_entrez <- group_split(de_bts2_all_tb_pos_entrez)



get_cell_marker <- function(x) {
  out <- enricher(x$entrezid, TERM2GENE=cell_markers, minGSSize=1)
  out
}

de_bts2_all_enricher_out <- lapply(de_bts2_all_tb_pos_entrez, get_cell_marker)
View(as.data.frame(de_bts2_all_enricher_out[[2]]))
View(as.data.frame(de_bts2_all_enricher_out[[3]]))
View(as.data.frame(de_bts2_all_enricher_out[[5]]))
View(as.data.frame(de_bts2_all_enricher_out[[6]]))

DT::datatable(as.data.frame(de_bts2_all_enricher_out[[1]])) #microglia
DT::datatable(as.data.frame(de_bts2_all_enricher_out[[2]])) #Astrocyte lineage
DT::datatable(as.data.frame(de_bts2_all_enricher_out[[3]])) #Astrocyte lineage
DT::datatable(as.data.frame(de_bts2_all_enricher_out[[4]])) #olidendrocyte
DT::datatable(as.data.frame(de_bts2_all_enricher_out[[5]])) #Astrocyte lineage
DT::datatable(as.data.frame(de_bts2_all_enricher_out[[6]])) #Astrocyte lineage #28 de genes

de_bts2_all_enricher_out_1 <- as.data.frame(de_bts2_all_enricher_out[[1]]) #microglia
de_bts2_all_enricher_out_2 <- as.data.frame(de_bts2_all_enricher_out[[2]]) #Astrocyte lineage
de_bts2_all_enricher_out_3 <- as.data.frame(de_bts2_all_enricher_out[[3]]) #Astrocyte lineage
de_bts2_all_enricher_out_4 <- as.data.frame(de_bts2_all_enricher_out[[4]]) #olidendrocyte
de_bts2_all_enricher_out_5 <- as.data.frame(de_bts2_all_enricher_out[[5]]) #Astrocyte lineage
de_bts2_all_enricher_out_6 <- as.data.frame(de_bts2_all_enricher_out[[6]]) #Astrocyte lineage #28 de genes


write.csv(de_bts2_all_enricher_out_1, file = 'de_bts2_all_enricher_out_1.csv')
write.csv(de_bts2_all_enricher_out_2, file = 'de_bts2_all_enricher_out_2.csv')
write.csv(de_bts2_all_enricher_out_3, file = 'de_bts2_all_enricher_out_3.csv')
write.csv(de_bts2_all_enricher_out_4, file = 'de_bts2_all_enricher_out_4.csv')
write.csv(de_bts2_all_enricher_out_5, file = 'de_bts2_all_enricher_out_5.csv')
write.csv(de_bts2_all_enricher_out_6, file = 'de_bts2_all_enricher_out_6.csv')



dim(de_bts2_all_tb_pos_entrez[[1]])
dim(de_bts2_all_tb_pos_entrez[[2]])
dim(de_bts2_all_tb_pos_entrez[[3]])
dim(de_bts2_all_tb_pos_entrez[[4]])
dim(de_bts2_all_tb_pos_entrez[[5]])
dim(de_bts2_all_tb_pos_entrez[[6]])
tibble(cluster1='microglia', cluster4='olidendocyte', others='astrocyte lineage')

#marker for de 2356
de_bts2_all_enricher_out
de_bts5_2_tb_entrez <- inner_join(de_bts5_2_tb, annotation_ahb, by = c('gene'='gene_name'))
de_bts2_6_tb_entrez <- inner_join(de_bts2_6_tb, annotation_ahb, by = c('gene'='gene_name'))
de_bts6_3_tb_entrez <- inner_join(de_bts6_3_tb, annotation_ahb, by = c('gene'='gene_name'))
de_bts5_6_tb_entrez <- inner_join(de_bts5_6_tb, annotation_ahb, by = c('gene'='gene_name'))

c('BCAN','CHI3L1','SERPINE1','C6orf15')

de_bts5_2_enricher_out <- get_cell_marker(de_bts5_2_tb_entrez) %>% data.frame() 
de_bts2_6_enricher_out <- get_cell_marker(de_bts2_6_tb_entrez) %>% data.frame()
de_bts6_3_enricher_out <- get_cell_marker(de_bts6_3_tb_entrez) %>% data.frame()
de_bts5_6_enricher_out <- get_cell_marker(de_bts5_6_tb_entrez) %>% data.frame()



#plot vague point gene change#####################################
test_bts2data@connect_cluster$vague_point_list
plot(test_bts2data@trajectory_plot$plot)


plot(eval_bts2data$plot)
plot(eval_bts2data_clu$plot)

#1_cname
cluster_1_cname <- test_bts2data@trajectory$result_list$`1`

#4
cluster_4_cname <- test_bts2data@trajectory$result_list$`4`

#2,3,5,6
cluster_2356_cname <- c(rownames(test_bts2data@trajectory$result_list$`2`), rownames(test_bts2data@trajectory$result_list$`3`), rownames(test_bts2data@trajectory$result_list$`5`), rownames(test_bts2data@trajectory$result_list$`6`))

#1_time
cluster_1_time <- test_bts2data@pseudotime$`1`

#4
cluster_4_time <- test_bts2data@pseudotime$`4`

#2,3,5,6
cluster_2356_time <- c(test_bts2data@pseudotime$`2`, test_bts2data@pseudotime$`3`, test_bts2data@pseudotime$`5`, test_bts2data@pseudotime$`6`)

#correlation gene with time
names(sort(cluster_2356_time))

bts2_for_tsne <- test_bts2data@normalized_data[test_bts2data@selected_feature[["selected_gene"]], ]
bts2_for_tsne_2356 <- bts2_for_tsne
colnames(bts2_for_tsne_2356) = as.character(1:length(colnames(bts2_for_tsne)))
bts2_for_tsne_2356_sort <- bts2_for_tsne_2356[, names(sort(cluster_2356_time))]


cor_2356 <- as.tibble(cor(sort(cluster_2356_time), t(bts2_for_tsne_2356_sort)))
cor_2356_name <- cor_2356 %>% sort() 


cor_2356_pos <- names(cor_2356_name[1981:2000])
cor_2356_neg <- names(cor_2356_name[1:20])

#var
var_2356 <- bts2_for_tsne_2356_sort %>% apply(1, var) %>% as.tibble(rownames = 'gene') %>% 
  arrange(desc(value)) %>% top_n(n = 500)

var_2356_name <- var_2356$gene

#enrich high cor genes

cor_2356_pos_entrez <- inner_join(data.frame(gene = cor_2356_pos), annotation_ahb, by = c('gene'='gene_name'))
cor_2356_neg_entrez <- inner_join(data.frame(gene = cor_2356_neg), annotation_ahb, by = c('gene'='gene_name'))

cor_2356_pos_enricher <- get_cell_marker(cor_2356_pos_entrez) %>% data.frame()
cor_2356_neg_enricher <- get_cell_marker(cor_2356_neg_entrez) %>% data.frame()

#enrich high var genes

var_2356_entrez <- inner_join(data.frame(gene = var_2356_name), annotation_ahb, by = c('gene'='gene_name'))

var_2356_enricher <- get_cell_marker(var_2356_entrez) %>% data.frame()

#gene plot############################
cluster_2356_vague_name <- c(test_bts2data@connect_cluster$vague_point_list$`2,5`,
  test_bts2data@connect_cluster$vague_point_list$`2,6`,
  test_bts2data@connect_cluster$vague_point_list$`3,6`) %>% unique()

cluster_2356_cell_name_index_2356 <- test_bts2data@distribution_estimation$cluster_index %in% c(2,3,5,6)
cluster_2356_cell_name <- names(test_bts2data@distribution_estimation$cluster_index[cluster_2356_cell_name_index_2356])

#gene data
line_plot_gene_data <- test_bts2data@normalized_data
colnames(line_plot_gene_data) <- as.character(1:length(colnames(line_plot_gene_data)))

#sorted time
line_plot_gene_data <- line_plot_gene_data[, as.numeric(names(sort(cluster_2356_time)))]

line_plot_gene_data <- t(line_plot_gene_data[test_bts2data@selected_feature$selected_gene, ]) %>% as.tibble(rownames = 'cell') %>% 
  add_column(time = sort(cluster_2356_time), .after = 'cell') %>% add_column(index = 1:length(sort(cluster_2356_time)), .after = 'cell')


vindex_2356 <- line_plot_gene_data_pos[line_plot_gene_data_pos$cell %in% cluster_2356_vague_name, ]$index
var_2356[var_2356$gene%in% c('SLC2A1' ,'NOTCH2' ,'VIM'    ,'METRN'  ,'GFAP'  ,'ID2' ,   'CTNNB1' ,'FABP7' , 'PMP2',   'NOTCH1'), ]

vindex_2_5 <- line_plot_gene_data_pos[line_plot_gene_data_pos$cell %in% test_bts2data@connect_cluster$vague_point_list$`2,5`, ]$index
vindex_2_6 <- line_plot_gene_data_pos[line_plot_gene_data_pos$cell %in% test_bts2data@connect_cluster$vague_point_list$`2,6`, ]$index
vindex_3_6 <- line_plot_gene_data_pos[line_plot_gene_data_pos$cell %in% test_bts2data@connect_cluster$vague_point_list$`3,6`, ]$index

entrez2gene(c(2670,7431,6513,3398,2173,4853,79006,4851,5375,1499))

g <- ggplot() + geom_point(data = line_plot_gene_data, aes(index, GFAP)) + 
  geom_smooth(data = line_plot_gene_data, aes(index, GFAP)) + 
  geom_smooth(data = line_plot_gene_data, aes(index, BCAN)) + 
  geom_vline(xintercept = vindex_2356, linetype="dotted", 
             color = "blue", size=1.5)

plot(g)  
plot(eval_bts2data_clu$plot)

c('BCAN','CHI3L1','SERPINE1','C6orf15')



#plot
pdf(file = paste0(path, 'gene_line_plot.pdf'),width = 10,height = 7)
gene_line_plot(gene_name = c('BCAN','CHI3L1','SERPINE1','C6orf15','GFAP'), 
               vague_cell_cluster = c('2_5', '3_6'),
               limit = TRUE, xlim = c(0,500), ylim = c(0, 200))

dev.off()


pdf(file = paste0(path, 'bts2data_trajectory_plot.pdf'),width = 10,height = 10)
plot(test_bts2data@trajectory_plot$plot)
dev.off()


pdf(file = paste0(path, 'bts2data_trajectory_clu_plot.pdf'),width = 10,height = 10)
plot(eval_bts2data_clu$plot)
dev.off()












annotation_ahb[annotation_ahb$entrezid %in% c(183,11341, 58473, 57447, 11170, 9415, 2670), ]$gene_name

com <- lapply(de_bts2_all_tb_pos_entrez, function(x) dim(x)[1])
com4<-sum(com[[4]])/sum(unlist(com))
com1<-sum(com[[1]])/sum(unlist(com))
como<-1-sum(com1,com4)
############################function
connect_cluster_vague_name<-function(object, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01){
  
  counted_possibility = object@point_possibility$counted_possibility
  n=ncol(counted_possibility)
  vague_ratio<-list()
  cluster_connection<-matrix(0, ncol = n, nrow = n)
  
  counted_possibility<-counted_possibility[!is.na(counted_possibility)[,1], ]
  vague_point_list = list()
  #i=first distribution for comparison and j= the second
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      
      ratio<-sum((counted_possibility[,i]+counted_possibility[,j])>=sum_cri & 
                   abs(counted_possibility[,i]-counted_possibility[,j])<=diff_cri)/
        sum((counted_possibility[,i]+counted_possibility[,j])>=sum_cri)
      
      if (is.na(ratio)) ratio = 0
      
      list_name<-paste('ratio',i,j, sep = '_')
      vague_ratio[[list_name]]<-ratio
      
      #determine cluster_trajectory matrix
      if (ratio>=vague_cri) cluster_connection[i,j]<-cluster_connection[j,i]<-1
      
      ###############################vague_point_vector
      vague_point_list_index <- (counted_possibility[,i]+counted_possibility[,j])>=sum_cri & 
        abs(counted_possibility[,i]-counted_possibility[,j])<=diff_cri
      vague_point_name <- rownames(counted_possibility)[vague_point_list_index]
      
      
      vague_point_list[[paste(i, j, sep = ',')]] <-  vague_point_name
      
      
      
    }
  }
  rownames(cluster_connection)<-colnames(cluster_connection)<-1:nrow(cluster_connection)
  result<-list(cluster_connection=cluster_connection, vague_ratio=vague_ratio, vague_point_list = vague_point_list)
  object@connect_cluster <- result
  return(object)
}


plot_trajectory_vague_name<-function(object, plot_title = '', width = 0.2, height = 0.2){
  
  pseudotime = object@pseudotime
  trajectory = object@trajectory
  connection_matrix = object@connect_cluster$cluster_connection
  vague_point_list = object@connect_cluster$vague_point_list
  
  
  g<-igraph::graph_from_adjacency_matrix(connection_matrix, mode = 'undirected')
  g_comp<-igraph::components(g)
  n <- tail(sort(table(g_comp$membership)), 1) + 1
  ordered_pseudotime<-pseudotime[[1]]
  ordered_trajectory<-trajectory$result_list[[1]]
  for (i in pseudotime[-1]) ordered_pseudotime<-rbind(as.matrix(i), as.matrix(ordered_pseudotime))
  for (i in trajectory$result_list[-1]) ordered_trajectory<-rbind(as.matrix(i), as.matrix(ordered_trajectory))
  
  ordered_pseudotime<-ordered_pseudotime[order(as.numeric(rownames(ordered_pseudotime))),]
  ordered_trajectory<-ordered_trajectory[order(as.numeric(rownames(ordered_trajectory))),]
  
  gg_data<-as.data.frame(ordered_trajectory)
  
  #vague point data
  x = unlist(test_bts2data@connect_cluster$vague_point_list)
  vague_point_gg_data <- gg_data[x[!duplicated(x)],]
  ###################
  g<-ggplot()
  
  for (i in trajectory[["lines_list"]]) {
    line <- as.data.frame(i$s[i$ord,])
    colnames(line) <- c('x', 'y')
    g <- g + geom_path(data = line, aes(x,y), size = 1.5)
  }
  
  #vague point data
  g <- g + geom_point(data = gg_data[-as.numeric(rownames(vague_point_gg_data)), ], aes(x,y, fill = ordered_pseudotime[-as.numeric(rownames(vague_point_gg_data))]), color = 'black', size=8, pch =21, stroke = 0.5, position = position_jitter(width = width, height = height)) +
    scale_fill_gradientn(name = 'Time', colors = colorRampPalette(c("#1b98e0", "red"))(5)) + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                         axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                                         axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  g <- g + geom_point(data = vague_point_gg_data, aes(x,y), fill = 'grey', color = 'black', size=8, pch =21, stroke = 0.5, position = position_jitter(width = width + 0.2, height = height + 0.2))
  #######################
  if (length(plot_title) >= 1) g <- g + ggtitle(label = plot_title) 
  result = list(ordered_pseudotime = ordered_pseudotime, ordered_trajectory = ordered_trajectory, plot = g)
  object@trajectory_plot<- result
  return(object)
}


entrez2gene <- function(entrez_id){
  out <- annotation_ahb[annotation_ahb$entrezid %in% entrez_id, ]$gene_name
  out
}


gene_line_plot <- function(gene_name, vague_cell_cluster, limit = FALSE, xlim, ylim){
  
  line_data <- line_plot_gene_data[, c( 'index', gene_name)] %>%
    pivot_longer(cols = -1, names_to = 'gene', values_to = 'value') 
  
  
  vline_data <- tibble(x = c(test_bts2data@connect_cluster$vague_point_list$`2,5`, 
                             test_bts2data@connect_cluster$vague_point_list$`2,6`,
                             test_bts2data@connect_cluster$vague_point_list$`3,6`),
                       index = c(rep('2_5', length(test_bts2data@connect_cluster$vague_point_list$`2,5`)),
                                 rep('2_6', length(test_bts2data@connect_cluster$vague_point_list$`2,6`)),
                                 rep('3_6', length(test_bts2data@connect_cluster$vague_point_list$`3,6`)))) %>%
    filter(index %in% vague_cell_cluster)
  
  
  
  
  vline_data <- inner_join(vline_data, line_plot_gene_data_pos, by = c('x'='cell')) %>% 
    select(x, index.x, index.y) %>%
    rename(index = index.x, vx = index.y)
  
  
  g <- ggplot()
  g <- g + 
    geom_vline(data = vline_data, aes(color = index, xintercept = vx), size = 1.5)
  
  g <- g + 
    geom_smooth(data = line_data, aes(x = index, y = value, color = gene), size = 2.5)
  g <- g + theme_classic()+xlab('cell index')+ylab('gene') + 
    theme(legend.title=element_text(size=18, face="bold"), legend.text=element_text(size=16),
          axis.text.y = element_text(size=18), axis.text.x = element_text(size=18),
          axis.title.y= element_text(size=18), axis.title.x= element_text(size=18))
  
  
  if (limit==TRUE) g <- g + coord_cartesian(xlim=xlim, ylim=ylim)
  plot(g)
}


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
