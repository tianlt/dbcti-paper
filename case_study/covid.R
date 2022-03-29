library(tidyverse)
library(GEOquery)
library(clusterProfiler)
library(dbcti)
path <- ''
path_case = ''

#annotation data
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>% 
  tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID = strsplit(geneID, ', '))
annotation_ahb <- read.csv('datasets/annotations_ahb.txt')


x=getGEO(filename=paste0(path, 'GSE171524_series_matrix.txt'), getGPL = FALSE) %>% data.frame() %>%
  rownames_to_column(var = 'name') %>% as.tibble() %>% select(name, title) %>% 
  mutate(index = substring(title, 1, 1)) %>% group_by(index)

ct <- x$name[1:7]
cv <- x$name[8:27]

L01cov = read.csv('GSM5226581_L01cov_raw_counts.csv.gz', row.names = 1)
dim(L01cov)

#y##########################
set.seed(122)
L01cov <- as.data.frame(L01cov)
test_L01cov <- create_object(L01cov, normalized = FALSE) 
test_L01cov <- normalize(test_L01cov, gene_cri = 1, cell_cri = 1)
test_L01cov <- select_var_feature(test_L01cov, use_normalized_data = TRUE, n = 2000)
test_L01cov <- tsneplot(test_L01cov, use_normalized_data = TRUE, specified_gene = FALSE, pca = TRUE, perplexity = 10)

test_L01cov <- contour_plot(test_L01cov)
test_L01cov <- distribution_estimation(test_L01cov, ndraw = 1000, expansion = 2.5, ... = c(1,4),2,3,6,5,7,9,8) 
test_L01cov <- point_possibility(test_L01cov, r = 5)
test_L01cov <- connect_cluster(test_L01cov, sum_cri = 0.8, diff_cri = 0.3, vague_cri = 0.05)
test_L01cov <- infer_trajectory(test_L01cov, iter_n =50)
test_L01cov <- calculate_pseudotime(test_L01cov, start_state_name = c('2', '3' ,'8'))
test_L01cov <- plot_trajectory(test_L01cov)
plot(test_L01cov@trajectory_plot$plot)

#evaluation

eval_L01cov_clu <- plot_index(index = factor(t(test_L01cov@distribution_estimation[["cluster_index"]])), trajectory = test_L01cov@trajectory, connection_matrix = test_L01cov@connect_cluster$cluster_connection)
plot(eval_L01cov_clu$plot)

pdf(file = paste0(path_case, 'L01cov_trajectory_clu_plot.pdf'),width = 10,height = 10)
plot(eval_L01cov_clu$plot)
dev.off()


#de
library(Seurat)
L01cov_seurat <- CreateSeuratObject(counts = L01cov)
L01cov_seurat <- NormalizeData(L01cov_seurat, normalization.method = "LogNormalize")
L01cov_seurat <- ScaleData(L01cov_seurat)
L01cov_seurat <- FindVariableFeatures(L01cov_seurat, selection.method = "vst", nfeatures = 2000)
L01cov_seurat <- RunPCA(L01cov_seurat, do.print=FALSE)
Idents(L01cov_seurat) <- factor((test_L01cov@distribution_estimation[["cluster_index"]]))

de_L01cov_all<-FindAllMarkers(L01cov_seurat, logfc.threshold = 0.1)
de_L01cov_1_2<-FindMarkers(L01cov_seurat, ident.1 = c(1,2),logfc.threshold = 0.1, only.pos = TRUE)
de_L01cov_8<-FindMarkers(L01cov_seurat, ident.1 = c(8),logfc.threshold = 0.1, only.pos = TRUE)
de_L01cov_3_7<-FindMarkers(L01cov_seurat, ident.1 = c(3,4,5,6,7),logfc.threshold = 0.1, only.pos = TRUE)
de_L01cov_3_4<-FindMarkers(L01cov_seurat, ident.1 = 3, ident.2 = 4, logfc.threshold = 0.1)
de_L01cov_2_3<-FindMarkers(L01cov_seurat, ident.1 = 2, ident.2 = 3, logfc.threshold = 0.1)

de_L01cov_all_tb_pos <- de_L01cov_all %>% 
  as_tibble() %>% group_by(cluster) %>% 
  filter(avg_logFC > 0.2 & p_val_adj < 0.05) 

#de to cell marker
de_L01cov_all_tb_pos_entrez <- inner_join(de_L01cov_all_tb_pos, annotation_ahb, by = c('gene'='gene_name')) %>%
  group_split()

get_cell_marker <- function(x) {
  out <- enricher(x$entrezid, TERM2GENE=cell_markers, minGSSize=1)
  out
}

de_L01cov_all_enricher_out <- lapply(de_L01cov_all_tb_pos_entrez, get_cell_marker)

#Immune cell lineage
#REF The Cellular Origin of Activated Fibroblasts in the Infarcted and Remodeling Myocardium
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[5]])) #Exhausted CD4+ T cell/Immune cell
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[2]])) #monocyte-dereived/MYELOID CELL/Macrophage/B cell/Exhausted CD4+ T cell/IMMUNE

#FIBEROBLAST/AT/EPILETHAIL CELL LINEAGE
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[3]])) #interstitial-derived/Exhausted CD8+ T cell/Macrophage/Specialist antigen presenting cell
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[4]])) #MKI67+ progenitor cell/Endothelial cell/Hemangioblast/at1/Mesenchymal/basal stem cell/progeniter
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[1]])) #fiberblast/myofiberblast
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[6]])) #at1 cell/Ionocyte cell
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[7]])) #epithelial cell/varianet clare cell

#Ciliated epithelial cell lineage
DT::datatable(as.data.frame(de_L01cov_all_enricher_out[[8]])) #Ciliated epithelial cell

find_marker <- function(x, marker) {
  out <- sum(marker %in% x$gene)
  out
}

#cell marker annotation
#from Combinations of differentiation markers distinguish subpopulations of alveolar epithelial cells in adult lung
#The use of alveolar epithelial type I cell-selective markers to investigate lung injury and repair
at1<-c('AQP5','HOPX','PDPN','CAV1','CAV2','CYP2B1')
at2<-c('NKX2.1','SFTPC')
#CD248 and integrin alpha-8 are candidate markers for differentiating lung fibroblast subtypes
#Collagen-producing lung cell atlas identifies multiple subsets with distinct localization and relevance to fibrosis
fib<-c('CD248','ITGA8','Pdgfra', 'Tcf21', 'Npnt', 'Col13a1','Cthrc1')#Cthrc1 pathology
fib<-toupper(fib)

#Identification of a Zeb1 expressing basal stem cell subpopulation in the prostate

bas<-c('TP63','KRT5','NGFR','ITGA6', 'CD117', 'CD133', 'CD44', 'Trop2', 'CD49f', 'Sca1')
bas<-toupper(bas)

#Novel Method for Isolation of Murine Clara Cell Secretory Protein-Expressing Cells with Traces of Stemness
cla<-c('CCSP','SCGB1A1')

#Monocyte-derived Alveolar Macrophages: The Dark Side of Lung Repair?
moam<-c('SIGLEC5','CD170','FCGR1A','CD64','ITGAX', 'CD11C','ADGRE1', 'MERTK')#'SIGLEC5','CD170' negative
tram
inm<-c('FCGR1A', 'CD64', 'ITGAM', 'CD11B', 'CD11C', 'ITGAX')

#Fibroblasts and Myofibroblasts: What are we talking about?
myofib<-c('VIM','ACTA2', 'PALLD', 'PXN','VCL','TNS1')
lapply(de_L01cov_all_tb_pos_entrez, find_marker, at1)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, at2)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, fib)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, bas)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, myofib)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, cla)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, moam)
lapply(de_L01cov_all_tb_pos_entrez, find_marker, inm)
bas[bas %in% de_L01cov_all_tb_pos_entrez[[4]]$gene]
bas[bas %in% de_L01cov_all_tb_pos_entrez[[3]]$gene]
fib[fib %in% de_L01cov_all_tb_pos_entrez[[1]]$gene]
moam[moam %in% de_L01cov_all_tb_pos_entrez[[2]]$gene]
moam[moam %in% de_L01cov_all_tb_pos_entrez[[3]]$gene]
inm[inm %in% de_L01cov_all_tb_pos_entrez[[2]]$gene]
inm[inm %in% de_L01cov_all_tb_pos_entrez[[3]]$gene]
at1[at1 %in% de_L01cov_all_tb_pos_entrez[[6]]$gene]

#composition
com <- lapply(de_L01cov_all_tb_pos_entrez, function(x) dim(x)[1])
com2_5<-sum(com[[2]],com[[5]])/sum(unlist(com))
com8<-sum(com[[8]])/sum(unlist(com))
com2_5<-1-sum(com2_5,com8)


#3 and 4 pathway
library(org.Hs.eg.db)
de_L01cov_3_4_tibble <- de_L01cov_3_4 %>% mutate(index = case_when(avg_logFC >= 0 ~ 1,
                                           avg_logFC < 0 ~ 0)) %>%
  filter(p_val_adj <= 0.05) %>% rownames_to_column('gene') %>% as.tibble() %>%
  group_by(index) %>% mutate(abs_log = abs(avg_logFC)) %>%
  top_n(n = 100, wt = abs_log) %>%   
  arrange(desc(abs_log), .by_group = TRUE) %>% group_split()

L01cov_bk_gene <- rownames(L01cov)
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = L01cov_bk_gene,  # data to use for retrieval
                                           columns = c("ENSEMBL", "ENTREZID","GENENAME"), # information to retreive for given data
                                           keytype = "SYMBOL") # type of data given in 'keys' argument



annotations_orgDb <- annotations_orgDb[which(duplicated(annotations_orgDb$SYMBOL) == FALSE), ]

# annotation_ahb <- read.csv('C:/Users/13598533/Documents/Functional_analysis/Functional_analysis/data/annotations_ahb.txt')

de_L01cov_3_4_ahb <- lapply(de_L01cov_3_4_tibble, inner_join, y = annotation_ahb, by = c('gene'='gene_name'))

L01cov_bk_gene_ahb <- tibble(gene = L01cov_bk_gene) %>%
  inner_join(y = annotation_ahb, by = c('gene'='gene_name'))
  

de_L01cov_3_4_genes <- lapply(de_L01cov_3_4_ahb, function(x) as.character(x$gene_id))
all_L01cov_genes <-  as.character(L01cov_bk_gene_ahb$gene_id)

enrichGO_map_3_4 <- function(x){
  ego<-enrichGO(gene = de_L01cov_3_4_genes[[x]],
                universe = all_L01cov_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'ALL',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrichGO_map_3_4_mf <- function(x){
  ego<-enrichGO(gene = de_L01cov_3_4_genes[[x]],
                universe = all_L01cov_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'MF',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrichGO_map_3_4_cc <- function(x){
  ego<-enrichGO(gene = de_L01cov_3_4_genes[[x]],
                universe = all_L01cov_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'CC',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrich_L01cov_3_4 <- map(c(1:2), enrichGO_map_3_4)

summary_enrich_L01cov_3_4 <- lapply(enrich_L01cov_3_4, data.frame)
pdf(file = paste0(path_case, 'cluster4_3de_enrich.pdf'),width = 11,height = 6)

dotplot(enrich_L01cov_3_4[[1]], showCategory=20) #3<4 stem de
dev.off()
pdf(file = paste0(path_case, 'cluster3_4de_enrich.pdf'),width = 11,height = 6)
dotplot(enrich_L01cov_3_4[[2]], showCategory=20) #3>4 interstitial macrophages de
dev.off()
#mf
enrich_L01cov_3_4_mf <- map(c(1:2), enrichGO_map_3_4_mf)

summary_enrich_L01cov_3_4_mf <- lapply(enrich_L01cov_3_4_mf, data.frame)
pdf(file = paste0(path_case, 'cluster4_3de_enrich_mf.pdf'),width = 11,height = 6)

dotplot(enrich_L01cov_3_4_mf[[1]], showCategory=20) #3<4 stem de
dev.off()
pdf(file = paste0(path_case, 'cluster3_4de_enrich_mf.pdf'),width = 11,height = 6)

dotplot(enrich_L01cov_3_4_mf[[2]], showCategory=20) #3>4 interstitial macrophages de
dev.off()
#cc
enrich_L01cov_3_4_cc <- map(c(1:2), enrichGO_map_3_4_cc)

summary_enrich_L01cov_3_4_cc <- lapply(enrich_L01cov_3_4_cc, data.frame)
pdf(file = paste0(path_case, 'cluster4_3de_enrich_cc.pdf'),width = 11,height = 6)

dotplot(enrich_L01cov_3_4_cc[[1]], showCategory=20) #3<4 stem de
dev.off()

pdf(file = paste0(path_case, 'cluster3_4de_enrich_cc.pdf'),width = 11,height = 6)

dotplot(enrich_L01cov_3_4_cc[[2]], showCategory=20) #3>4 interstitial macrophages de
dev.off()
#2 and 3 pathway
de_L01cov_2_3_tibble <- de_L01cov_2_3 %>% mutate(index = case_when(avg_logFC >= 0 ~ 1,
                                                                   avg_logFC < 0 ~ 0)) %>%
  filter(p_val_adj <= 0.05) %>% rownames_to_column('gene') %>% as.tibble() %>%
  group_by(index) %>% mutate(abs_log = abs(avg_logFC)) %>%
  top_n(n = 100, wt = abs_log) %>%   
  arrange(desc(abs_log), .by_group = TRUE) %>% group_split()


de_L01cov_2_3_ahb <- lapply(de_L01cov_2_3_tibble, inner_join, y = annotation_ahb, by = c('gene'='gene_name'))

L01cov_bk_gene_ahb <- tibble(gene = L01cov_bk_gene) %>%
  inner_join(y = annotation_ahb, by = c('gene'='gene_name'))


de_L01cov_2_3_genes <- lapply(de_L01cov_2_3_ahb, function(x) as.character(x$gene_id))
all_L01cov_genes <-  as.character(L01cov_bk_gene_ahb$gene_id)

enrichGO_map_2_3 <- function(x){
  ego<-enrichGO(gene = de_L01cov_2_3_genes[[x]],
                universe = all_L01cov_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'ALL',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrichGO_map_2_3_cc <- function(x){
  ego<-enrichGO(gene = de_L01cov_2_3_genes[[x]],
                universe = all_L01cov_genes,
                keyType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'CC',
                qvalueCutoff = 0.05,
                readable = TRUE)
  ego
}

enrich_L01cov_2_3 <- map(c(1:2), enrichGO_map_2_3)

summary_enrich_L01cov_2_3 <- lapply(enrich_L01cov_2_3, data.frame)
#plot
pdf(file = paste0(path_case, 'cluster3_2de_enrich.pdf'),width = 11,height = 6)
dotplot(enrich_L01cov_2_3[[1]], showCategory=20) #2<3 interstitial macrophages de
dev.off()

#plot
pdf(file = paste0(path_case, 'cluster2_3de_enrich.pdf'),width = 11,height = 6)
dotplot(enrich_L01cov_2_3[[2]], showCategory=20) #2>3 monocyte macrophages de
dev.off()


enrich_L01cov_2_3_cc <- map(c(1:2), enrichGO_map_2_3_cc)

summary_enrich_L01cov_2_3_cc <- lapply(enrich_L01cov_2_3_cc, data.frame)

#plot
pdf(file = paste0(path_case, 'cluster3_2de_enrich_cc.pdf'),width = 11,height = 6)
dotplot(enrich_L01cov_2_3_cc[[1]], showCategory=20) #2<3 interstitial macrophages de
dev.off()
#plot
pdf(file = paste0(path_case, 'cluster2_3de_enrich_cc.pdf'),width = 11,height = 6)
dotplot(enrich_L01cov_2_3_cc[[2]], showCategory=20) #2>3 monocyte macrophages de
dev.off()

#heatmap 3-4
data_3_4_gene <- de_L01cov_3_4 %>% mutate(index = case_when(avg_logFC >= 0 ~ 1,
                                                                   avg_logFC < 0 ~ 0)) %>%
  filter(p_val_adj <= 0.05) %>% rownames_to_column('gene') %>% as.tibble() %>%
  group_by(index) %>% mutate(abs_log = abs(avg_logFC)) %>%
  top_n(n = 25, wt = abs_log) %>%   
  arrange(desc(abs_log), .by_group = TRUE)  %>% ungroup() %>%
  pull(gene) %>% unique()
 

data_3_cell <- test_L01cov@distribution_estimation[["cluster_index"]] == 3 
data_4_cell <- test_L01cov@distribution_estimation[["cluster_index"]] == 4
# data_3_4 <- L01cov[data_3_4_gene, data_3_cell] %>% 
#   cbind(L01cov[data_3_4_gene, data_4_cell]) %>%
#   #apply(1, function(x) (x-mean(x))/sd(x)) %>%
#   #t() %>%
#   apply(1, scales::rescale, to = c(-1,1)) %>%
#   t() %>% data.frame() %>%
#   rename_all(~as.character(1:sum(data_3_4_cell))) %>%
#   rownames_to_column('gene') %>% as_tibble() %>%
#   pivot_longer(cols=-1, names_to = "cell", values_to = "value") 
  
data_3_4 <- L01cov[data_3_4_gene, data_3_cell] %>% 
  cbind(L01cov[data_3_4_gene, data_4_cell]) %>% 
  log1p() %>%
  apply(1, function(x) (x-mean(x))/sd(x)) %>% 
  t() %>% 
  data.frame() %>% 
  rename_all(~as.character(1:sum(data_3_4_cell))) %>%
  rownames_to_column('gene') %>% as_tibble() %>%
  pivot_longer(cols=-1, names_to = "cell", values_to = "value") 

data_3_4$gene <- factor(data_3_4$gene, levels = unique(data_3_4$gene))
data_3_4$cell <- factor(data_3_4$cell, levels = unique(data_3_4$cell))
#plot

pdf(file = paste0(path_case, 'data_3_4_heatmap.pdf'),width = 7,height = 7)
ggplot(data_3_4,aes(x=cell,y=gene,fill=value))+
  geom_tile() +  
  scale_fill_gradientn(colours = c("black", "cyan", "red")) +
  theme(axis.title.x=element_text(size = 14, face = "bold"),
        axis.title.y=element_text(size = 14, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=14, face="bold"),
        legend.text=element_text(size=14),
        axis.text.y=element_text(face = "bold")) 
  
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
