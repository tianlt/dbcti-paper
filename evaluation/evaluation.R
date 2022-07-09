path = 'C:/Users/13598533/OneDrive - UTS/Desktop/cell_cycle_trajectory/DBCTI/dbcti_package/code_upload/dbcti-paper/simulation_and_real_data/new_code'
setwd(path)
library(dbcti)
library(Rtsne)
library(phateR)
library(dca)
library(ggplot2)
library(tidyverse)
####################
test_leng<-readRDS('test_hescmt.rds')
test_nestorowa<-readRDS('test_nestorowa_batch1_type_filter_no_stem.rds')
test_kowalczyk <- readRDS('test_stem_mouse_C57BL6_old.rds')
test_camp<-readRDS('test_camp1_batch1.rds')
test_yan<-readRDS('test_yan.rds')

# leng_data <- readRDS('hescmt_for_new.rds')
leng_index <- readRDS('hescsmsinfo_cycle_index_for_new.rds')

# nestorowa_data <- readRDS('nestorowa_data_batch1_type_filter_no_stem_for_new.rds')
nestorowa_index <- readRDS('nestorowa_index_batch1_type_filter_no_stem_for_new.rds')
kowalczyk_index <- readRDS('stem_mouse_C57BL6_index_old_for_new.rds')
camp_index <- readRDS('camp1_index_batch1_for_new.rds')
yan_index <- readRDS('yan_index_for_new.rds')


# color by reference
plot(test_leng@tsne_data$Y, col = leng_index)
plot(test_nestorowa@tsne_data$Y, col = as.factor(nestorowa_index))
plot(test_kowalczyk@tsne_data$Y, col = as.factor(kowalczyk_index))
plot(test_camp@tsne_data$Y, col = as.factor(camp_index))
plot(test_yan@tsne_data$Y, col = as.factor(yan_index))


saveRDS(test_leng@tsne_data$Y,'leng_data_dbcti_tsne.rds')
saveRDS(test_nestorowa@tsne_data$Y,'nestorowa_data_dbcti_tsne.rds')
saveRDS(test_kowalczyk@tsne_data$Y,'kowalczyk_data_dbcti_tsne.rds')
saveRDS(test_camp@tsne_data$Y,'camp_data_dbcti_tsne.rds')
saveRDS(test_yan@tsne_data$Y,'yan_data_dbcti_tsne.rds')
# color by dbcti cluster
plot(test_leng@tsne_data$Y, col = test_leng@distribution_estimation[["cluster_index"]])
plot(test_nestorowa@tsne_data$Y, col = test_nestorowa@distribution_estimation[["cluster_index"]])
plot(test_kowalczyk@tsne_data$Y, col = test_kowalczyk@distribution_estimation[["cluster_index"]])
plot(test_camp@tsne_data$Y, col = test_camp@distribution_estimation[["cluster_index"]])
plot(test_yan@tsne_data$Y, col = test_yan@distribution_estimation[["cluster_index"]])

leng_data_gene = test_leng@selected_feature[["selected_gene"]]
leng_data_to_dr = test_leng@normalized_data[leng_data_gene,]
nestorowa_data_gene = test_nestorowa@selected_feature$selected_gene
nestorowa_data_to_dr = test_nestorowa@normalized_data[nestorowa_data_gene,]
kowalczyk_data_gene = test_kowalczyk@specified_gene
kowalczyk_data_to_dr = test_kowalczyk@normalized_data[kowalczyk_data_gene,]
camp_data_gene = test_camp@selected_feature$selected_gene
camp_data_to_dr = test_camp@normalized_data[camp_data_gene,]
yan_data_gene = test_yan@selected_feature$selected_gene
yan_data_to_dr = test_yan@normalized_data[yan_data_gene,]


saveRDS(leng_data_to_dr,'leng_data_to_dr.rds')
saveRDS(nestorowa_data_to_dr,'nestorowa_data_to_dr.rds')
saveRDS(kowalczyk_data_to_dr,'kowalczyk_data_to_dr.rds')
saveRDS(camp_data_to_dr,'camp_data_to_dr.rds')
saveRDS(yan_data_to_dr,'yan_data_to_dr.rds')


# from ihpc ####################################################################
library(phateR)
library(DCA)
library(NbClust)
library(NMF)

leng_data_to_dr <- readRDS('leng_data_to_dr.rds')
nestorowa_data_to_dr <- readRDS('nestorowa_data_to_dr.rds')
kowalczyk_data_to_dr <- readRDS('kowalczyk_data_to_dr.rds')
camp_data_to_dr <- readRDS('camp_data_to_dr.rds')
yan_data_to_dr <- readRDS('yan_data_to_dr.rds')

leng_data_dbcti_tsne <- readRDS('leng_data_dbcti_tsne.rds')
nestorowa_data_dbcti_tsne <- readRDS('nestorowa_data_dbcti_tsne.rds')
kowalczyk_data_dbcti_tsne <- readRDS('kowalczyk_data_dbcti_tsne.rds')
camp_data_dbcti_tsne <- readRDS('camp_data_dbcti_tsne.rds')
yan_data_dbcti_tsne <- readRDS('yan_data_dbcti_tsne.rds')


leng_index <- readRDS('hescsmsinfo_cycle_index_for_new.rds')
nestorowa_index <- readRDS('nestorowa_index_batch1_type_filter_no_stem_for_new.rds')
nestorowa_integrated_index <- nestorowa_index
for (i in 1:length(nestorowa_integrated_index)) {
  if (nestorowa_integrated_index[i] == 'MPP' | nestorowa_integrated_index[i] == 'LMPP'){
    nestorowa_integrated_index[i] <- 'MPP/LMPP'
  }
}
kowalczyk_index <- readRDS('stem_mouse_C57BL6_index_old_for_new.rds')
camp_index <- readRDS('camp1_index_batch1_for_new.rds')
camp_index <- unlist(c(camp_index))
yan_index <- readRDS('yan_index_for_new.rds')
yan_index <- unlist(c(yan_index))



# phate##################
leng_phate <- phate(t(leng_data_to_dr))
nestorowa_phate <- phate(t(nestorowa_data_to_dr))
kowalczyk_phate <- phate(t(kowalczyk_data_to_dr))
camp_phate <- phate(t(camp_data_to_dr))
yan_phate <- phate(t(yan_data_to_dr))

leng_phate_data <- leng_phate$embedding
nestorowa_phate_data <- nestorowa_phate$embedding
kowalczyk_phate_data <- kowalczyk_phate$embedding
camp_phate_data <- camp_phate$embedding
yan_phate_data <- yan_phate$embedding


plot(leng_phate$embedding, col = leng_index)
plot(nestorowa_phate$embedding, col = as.factor(nestorowa_index))
plot(kowalczyk_phate$embedding, col = as.factor(kowalczyk_index))
plot(camp_phate$embedding, col = as.factor(camp_index))
plot(yan_phate$embedding, col = as.factor(yan_index))

saveRDS(leng_phate, 'leng_phate.rds')
saveRDS(nestorowa_phate, 'nestorowa_phate.rds')
saveRDS(kowalczyk_phate, 'kowalczyk_phate.rds')
saveRDS(camp_phate, 'camp_phate.rds')
saveRDS(yan_phate, 'yan_phate.rds')
##########################
# dca
leng_dca <- dca(as.matrix(leng_data_to_dr), n.fac = 2)
nestorowa_dca <- dca(as.matrix(nestorowa_data_to_dr), n.fac = 2)
kowalczyk_dca <- dca(as.matrix(kowalczyk_data_to_dr), n.fac = 2)
camp_dca <- dca(as.matrix(camp_data_to_dr), n.fac = 2)
yan_dca <- dca(as.matrix(yan_data_to_dr), n.fac = 2)

leng_dca_data <- leng_dca$fac
nestorowa_dca_data <- nestorowa_dca$fac
kowalczyk_dca_data <- kowalczyk_dca$fac
camp_dca_data <- camp_dca$fac
yan_dca_data <- yan_dca$fac

plot(leng_dca$fac, col = leng_index)
plot(nestorowa_dca$fac, col = as.factor(nestorowa_index))
plot(kowalczyk_dca$fac, col = as.factor(kowalczyk_index))
plot(camp_dca$fac, col = as.factor(camp_index))
plot(yan_dca$fac, col = as.factor(yan_index))

saveRDS(leng_dca, 'leng_dca.rds')
saveRDS(nestorowa_dca, 'nestorowa_dca.rds')
saveRDS(kowalczyk_dca, 'kowalczyk_dca.rds')
saveRDS(camp_dca, 'camp_dca.rds')
saveRDS(yan_dca, 'yan_dca.rds')
##########################


# calculate purity
purity_kmeans <- function(data, index){
  kmeans_res = NbClust::NbClust(data, method = 'kmeans')$Best.partition
  res <- NMF::purity(as.factor(kmeans_res), as.factor(index))
  return(res)
}

names(nestorowa_index) <- NULL
names(nestorowa_integrated_index) <- NULL
names(kowalczyk_index) <- NULL
names(camp_index) <- NULL
names(yan_index) <- NULL

leng_phate_purity <- purity_kmeans(leng_phate$embedding, leng_index)
nestorowa_phate_purity <- purity_kmeans(nestorowa_phate$embedding, nestorowa_index)
kowalczyk_phate_purity <- purity_kmeans(kowalczyk_phate$embedding, kowalczyk_index)
camp_phate_purity <- purity_kmeans(camp_phate$embedding, camp_index)
yan_phate_purity <- purity_kmeans(yan_phate$embedding, yan_index)

leng_dca_purity <- purity_kmeans(leng_dca$fac, leng_index)
nestorowa_dca_purity <- purity_kmeans(nestorowa_dca$fac, nestorowa_index)
kowalczyk_dca_purity <- purity_kmeans(kowalczyk_dca$fac, kowalczyk_index)
camp_dca_purity <- purity_kmeans(camp_dca$fac, camp_index)
yan_dca_purity <- purity_kmeans(yan_dca$fac, yan_index)

leng_dbcti_purity <- purity_kmeans(leng_data_dbcti_tsne, leng_index)
nestorowa_dbcti_purity <- purity_kmeans(nestorowa_data_dbcti_tsne, nestorowa_index)
kowalczyk_dbcti_purity <- purity_kmeans(kowalczyk_data_dbcti_tsne, kowalczyk_index)
camp_dbcti_purity <- purity_kmeans(camp_data_dbcti_tsne, camp_index)
yan_dbcti_purity <- purity_kmeans(yan_data_dbcti_tsne, yan_index)

purity_df <- data.frame(leng_phate_purity, nestorowa_phate_purity,
                        leng_dca_purity, nestorowa_dca_purity,
                        leng_dbcti_purity, nestorowa_dbcti_purity)

saveRDS(purity_df, 'purity_df.rds')

# for nestorowa_integrated_index
nestorowa_integrated_phate_purity <- purity_kmeans(nestorowa_phate$embedding, nestorowa_integrated_index)
nestorowa_integrated_dca_purity <- purity_kmeans(nestorowa_dca$fac, nestorowa_integrated_index)
nestorowa_integrated_dbcti_purity <- purity_kmeans(nestorowa_data_dbcti_tsne, nestorowa_integrated_index)

purity_integrated_df <- data.frame(leng_phate_purity, nestorowa_integrated_phate_purity,
                                   leng_dca_purity, nestorowa_integrated_dca_purity,
                                   leng_dbcti_purity, nestorowa_integrated_dbcti_purity)

saveRDS(nestorowa_integrated_index, 'nestorowa_integrated_index.rds')
saveRDS(purity_integrated_df, 'purity_integrated_df.rds')

#########################################################################
#added for kowalczyk etc.
purity_adf <- data.frame(leng_phate_purity, nestorowa_phate_purity, kowalczyk_phate_purity, camp_phate_purity, yan_phate_purity,
                         leng_dca_purity, nestorowa_dca_purity,kowalczyk_dca_purity, camp_dca_purity, yan_dca_purity,
                         leng_dbcti_purity, nestorowa_dbcti_purity, kowalczyk_dbcti_purity, camp_dbcti_purity, yan_dbcti_purity)
saveRDS(purity_adf , 'purity_adf.rds')

purity_integrated_adf <- data.frame(leng_phate_purity, nestorowa_integrated_phate_purity, kowalczyk_phate_purity, camp_phate_purity, yan_phate_purity,
                                    leng_dca_purity, nestorowa_integrated_dca_purity,kowalczyk_dca_purity, camp_dca_purity, yan_dca_purity,
                                    leng_dbcti_purity, nestorowa_integrated_dbcti_purity, kowalczyk_dbcti_purity, camp_dbcti_purity, yan_dbcti_purity)
saveRDS(purity_integrated_adf , 'purity_integrated_adf.rds')
#from ihpc########################################################################

leng_phate <- readRDS('leng_phate.rds')
nestorowa_phate <- readRDS('nestorowa_phate.rds')
kowalczyk_phate <- readRDS('kowalczyk_phate.rds')
camp_phate <- readRDS('camp_phate.rds')
yan_phate <- readRDS('yan_phate.rds')

leng_dca <- readRDS('leng_dca.rds')
nestorowa_dca <- readRDS('nestorowa_dca.rds')
kowalczyk_dca <- readRDS('kowalczyk_dca.rds')
camp_dca <- readRDS('camp_dca.rds')
yan_dca <- readRDS('yan_dca.rds')

purity_df <- readRDS('purity_df.rds')
purity_integrated_df <- readRDS('purity_integrated_df.rds')
purity_adf <- readRDS('purity_adf.rds')
purity_integrated_adf <- readRDS('purity_integrated_adf.rds')

nestorowa_integrated_index <- readRDS('nestorowa_integrated_index.rds')
#unintergrated
leng_phate_ggdf <- as_tibble(leng_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = leng_index)
nestorowa_phate_ggdf <- as_tibble(nestorowa_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = as.factor(nestorowa_index))
kowalczyk_phate_ggdf <- as_tibble(kowalczyk_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = as.factor(kowalczyk_index))
camp_phate_ggdf <- as_tibble(camp_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = as.factor(camp_index))
yan_phate_ggdf <- as_tibble(yan_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = as.factor(yan_index))

leng_dca_ggdf <- tibble(x = leng_dca$fac[,1], y = leng_dca$fac[,2]) %>%
  add_column(ref_index = leng_index)
nestorowa_dca_ggdf <- tibble(x = nestorowa_dca$fac[,1], y = nestorowa_dca$fac[,2]) %>%
  add_column(ref_index = as.factor(nestorowa_index))
kowalczyk_dca_ggdf <- tibble(x = kowalczyk_dca$fac[,1], y = kowalczyk_dca$fac[,2]) %>%
  add_column(ref_index = as.factor(kowalczyk_index))
camp_dca_ggdf <- tibble(x = camp_dca$fac[,1], y = camp_dca$fac[,2]) %>%
  add_column(ref_index = as.factor(camp_index))
yan_dca_ggdf <- tibble(x = yan_dca$fac[,1], y = yan_dca$fac[,2]) %>%
  add_column(ref_index = as.factor(yan_index))

leng_tsne_ggdf <- tibble(x = test_leng@tsne_data$Y[,1], y = test_leng@tsne_data$Y[,2]) %>%
  add_column(ref_index = leng_index)
nestorowa_tsne_ggdf <- tibble(x = test_nestorowa@tsne_data$Y[,1], y = test_nestorowa@tsne_data$Y[,2]) %>%
  add_column(ref_index = as.factor(nestorowa_index))
kowalczyk_tsne_ggdf <- tibble(x = test_kowalczyk@tsne_data$Y[,1], y = test_kowalczyk@tsne_data$Y[,2]) %>%
  add_column(ref_index = as.factor(kowalczyk_index))
camp_tsne_ggdf <- tibble(x = test_camp@tsne_data$Y[,1], y = test_camp@tsne_data$Y[,2]) %>%
  add_column(ref_index = as.factor(camp_index))
yan_tsne_ggdf <- tibble(x = test_yan@tsne_data$Y[,1], y = test_yan@tsne_data$Y[,2]) %>%
  add_column(ref_index = as.factor(yan_index))

#integrated
nestorowa_integrated_phate_ggdf <- as_tibble(nestorowa_phate$embedding) %>% rename(x = PHATE1, y = PHATE2) %>%
  add_column(ref_index = as.factor(nestorowa_integrated_index))
nestorowa_integrated_dca_ggdf <- tibble(x = nestorowa_dca$fac[,1], y = nestorowa_dca$fac[,2]) %>%
  add_column(ref_index = as.factor(nestorowa_integrated_index))
nestorowa_integrated_tsne_ggdf <- tibble(x = test_nestorowa@tsne_data$Y[,1], y = test_nestorowa@tsne_data$Y[,2]) %>%
  add_column(ref_index = as.factor(nestorowa_integrated_index))


# plot(leng_dca$fac, col = leng_index)
# plot(nestorowa_dca$fac, col = as.factor(nestorowa_index))
# plot(leng_phate$embedding, col = leng_index)
# plot(nestorowa_phate$embedding, col = as.factor(nestorowa_index))
# plot(test_leng@tsne_data$Y, col = leng_index)
# plot(test_nestorowa@tsne_data$Y, col = as.factor(nestorowa_index))



# plot
# rd plot
dr_plot <- function(data){
  g <- ggplot() + geom_point(data = data, aes(x,y, color = ref_index), size=3) + theme_classic() + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                    axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                    axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))
  return(g)
}

leng_phate_plot <- dr_plot(leng_phate_ggdf)
nestorowa_phate_plot <- dr_plot(nestorowa_phate_ggdf)
kowalczyk_phate_plot <- dr_plot(kowalczyk_phate_ggdf)
camp_phate_plot <- dr_plot(camp_phate_ggdf)
yan_phate_plot <- dr_plot(yan_phate_ggdf)

leng_dca_plot <- dr_plot(leng_dca_ggdf)
nestorowa_dca_plot <- dr_plot(nestorowa_dca_ggdf)
kowalczyk_dca_plot <- dr_plot(kowalczyk_dca_ggdf)
camp_dca_plot <- dr_plot(camp_dca_ggdf)
yan_dca_plot <- dr_plot(yan_dca_ggdf)

leng_tsne_plot <- dr_plot(leng_tsne_ggdf)
nestorowa_tsne_plot <- dr_plot(nestorowa_tsne_ggdf)
kowalczyk_tsne_plot <- dr_plot(kowalczyk_tsne_ggdf)
camp_tsne_plot <- dr_plot(camp_tsne_ggdf)
yan_tsne_plot <- dr_plot(yan_tsne_ggdf)

nestorowa_integrated_phate_plot <- dr_plot(nestorowa_integrated_phate_ggdf)
nestorowa_integrated_dca_plot <- dr_plot(nestorowa_integrated_dca_ggdf)
nestorowa_integrated_tsne_plot <- dr_plot(nestorowa_integrated_tsne_ggdf)

pdf(file = 'nar_rev_diagram/leng_phate_plot.pdf',width = 6,height = 6)
plot(leng_phate_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/nestorowa_phate_plot.pdf',width = 6,height = 6)
plot(nestorowa_phate_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/kowalczyk_phate_plot.pdf',width = 6,height = 6)
plot(kowalczyk_phate_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/camp_phate_plot.pdf',width = 6,height = 6)
plot(camp_phate_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/yan_phate_plot.pdf',width = 6,height = 6)
plot(yan_phate_plot)#
dev.off()
#
pdf(file = 'nar_rev_diagram/leng_dca_plot.pdf',width = 6,height = 6)
plot(leng_dca_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/nestorowa_dca_plot.pdf',width = 6,height = 6)
plot(nestorowa_dca_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/kowalczyk_dca_plot.pdf',width = 6,height = 6)
plot(kowalczyk_dca_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/camp_dca_plot.pdf',width = 6,height = 6)
plot(camp_dca_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/yan_dca_plot.pdf',width = 6,height = 6)
plot(yan_dca_plot)#
dev.off()
#
pdf(file = 'nar_rev_diagram/leng_tsne_plot.pdf',width = 6,height = 6)
plot(leng_tsne_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/nestorowa_tsne_plot.pdf',width = 6,height = 6)
plot(nestorowa_tsne_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/kowalczyk_tsne_plot.pdf',width = 6,height = 6)
plot(kowalczyk_tsne_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/camp_tsne_plot.pdf',width = 6,height = 6)
plot(camp_tsne_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/yan_tsne_plot.pdf',width = 6,height = 6)
plot(yan_tsne_plot)#
dev.off()
#
pdf(file = 'nar_rev_diagram/nestorowa_integrated_phate_plot.pdf',width = 6,height = 6)
plot(nestorowa_integrated_phate_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/nestorowa_integrated_dca_plot.pdf',width = 6,height = 6)
plot(nestorowa_integrated_dca_plot)#
dev.off()

pdf(file = 'nar_rev_diagram/nestorowa_integrated_tsne_plot.pdf',width = 6,height = 6)
plot(nestorowa_integrated_tsne_plot)#
dev.off()

#data
data_name <- gsub('_purity', '', colnames(purity_df))
data <- lapply(str_split(data_name,'_'), function(x) x[1]) %>% unlist() %>%
  as.factor()
dr <- as.factor(c('Phate', 'Phate', 'DCA', 'DCA', 'TSNE', 'TSNE'))

#a data
data_aname <- gsub('_purity', '', colnames(purity_adf))
adata <- lapply(str_split(data_aname,'_'), function(x) x[1]) %>% unlist() %>%
  as.factor()
adr <- as.factor(c('Phate', 'Phate', 'Phate', 'Phate', 'Phate', 
                   'DCA', 'DCA', 'DCA', 'DCA', 'DCA', 
                   'TSNE', 'TSNE', 'TSNE', 'TSNE', 'TSNE'))

#data
# unintegrated purity_integrated_ggdf
value <- c(t(purity_df[1,]))
purity_ggdf <- tibble(value = value, data_name = data_name, data = data, dr = dr)

g_purity_ggdf <- ggplot() + geom_col(data = purity_ggdf, aes(dr, value, fill = data), position = 'dodge') +
  ylab('Score') + 
  theme_classic() +
  # geom_text(data = purity_ggdf, aes(dr, value, fill = data, label = data),position = position_dodge(.9), vjust = -0.2, size = 5.3) +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), plot.title = element_text(face = 'bold', size=22), legend.title=element_text(size=22, face="bold"),legend.text=element_text(size=20))

pdf(file = 'nar_rev_diagram/g_purity_ggdf.pdf',width = 8,height = 6)
plot(g_purity_ggdf)#
dev.off()

# integrated purity_integrated_ggdf
value_integrated <- c(t(purity_integrated_df[1,]))
purity_integrated_ggdf <- tibble(value = value_integrated, data_name = data_name, data = data, dr = dr)
g_purity_integrated_ggdf <- ggplot() + geom_col(data = purity_integrated_ggdf, aes(dr, value, fill = data), position = 'dodge') +
  ylab('Score') + 
  theme_classic() +
  # geom_text(data = purity_ggdf, aes(dr, value, fill = data, label = data),position = position_dodge(.9), vjust = -0.2, size = 5.3) +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), plot.title = element_text(face = 'bold', size=22), legend.title=element_text(size=22, face="bold"),legend.text=element_text(size=20))

pdf(file = 'nar_rev_diagram/g_purity_integrated_ggdf.pdf',width = 8,height = 6)
plot(g_purity_integrated_ggdf)#
dev.off()

#a data
# unintegrated purity_integrated_ggdf
avalue <- c(t(purity_adf[1,]))
purity_aggdf <- tibble(value = avalue, data_name = data_aname, data = adata, dr = adr)

g_purity_aggdf <- ggplot() + geom_col(data = purity_aggdf, aes(dr, value, fill = data), position = 'dodge') +
  ylab('Score') + 
  theme_classic() +
  # geom_text(data = purity_ggdf, aes(dr, value, fill = data, label = data),position = position_dodge(.9), vjust = -0.2, size = 5.3) +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), plot.title = element_text(face = 'bold', size=22), legend.title=element_text(size=22, face="bold"),legend.text=element_text(size=20))

pdf(file = 'nar_rev_diagram/g_purity_aggdf.pdf',width = 8,height = 6)
plot(g_purity_aggdf)#
dev.off()

# integrated purity_integrated_ggdf
avalue_integrated <- c(t(purity_integrated_adf[1,]))
purity_integrated_aggdf <- tibble(value = avalue_integrated, data_name = data_aname, data = adata, dr = adr)
g_purity_integrated_aggdf <- ggplot() + geom_col(data = purity_integrated_aggdf, aes(dr, value, fill = data), position = 'dodge') +
  ylab('Score') + 
  theme_classic() +
  # geom_text(data = purity_ggdf, aes(dr, value, fill = data, label = data),position = position_dodge(.9), vjust = -0.2, size = 5.3) +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), plot.title = element_text(face = 'bold', size=22), legend.title=element_text(size=22, face="bold"),legend.text=element_text(size=20))

pdf(file = 'nar_rev_diagram/g_purity_integrated_aggdf.pdf',width = 8,height = 6)
plot(g_purity_integrated_aggdf)#
dev.off()

purity_integrated_aggdf %>% group_by(dr) %>%
  select(value) %>% summarise_at(vars(value),              
                                 list(name = mean))       
#r, cdiff, csum parameter#################################################################################
#parameter r########################################
#leng parameter r
leng_connection_m <- test_leng@connect_cluster$cluster_connection
leng_point_possibility_df <- 
  tibble(leng_ds_1 = test_leng@point_possibility[["counted_possibility"]][["ds_1"]],
         leng_ds_2 = test_leng@point_possibility[["counted_possibility"]][["ds_2"]],
         leng_ds_3 = test_leng@point_possibility[["counted_possibility"]][["ds_3"]])



test_para_leng <- test_leng
test_para_leng <- point_possibility(test_para_leng, r = 2)
test_para_leng <- connect_cluster(test_para_leng, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)
test_para_leng@connect_cluster$cluster_connection

calculate_l_distance_r <- function(data, r){
  data <- point_possibility(data, r = r)
  data_tb <- tibble(leng_ds_1 = data@point_possibility[["counted_possibility"]][["ds_1"]],
                    leng_ds_2 = data@point_possibility[["counted_possibility"]][["ds_2"]],
                    leng_ds_3 = data@point_possibility[["counted_possibility"]][["ds_3"]])
  diff_tb <- leng_point_possibility_df - data_tb  
  distance <- diff_tb %>% mutate_all(function(x) x**2) %>% sum()
  return(distance)
}

distance_leng_r_1 <- calculate_l_distance_r(test_para_leng, r = 1)
distance_leng_r_1.5 <- calculate_l_distance_r(test_para_leng, r = 1.5)
distance_leng_r_2 <- calculate_l_distance_r(test_para_leng, r = 2)
distance_leng_r_2.5 <- calculate_l_distance_r(test_para_leng, r = 2.5)
distance_leng_r_3 <- calculate_l_distance_r(test_para_leng, r = 3)



avg_l_dis_tb <- tibble(l_r_1 = distance_leng_r_1/length(leng_index), 
                         l_r_1.5 = distance_leng_r_1.5/length(leng_index), 
                         l_r_2 = distance_leng_r_2/length(leng_index), 
                         l_r_2.5 = distance_leng_r_2.5/length(leng_index), 
                         l_r_3 = distance_leng_r_3/length(leng_index))

#nestorowa parameter r
nestorowa_connection_m <- test_nestorowa@connect_cluster$cluster_connection

nestorowa_point_possibility_df <- 
  tibble(nestorowa_ds_1 = test_nestorowa@point_possibility[["counted_possibility"]][["ds_1"]],
         nestorowa_ds_2 = test_nestorowa@point_possibility[["counted_possibility"]][["ds_2"]],
         nestorowa_ds_3 = test_nestorowa@point_possibility[["counted_possibility"]][["ds_3"]],
         nestorowa_ds_4 = test_nestorowa@point_possibility[["counted_possibility"]][["ds_4"]])



test_para_nestorowa <- test_nestorowa
test_para_nestorowa <- point_possibility(test_para_nestorowa, r = 2)
test_para_nestorowa <- connect_cluster(test_para_nestorowa, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)

test_para_nestorowa@connect_cluster$cluster_connection

calculate_n_distance_r <- function(data, r){
  data <- point_possibility(data, r = r)
  data_tb <- tibble(nestorowa_ds_1 = data@point_possibility[["counted_possibility"]][["ds_1"]],
                    nestorowa_ds_2 = data@point_possibility[["counted_possibility"]][["ds_2"]],
                    nestorowa_ds_3 = data@point_possibility[["counted_possibility"]][["ds_3"]],
                    nestorowa_ds_4 = data@point_possibility[["counted_possibility"]][["ds_4"]])
  diff_tb <- nestorowa_point_possibility_df - data_tb  
  distance <- diff_tb %>% mutate_all(function(x) x**2) %>% sum()
  return(distance)
}

distance_nestorowa_r_1 <- calculate_n_distance_r(test_para_nestorowa, r = 1)
distance_nestorowa_r_1.5 <- calculate_n_distance_r(test_para_nestorowa, r = 1.5)
distance_nestorowa_r_2 <- calculate_n_distance_r(test_para_nestorowa, r = 2)
distance_nestorowa_r_2.5 <- calculate_n_distance_r(test_para_nestorowa, r = 2.5)
distance_nestorowa_r_3 <- calculate_n_distance_r(test_para_nestorowa, r = 3)


avg_n_dis_tb <- tibble(l_r_1 = distance_nestorowa_r_1/length(nestorowa_index), 
                       l_r_1.5 = distance_nestorowa_r_1.5/length(nestorowa_index), 
                       l_r_2 = distance_nestorowa_r_2/length(nestorowa_index), 
                       l_r_2.5 = distance_nestorowa_r_2.5/length(nestorowa_index), 
                       l_r_3 = distance_nestorowa_r_3/length(nestorowa_index))

#kowalczyk parameter r
kowalczyk_connection_m <- test_kowalczyk@connect_cluster$cluster_connection
kowalczyk_point_possibility_df <- 
  tibble(kowalczyk_ds_1 = test_kowalczyk@point_possibility[["counted_possibility"]][["ds_1"]],
         kowalczyk_ds_2 = test_kowalczyk@point_possibility[["counted_possibility"]][["ds_2"]],
         kowalczyk_ds_3 = test_kowalczyk@point_possibility[["counted_possibility"]][["ds_3"]])



test_para_kowalczyk <- test_kowalczyk
test_para_kowalczyk <- point_possibility(test_para_kowalczyk, r = 2)
test_para_kowalczyk <- connect_cluster(test_para_kowalczyk, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)

test_para_kowalczyk@connect_cluster$cluster_connection

calculate_k_distance_r <- function(data, r){
  data <- point_possibility(data, r = r)
  data_tb <- tibble(kowalczyk_ds_1 = data@point_possibility[["counted_possibility"]][["ds_1"]],
                    kowalczyk_ds_2 = data@point_possibility[["counted_possibility"]][["ds_2"]],
                    kowalczyk_ds_3 = data@point_possibility[["counted_possibility"]][["ds_3"]])
  diff_tb <- kowalczyk_point_possibility_df - data_tb  
  distance <- diff_tb %>% mutate_all(function(x) x**2) %>% sum()
  return(distance)
}

distance_kowalczyk_r_1 <- calculate_k_distance_r(test_para_kowalczyk, r = 1)
distance_kowalczyk_r_1.5 <- calculate_k_distance_r(test_para_kowalczyk, r = 1.5)
distance_kowalczyk_r_2 <- calculate_k_distance_r(test_para_kowalczyk, r = 2)
distance_kowalczyk_r_2.5 <- calculate_k_distance_r(test_para_kowalczyk, r = 2.5)
distance_kowalczyk_r_3 <- calculate_k_distance_r(test_para_kowalczyk, r = 3)


avg_k_dis_tb <- tibble(l_r_1 = distance_kowalczyk_r_1/length(kowalczyk_index), 
                       l_r_1.5 = distance_kowalczyk_r_1.5/length(kowalczyk_index), 
                       l_r_2 = distance_kowalczyk_r_2/length(kowalczyk_index), 
                       l_r_2.5 = distance_kowalczyk_r_2.5/length(kowalczyk_index), 
                       l_r_3 = distance_kowalczyk_r_3/length(kowalczyk_index))

#camp parameter r
camp_connection_m <- test_camp@connect_cluster$cluster_connection
camp_point_possibility_df <- 
  tibble(camp_ds_1 = test_camp@point_possibility[["counted_possibility"]][["ds_1"]],
         camp_ds_2 = test_camp@point_possibility[["counted_possibility"]][["ds_2"]],
         camp_ds_3 = test_camp@point_possibility[["counted_possibility"]][["ds_3"]],
         camp_ds_4 = test_camp@point_possibility[["counted_possibility"]][["ds_4"]])



test_para_camp <- test_camp
test_para_camp <- point_possibility(test_para_camp, r = 2)
test_para_camp <- connect_cluster(test_para_camp, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)

test_para_camp@connect_cluster$cluster_connection

calculate_c_distance_r <- function(data, r){
  data <- point_possibility(data, r = r)
  data_tb <- tibble(camp_ds_1 = data@point_possibility[["counted_possibility"]][["ds_1"]],
                    camp_ds_2 = data@point_possibility[["counted_possibility"]][["ds_2"]],
                    camp_ds_3 = data@point_possibility[["counted_possibility"]][["ds_3"]],
                    camp_ds_4 = data@point_possibility[["counted_possibility"]][["ds_4"]])
  diff_tb <- camp_point_possibility_df - data_tb  
  distance <- diff_tb %>% mutate_all(function(x) x**2) %>% sum()
  return(distance)
}

distance_camp_r_1 <- calculate_c_distance_r(test_para_camp, r = 1)
distance_camp_r_1.5 <- calculate_c_distance_r(test_para_camp, r = 1.5)
distance_camp_r_2 <- calculate_c_distance_r(test_para_camp, r = 2)
distance_camp_r_2.5 <- calculate_c_distance_r(test_para_camp, r = 2.5)
distance_camp_r_3 <- calculate_c_distance_r(test_para_camp, r = 3)


avg_c_dis_tb <- tibble(l_r_1 = distance_camp_r_1/length(camp_index), 
                       l_r_1.5 = distance_camp_r_1.5/length(camp_index), 
                       l_r_2 = distance_camp_r_2/length(camp_index), 
                       l_r_2.5 = distance_camp_r_2.5/length(camp_index), 
                       l_r_3 = distance_camp_r_3/length(camp_index))

#yan parameter r
yan_connection_m <- test_yan@connect_cluster$cluster_connection
yan_point_possibility_df <- 
  tibble(yan_ds_1 = test_yan@point_possibility[["counted_possibility"]][["ds_1"]],
         yan_ds_2 = test_yan@point_possibility[["counted_possibility"]][["ds_2"]],
         yan_ds_3 = test_yan@point_possibility[["counted_possibility"]][["ds_3"]],
         yan_ds_4 = test_yan@point_possibility[["counted_possibility"]][["ds_4"]])



test_para_yan <- test_yan
test_para_yan <- point_possibility(test_para_yan, r = 2)
test_para_yan <- connect_cluster(test_para_yan, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01)

test_para_yan@connect_cluster$cluster_connection

calculate_y_distance_r <- function(data, r){
  data <- point_possibility(data, r = r)
  data_tb <- tibble(yan_ds_1 = data@point_possibility[["counted_possibility"]][["ds_1"]],
                    yan_ds_2 = data@point_possibility[["counted_possibility"]][["ds_2"]],
                    yan_ds_3 = data@point_possibility[["counted_possibility"]][["ds_3"]],
                    yan_ds_4 = data@point_possibility[["counted_possibility"]][["ds_4"]])
  diff_tb <- yan_point_possibility_df - data_tb  
  distance <- diff_tb %>% mutate_all(function(x) x**2) %>% sum()
  return(distance)
}

distance_yan_r_1 <- calculate_y_distance_r(test_para_yan, r = 1)
distance_yan_r_1.5 <- calculate_y_distance_r(test_para_yan, r = 1.5)
distance_yan_r_2 <- calculate_y_distance_r(test_para_yan, r = 2)
distance_yan_r_2.5 <- calculate_y_distance_r(test_para_yan, r = 2.5)
distance_yan_r_3 <- calculate_y_distance_r(test_para_yan, r = 3)


avg_y_dis_tb <- tibble(l_r_1 = distance_yan_r_1/length(yan_index), 
                       l_r_1.5 = distance_yan_r_1.5/length(yan_index), 
                       l_r_2 = distance_yan_r_2/length(yan_index), 
                       l_r_2.5 = distance_yan_r_2.5/length(yan_index), 
                       l_r_3 = distance_yan_r_3/length(yan_index))
#plot
r_value <- str_split(rownames(t(avg_l_dis_tb)), '_') %>% 
            lapply(function(x) x[3]) %>%
            unlist() %>%
            as.numeric() 


avg_l_dis_ggtb <- t(avg_l_dis_tb) %>% tibble() %>%
  add_column(r_parameter = r_value, data = as.factor(rep('leng', 5))
  ) %>% rename(distance=1)



avg_n_dis_ggtb <- t(avg_n_dis_tb) %>% tibble() %>%
  add_column(r_parameter = r_value, data = as.factor(rep('nestorowa', 5))
  ) %>% rename(distance=1)



avg_k_dis_ggtb <- t(avg_k_dis_tb) %>% tibble() %>%
  add_column(r_parameter = r_value, data = as.factor(rep('kowalczyk', 5))
  ) %>% rename(distance=1)



avg_c_dis_ggtb <- t(avg_c_dis_tb) %>% tibble() %>%
  add_column(r_parameter = r_value, data = as.factor(rep('camp', 5))
  ) %>% rename(distance=1)



avg_y_dis_ggtb <- t(avg_y_dis_tb) %>% tibble() %>%
  add_column(r_parameter = r_value, data = as.factor(rep('yan', 5))
  ) %>% rename(distance=1)




#add

avg_dis_aggtb <- avg_l_dis_ggtb %>% 
  bind_rows(avg_n_dis_ggtb) %>% 
  bind_rows(avg_k_dis_ggtb) %>% 
  bind_rows(avg_c_dis_ggtb) %>% 
  bind_rows(avg_y_dis_ggtb)


avg_dis_aggtb_plot <- ggplot() +  geom_line(data = avg_dis_aggtb, aes(r_parameter, distance, fill = data, color = data), size = 3) +
  # geom_point(data = total_dis_ggtb, aes(r_parameter, distance, fill = data, color = data), size = 5) +
  xlab('r value') + 
  ylab('Distance') + 
  ylim(0,0.1) +
  theme_classic() +
  # geom_text(data = purity_ggdf, aes(dr, value, fill = data, label = data),position = position_dodge(.9), vjust = -0.2, size = 5.3) +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_text(face = 'bold', size = 22), title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), plot.title = element_text(face = 'bold', size=22), legend.title=element_text(size=22, face="bold"),legend.text=element_text(size=20))

pdf(file = 'nar_rev_diagram/avg_dis_aggtb_plot.pdf',width = 8,height = 6)
plot(avg_dis_aggtb_plot)#
dev.off()
#parameter c_sum, c_diff########################################


leng_connection_mt <- test_para_leng@connect_cluster$cluster_connection
nestorowa_connection_mt <- test_para_nestorowa@connect_cluster$cluster_connection
kowalczyk_connection_mt <- test_para_kowalczyk@connect_cluster$cluster_connection
camp_connection_mt <- test_para_camp@connect_cluster$cluster_connection
yan_connection_mt <- test_para_yan@connect_cluster$cluster_connection

calculate_distance_csum_cdiff <- function(data, csum, cdiff, mt_cri){
  data <- connect_cluster(data, sum_cri = csum, diff_cri = cdiff)
  data <- data@connect_cluster$cluster_connection
  mt_cri <- mt_cri
  distance <- sum(mt_cri == data) / (dim(mt_cri)[1] ** 2)
  return(distance)
}




c_sum <- 0.8
c_diff <- 0.5
c_sum_range <- seq(c_sum-0.15, c_sum+0.15, 0.05)
c_diff_range <- seq(c_diff-0.15, c_diff+0.15, 0.05)

c_sum_diff_tb_l <- tibble(.rows = 1)
c_sum_diff_tb_n <- tibble(.rows = 1)
for (i in c_sum_range) {
  for (j in c_diff_range) {
    score_l <- calculate_distance_csum_cdiff(test_para_leng, csum = i, cdiff = j, mt_cri = leng_connection_mt)
    score_n <- calculate_distance_csum_cdiff(test_para_nestorowa, csum = i, cdiff = j, mt_cri = nestorowa_connection_mt)
    name_l <- paste0('leng_csum_', i, '_cdiff_', j)
    name_n <- paste0('nestorowa_csum_', i, '_cdiff_', j)

    c_sum_diff_tb_l <- c_sum_diff_tb_l %>% 
      add_column(!!name_l := score_l)
    c_sum_diff_tb_n <- c_sum_diff_tb_n %>% 
      add_column(!!name_n := score_n)
  }
}

#add
c_sum_diff_tb_k <- tibble(.rows = 1)
c_sum_diff_tb_c <- tibble(.rows = 1)
c_sum_diff_tb_y <- tibble(.rows = 1)
#k
for (i in c_sum_range) {
  for (j in c_diff_range) {
    score_k <- calculate_distance_csum_cdiff(test_para_kowalczyk, csum = i, cdiff = j, mt_cri = kowalczyk_connection_mt)
    name_k <- paste0('kowalczyk_csum_', i, '_cdiff_', j)

    c_sum_diff_tb_k <- c_sum_diff_tb_k %>% 
      add_column(!!name_k := score_k)
    
  }
}
#c,y
for (i in c_sum_range) {
  for (j in c_diff_range) {
    score_c <- calculate_distance_csum_cdiff(test_para_camp, csum = i, cdiff = j, mt_cri = camp_connection_mt)
    score_y <- calculate_distance_csum_cdiff(test_para_yan, csum = i, cdiff = j, mt_cri = yan_connection_mt)
    
    name_c <- paste0('camp_csum_', i, '_cdiff_', j)
    name_y <- paste0('yan_csum_', i, '_cdiff_', j)
    
    c_sum_diff_tb_c <- c_sum_diff_tb_c %>% 
      add_column(!!name_c := score_c)
    c_sum_diff_tb_y <- c_sum_diff_tb_y %>% 
      add_column(!!name_y := score_y)
    
  }
}
# plot
c_sum_diff_tb_l_index <- colnames(c_sum_diff_tb_l) %>% str_split('_')
c_sum_diff_tb_l_ggtb <- c_sum_diff_tb_l %>% t() %>% as_tibble() %>%
  rename(value = 1) %>% add_column(data = c_sum_diff_tb_l_index %>% 
                                     lapply(function(x) x[1]) %>% unlist(),
                                   C_sum = c_sum_diff_tb_l_index %>% 
                                      lapply(function(x) x[3]) %>% unlist() %>% as.numeric(),
                                   C_diff = c_sum_diff_tb_l_index %>% 
                                     lapply(function(x) x[5]) %>% unlist() %>% as.numeric())
c_sum_diff_tb_n_index <- colnames(c_sum_diff_tb_n) %>% str_split('_')
c_sum_diff_tb_n_ggtb <- c_sum_diff_tb_n %>% t() %>% as_tibble() %>%
  rename(value = 1) %>% add_column(data = c_sum_diff_tb_n_index %>% 
                                     lapply(function(x) x[1]) %>% unlist(),
                                   C_sum = c_sum_diff_tb_n_index %>% 
                                     lapply(function(x) x[3]) %>% unlist() %>% as.numeric(),
                                   C_diff = c_sum_diff_tb_n_index %>% 
                                     lapply(function(x) x[5]) %>% unlist() %>% as.numeric())
c_sum_diff_tb_k_index <- colnames(c_sum_diff_tb_k) %>% str_split('_')
c_sum_diff_tb_k_ggtb <- c_sum_diff_tb_k %>% t() %>% as_tibble() %>%
  rename(value = 1) %>% add_column(data = c_sum_diff_tb_k_index %>% 
                                     lapply(function(x) x[1]) %>% unlist(),
                                   C_sum = c_sum_diff_tb_k_index %>% 
                                     lapply(function(x) x[3]) %>% unlist() %>% as.numeric(),
                                   C_diff = c_sum_diff_tb_k_index %>% 
                                     lapply(function(x) x[5]) %>% unlist() %>% as.numeric())
c_sum_diff_tb_c_index <- colnames(c_sum_diff_tb_c) %>% str_split('_')
c_sum_diff_tb_c_ggtb <- c_sum_diff_tb_c %>% t() %>% as_tibble() %>%
  rename(value = 1) %>% add_column(data = c_sum_diff_tb_c_index %>% 
                                     lapply(function(x) x[1]) %>% unlist(),
                                   C_sum = c_sum_diff_tb_c_index %>% 
                                     lapply(function(x) x[3]) %>% unlist() %>% as.numeric(),
                                   C_diff = c_sum_diff_tb_c_index %>% 
                                     lapply(function(x) x[5]) %>% unlist() %>% as.numeric())
c_sum_diff_tb_y_index <- colnames(c_sum_diff_tb_y) %>% str_split('_')
c_sum_diff_tb_y_ggtb <- c_sum_diff_tb_y %>% t() %>% as_tibble() %>%
  rename(value = 1) %>% add_column(data = c_sum_diff_tb_y_index %>% 
                                     lapply(function(x) x[1]) %>% unlist(),
                                   C_sum = c_sum_diff_tb_y_index %>% 
                                     lapply(function(x) x[3]) %>% unlist() %>% as.numeric(),
                                   C_diff = c_sum_diff_tb_y_index %>% 
                                     lapply(function(x) x[5]) %>% unlist() %>% as.numeric())


c_sum_diff_tb_l_ggtb <- c_sum_diff_tb_l_ggtb %>% mutate(value = 1-value)
c_sum_diff_tb_n_ggtb <- c_sum_diff_tb_n_ggtb %>% mutate(value = 1-value)
c_sum_diff_tb_k_ggtb <- c_sum_diff_tb_k_ggtb %>% mutate(value = 1-value)
c_sum_diff_tb_c_ggtb <- c_sum_diff_tb_c_ggtb %>% mutate(value = 1-value)
c_sum_diff_tb_y_ggtb <- c_sum_diff_tb_y_ggtb %>% mutate(value = 1-value)

plot_c_sum_diff_tb <- function(data){
  a <- ggplot(data = data, aes(x = C_sum, y = C_diff)) +
    geom_tile(colour="black", size = 1.5, aes(fill = value, width=0.05, height=0.05))+
    geom_text(aes(label = value), color = "white", fontface = "bold", size = 6)+
    scale_x_continuous(expand = c(0,0), breaks = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95)) +
    
    scale_y_continuous(expand = c(0,0), breaks = c(0.35,0.4,0.45,0.5,0.55,0.6,0.65)) +
    # scale_fill_gradientn("Distance",limits = c(0, 1), colours = rainbow(3)) +
    # scale_fill_gradientn("Distance",limits = c(0, 1), colours = topo.colors(5)) +
    viridis::scale_fill_viridis("Distance",limits = c(0, 1), option = 'D') +
    # scale_colour_viridis_d("Distance",limits = c(0, 1)) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5), 
          axis.title.x = element_text(face="bold", size = 20),
          axis.title.y = element_text(face="bold", size = 20),
          axis.text = element_text(face="bold", size = 20),
          legend.title = element_text(face="bold", size = 20),
          legend.text=element_text(size=20))
  return(a)
}


c_sum_diff_tb_l_ggtb_plot <- plot_c_sum_diff_tb(c_sum_diff_tb_l_ggtb %>% mutate(value = round(value, digits = 2)))
c_sum_diff_tb_n_ggtb_plot <- plot_c_sum_diff_tb(c_sum_diff_tb_n_ggtb %>% mutate(value = round(value, digits = 2)))
c_sum_diff_tb_k_ggtb_plot <- plot_c_sum_diff_tb(c_sum_diff_tb_k_ggtb %>% mutate(value = round(value, digits = 2)))
c_sum_diff_tb_c_ggtb_plot <- plot_c_sum_diff_tb(c_sum_diff_tb_c_ggtb %>% mutate(value = round(value, digits = 2)))
c_sum_diff_tb_y_ggtb_plot <- plot_c_sum_diff_tb(c_sum_diff_tb_y_ggtb %>% mutate(value = round(value, digits = 2)))

pdf(file = 'nar_rev_diagram/c_sum_diff_tb_l_ggtb_plot.pdf',width = 7,height = 6)
plot(c_sum_diff_tb_l_ggtb_plot)
dev.off()

pdf(file = 'nar_rev_diagram/c_sum_diff_tb_n_ggtb_plot.pdf',width = 7,height = 6)
plot(c_sum_diff_tb_n_ggtb_plot)
dev.off()

pdf(file = 'nar_rev_diagram/c_sum_diff_tb_k_ggtb_plot.pdf',width = 7,height = 6)
plot(c_sum_diff_tb_k_ggtb_plot)
dev.off()

pdf(file = 'nar_rev_diagram/c_sum_diff_tb_c_ggtb_plot.pdf',width = 7,height = 6)
plot(c_sum_diff_tb_c_ggtb_plot)
dev.off()

pdf(file = 'nar_rev_diagram/c_sum_diff_tb_y_ggtb_plot.pdf',width = 7,height = 6)
plot(c_sum_diff_tb_y_ggtb_plot)
dev.off()

#mutrans############################################################################
#output ref index
#simulated
data_name <-  c('bc_0','bc_2','bc_9', 'bt_0','bt_2','bt_5', 'l_0','l_1','l_2', 'd_0','d_2','d_9', 'c_0','c_2','c_3')
file_name <- c('out_bifurcating_cycle','bc_2','bc_9', 'out_binary_tree','bt_2','bt_5', 'out_linear','l_1','l_2', 'out_disconnected','d_2','d_9', 'out_cycle','c_2','c_3')


for (i in file_name) {
  name <- i
  assign(name, readRDS(paste0(i, '.rds'))) 
  index_name <- paste0(i, '$dataset$progressions$from')
  index <- eval(parse(text = index_name))
  write.csv(index, paste0(data_name[which(file_name == i)], '_index.csv'))
}

#real index to csv
leng_index <- readRDS('hescsmsinfo_cycle_index_for_new.rds')
nestorowa_index <- readRDS('nestorowa_index_batch1_type_filter_no_stem_for_new.rds')
kowalczyk_index <- readRDS('stem_mouse_C57BL6_index_old_for_new.rds')
yan_index <- readRDS('yan_index_for_new.rds')
camp_index <- readRDS('camp1_index_batch1_for_new.rds')
write.csv(leng_index, 'leng_index.csv')
write.csv(nestorowa_index, 'nestorowa_index.csv')
write.csv(kowalczyk_index, 'kowalczyk_index.csv')
write.csv(yan_index, 'yan_index.csv')
write.csv(camp_index, 'camp_index.csv')

# calculate purity using inferred index########################
purity_df_simulated <- readRDS('purity_df_simulated.rds')
data_name_mutrans_simulated <-  c('bc_0','bc_2','bc_9', 'bt_0','bt_2','bt_5', 'l_0','l_1','l_2', 'd_0','d_2','d_9', 'c_0','c_2','c_3')

#simulated
for (i in data_name_mutrans_simulated) {
  #mutrans index
  a <- read_table(paste0('mutrans/mutrans_', i, '_index.csv')) %>% 
    colnames() %>% str_split(',') %>% unlist() %>% as.factor()
  # ref index
  b <- read.csv(paste0(i, '_index.csv'))[,2]
  purity_df_simulated[purity_df_simulated$.id == i, 'model_mutrans'] <- 
    NMF::purity(a, b)
  
}
saveRDS(purity_df_simulated, 'purity_df_simulated_nar_rev.rds')

#real
leng_index <- leng_index 
nestorowa_index <- nestorowa_index %>% as.factor()
kowalczyk_index <- kowalczyk_index %>% as.factor()
yan_index <- yan_index %>% c() %>% unlist() %>% as.factor()
camp_index <- camp_index %>% c() %>% unlist() %>% as.factor()

data_name_mutrans_real <- c('leng', 'nestorowa', 'kowalczyk_o', 'yan', 'camp')
purity_df_real <- readRDS('purity_df_real.rds')

for (i in data_name_mutrans_real) {
  #mutrans index
  a <- read_table(paste0('mutrans/mutrans_', i, '_data_index.csv')) %>% 
    colnames() %>% str_split(',') %>% unlist() %>% as.factor()
  # ref index
  if (i == 'kowalczyk_o') {
    b <- eval(parse(text = paste0('kowalczyk', '_index')))
  } else {
    b <- eval(parse(text = paste0(i, '_index')))
  }
  purity_df_real[purity_df_real$.id == paste0(i,'_data'), 'model_mutrans'] <- 
    NMF::purity(a, b)
}

saveRDS(purity_df_real, 'purity_df_real_nar_rev.rds')

#summary data
purity_df_simulated <- purity_df_simulated %>% 
  add_column(topology = unlist(lapply(str_split(purity_df_simulated$.id, '_'), function(x) x[1]))) 

purity_df_simulated$model_monocle <- as.numeric(purity_df_simulated$model_monocle)
purity_df_simulated$model_dbcti <- as.numeric(purity_df_simulated$model_dbcti)
purity_df_real$model_monocle <- as.numeric(purity_df_real$model_monocle)
purity_df_real$model_dbcti <- as.numeric(purity_df_real$model_dbcti)

purity_df_simulated_summary <- purity_df_simulated %>% group_by(topology) %>% 
  summarise(slingshot = mean(model_slingshot), paga = mean(model_paga),tscan = mean(model_tscan),monocle2 = mean(model_monocle),mutrans = mean(model_mutrans),dbcti = mean(model_dbcti))
purity_df_real_summary <- purity_df_real %>% group_by(.id) %>% 
  summarise(slingshot = mean(model_slingshot), paga = mean(model_paga),tscan = mean(model_tscan),monocle2 = mean(model_monocle),mutrans = mean(model_mutrans),dbcti = mean(model_dbcti))

purity_df_real_summary$.id <- c('Camp_data', 'Kowalczyk_data', 'Leng_data', 'Nestorowa_data', 'Yan_data')

saveRDS(purity_df_simulated_summary, 'purity_df_simulated_summary_nar_rev.rds')
saveRDS(purity_df_real_summary, 'purity_df_real_summary_nar_rev.rds')
#plot metrics ########################
purity_simulated_ggplot <- purity_df_simulated_summary %>% pivot_longer(cols = !topology,names_to = 'Tool', values_to = 'Score')
purity_real_ggplot <- purity_df_real_summary %>% pivot_longer(cols = !.id,names_to = 'Tool', values_to = 'Score')


purity_simulated_score_boxplot <- ggplot(purity_simulated_ggplot, aes(x = Tool, y = Score, fill = Tool), color = 'black')+
  geom_boxplot() + 
  theme_classic() +
  ylab('Score') +
  ggtitle(label = 'Purity score boxplot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), legend.position = 'none', plot.title = element_text(face = 'bold', size=22))
pdf(file = 'purity_simulated_score_boxplot_nar_rev.pdf',width = 11,height = 7)
plot(purity_simulated_score_boxplot)
dev.off()

purity_real_score_boxplot <- ggplot(purity_real_ggplot, aes(x = Tool, y = Score, fill = Tool), color = 'black')+
  geom_boxplot() + 
  theme_classic() +
  ylab('Score') +
  ggtitle(label = 'Purity score boxplot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), legend.position = 'none', plot.title = element_text(face = 'bold', size=22))

pdf(file = 'purity_real_score_boxplot_nar_rev.pdf',width = 11,height = 7)
plot(purity_real_score_boxplot)
dev.off()





