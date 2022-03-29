# hescmt.rds
# hescsmsinfo_cycle_index.rds
# stem_mouse_C57BL6_data_young.rds
# stem_mouse_C57BL6_index_young.rds
# stem_mouse_C57BL6_data_old.rds
# stem_mouse_C57BL6_index_old.rds
# camp1_data_batch1.rds
# camp1_index_batch1.rds
# nestorowa_data_batch1_type_filter_no_stem.rds
# nestorowa_index_batch1_type_filter_no_stem.rds
# yan_data.rds
# yan_index.rds


#data hesc#########################
hescmt<-read.csv('datasets/GSE64016_H1andFUCCI_normalized_EC_hesc.csv',row.names = 1)
hescsm<-read.csv('datasets/GSE64016_series_matrix_hesc.txt',sep = '\t',col.names = paste('v',1:600,sep = ''),fill = TRUE,header = FALSE)

#fiilter out gene with no expression
hescmt<-hescmt[apply(hescmt, 1, sum)!=0,]

#hescsmname<-unlist(strsplit(as.matrix(hescsm[12,2]),' '))
hescsmname<-as.vector(as.matrix(hescsm[c(27),c(-1,-462:-600)]))
hescsminfo<-hescsm[c(34,36,37),c(-1,-462:-600)]
colnames(hescsminfo)<-hescsmname
saveRDS(hescmt,'datasets/hescmt.rds')
saveRDS(hescsminfo,'datasets/hescsminfo.rds')
hescsmsinfo_cycle_index<-as.factor(unlist(lapply(str_split(colnames(hescmt),'_'), function(x) x[1])))
saveRDS(hescsmsinfo_cycle_index,'datasets/hescsmsinfo_cycle_index.rds')



#stem cell(mouse) (normalized value) (geo59114)#####
#C57BL6
stem_mouse_C57BL6_data<-read.csv('datasets/GSE59114_C57BL6_GEO_all.csv',sep = '\t')
gene_name_C57BL6 = stem_mouse_C57BL6_data[,1]
gene_name_C57BL6 <- gsub('^.|.$', '', gene_name_C57BL6)
d_index = duplicated(gene_name_C57BL6)

stem_mouse_C57BL6_data = stem_mouse_C57BL6_data[!d_index, ]
gene_name_C57BL6 = gene_name_C57BL6[!d_index]

ng_index <- str_detect(gene_name_C57BL6, '\\.')
gene_name_C57BL6 = gene_name_C57BL6[!ng_index]
stem_mouse_C57BL6_data = stem_mouse_C57BL6_data[!ng_index, ]
stem_mouse_C57BL6_data = stem_mouse_C57BL6_data[,-1]
rownames(stem_mouse_C57BL6_data) = gene_name_C57BL6

#index
stem_mouse_C57BL6_index = colnames(stem_mouse_C57BL6_data)
stem_mouse_C57BL6_index_age = unlist(lapply(str_split(stem_mouse_C57BL6_index, '_'), function(x) x[1]))
stem_mouse_C57BL6_index_cell_type = unlist(lapply(str_split(stem_mouse_C57BL6_index, '_'), function(x) x[2]))
stem_mouse_C57BL6_index = rbind(stem_mouse_C57BL6_index_age, stem_mouse_C57BL6_index_cell_type)


#young
stem_mouse_C57BL6_data_young = stem_mouse_C57BL6_data[,stem_mouse_C57BL6_index[1,]=='young']
stem_mouse_C57BL6_index_young = stem_mouse_C57BL6_index[2,stem_mouse_C57BL6_index[1,]=='young']
saveRDS(stem_mouse_C57BL6_data_young,'datasets/stem_mouse_C57BL6_data_young.rds')
saveRDS(stem_mouse_C57BL6_index_young,'datasets/stem_mouse_C57BL6_index_young.rds')

#old
stem_mouse_C57BL6_data_old = stem_mouse_C57BL6_data[,stem_mouse_C57BL6_index[1,]=='old']
stem_mouse_C57BL6_index_old = stem_mouse_C57BL6_index[2,stem_mouse_C57BL6_index[1,]=='old']
saveRDS(stem_mouse_C57BL6_data_old,'datasets/stem_mouse_C57BL6_data_old.rds')
saveRDS(stem_mouse_C57BL6_index_old,'datasets/stem_mouse_C57BL6_index_old.rds')






#camp1 human liver(GSE81252) (hemberg lab github)#########################
camp1 <- readRDS('datasets/camp1.rds')
camp1_data = assay(camp1)
camp1_index = data.frame(t(data.frame(camp1@colData@listData)), stringsAsFactors = FALSE)

camp1_data_batch1 = camp1_data[,camp1_index[5, ] == 1]
camp1_index_batch1 = camp1_index[,camp1_index[5, ] == 1]



saveRDS(camp1_data_batch1,'datasets/camp1_data_batch1.rds')
saveRDS(camp1_index_batch1,'datasets/camp1_index_batch1.rds')




#Nestorowa mouse embryo(SRP041736) (hemberg lab github)#########################
nestorowa <- readRDS('datasets/nestorowa.rds')
nestorowa_data = assay(nestorowa)
nestorowa_index = data.frame(t(data.frame(nestorowa@colData@listData)), stringsAsFactors = FALSE)
nestorowa_index = nestorowa_index[c(1,5,6), ]
nestorowa_index = apply(nestorowa_index, 2, as.character)
nestorowa_batch_na = is.na(nestorowa_index[1, ])
                           
nestorowa_index = nestorowa_index[, !nestorowa_batch_na]
nestorowa_data = nestorowa_data[, !nestorowa_batch_na]

nestorowa_data_batch1 = nestorowa_data[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High']
nestorowa_data_batch2 = nestorowa_data[, nestorowa_index[1, ] == 2 & nestorowa_index[3, ] == 'High']
nestorowa_index_batch1 = nestorowa_index[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High']
nestorowa_index_batch2 = nestorowa_index[, nestorowa_index[1, ] == 2 & nestorowa_index[3, ] == 'High']


  
#filter
nestorowa_type_filter_index = nestorowa_index[2, ] == 'Unknown' | nestorowa_index[2, ] == 'ESLAM' | nestorowa_index[2, ] == 'LTHSC'
nestorowa_data_batch1_type_filter = nestorowa_data[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_index]
nestorowa_data_batch2_type_filter = nestorowa_data[, nestorowa_index[1, ] == 2 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_index]
nestorowa_index_batch1_type_filter = nestorowa_index[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_index]
nestorowa_index_batch2_type_filter = nestorowa_index[, nestorowa_index[1, ] == 2 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_index]



#no stem
nestorowa_type_filter_no_stem_index = nestorowa_index[2, ] == 'Unknown' | nestorowa_index[2, ] == 'ESLAM' | nestorowa_index[2, ] == 'LTHSC' | nestorowa_index[2, ] == 'STHSC'
nestorowa_data_batch1_type_filter_no_stem = nestorowa_data[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_no_stem_index]
nestorowa_index_batch1_type_filter_no_stem = nestorowa_index[, nestorowa_index[1, ] == 1 & nestorowa_index[3, ] == 'High' & !nestorowa_type_filter_no_stem_index]

saveRDS(nestorowa_data_batch1_type_filter_no_stem,'datasets/nestorowa_data_batch1_type_filter_no_stem.rds')
saveRDS(nestorowa_index_batch1_type_filter_no_stem,'datasets/nestorowa_index_batch1_type_filter_no_stem.rds')


#yan human embryo(gse36552) (hemberg lab github) (normalized)#########################
yan <- readRDS('datasets/yan.rds')
yan_data = assay(yan)
yan_index = data.frame(t(data.frame(yan@colData@listData)), stringsAsFactors = FALSE)
saveRDS(yan_data,'datasets/yan_data.rds')
saveRDS(yan_index,'datasets/yan_index.rds')




