setwd('..')
library(dynplot)
library(tidyverse)
library(dynwrap)
library(dbcti)


simulated_data_tool <- readRDS('datasets/simulated_data_tool.rds')
simulated_data_tool_1 <- readRDS('datasets/simulated_data_tool_1.rds')
real_data_tool <- readRDS('datasets/real_data_tool.rds')
real_data_tool_1 <- readRDS('datasets/real_data_tool_1.rds')



#plot

plot_data <- function(data){
  for (i in names(data)){
    for(j in names(data[[i]])){
      try({
        pdf(file = paste0('diagram/', i, '_', j, '.pdf'),width = 6,height = 6)
        print(plot_dimred(data[[i]][[j]]))
        dev.off()
      }
      )
      
    }
  }
}

for(j in names(data[[i]])){
  print(j)
}


model_paga_1 <- read_rds('datasets/model_paga_1.rds')

pdf(file = paste0('diagram/', 'camp_data', '_', 'model_paga', '.pdf'),width = 6,height = 6)
plot_dimred(model_paga_1)
dev.off()


plot_data(simulated_data_tool)
plot_data(simulated_data_tool_1)
plot_data(real_data_tool)
plot_data(real_data_tool_1)

#combine list for metric evaluation
combined_list_simulated <- simulated_data_tool
for (i in names(combined_list_simulated)) {
  combined_list_simulated[[i]][['model_tscan']] <- simulated_data_tool_1[[i]][['model_tscan']]
  combined_list_simulated[[i]] <- combined_list_simulated[[i]][-3]
}
combined_list_real <- real_data_tool
for (i in names(combined_list_real)) {
  combined_list_real[[i]][['model_tscan']] <- real_data_tool_1[[i]][['model_tscan']]
  combined_list_real[[i]] <- combined_list_real[[i]][-3]
}

saveRDS(combined_list_simulated, 'combined_list_simulated.rds')
saveRDS(combined_list_real, 'combined_list_real.rds')

#get the index for metric
yan_data_monocle <- readRDS('datasets/yan_data_monocle.rds')
nestorowa_data_monocle <- readRDS('datasets/nestorowa_data_batch1_type_filter_no_stem_monocle.rds')
camp_data_monocle <- readRDS('datasets/camp1_data_batch1_monocle.rds')
kowalczyk_data_monocle <- readRDS('datasets/stem_mouse_C57BL6_old_monocle.rds')
leng_data_monocle <- readRDS('datasets/hesc_monocle.rds')

yan_index <- yan_data_monocle@phenoData@data[["index"]]
nestorowa_index <- nestorowa_data_monocle@phenoData@data[["index"]]
camp_index <- camp_data_monocle@phenoData@data[["index"]]
kowalczyk_index <- kowalczyk_data_monocle@phenoData@data[["index"]]
leng_index <- leng_data_monocle@phenoData@data[["index"]]

saveRDS(yan_index, 'datasets/yan_index.rds')
saveRDS(nestorowa_index, 'datasets/nestorowa_index.rds')
saveRDS(camp_index, 'datasets/camp_index.rds')
saveRDS(kowalczyk_index, 'datasets/kowalczyk_index.rds')
saveRDS(leng_index, 'datasets/leng_index.rds')

###################################################################
purity_df_simulated <- readRDS('datasets/purity_df_simulated.rds')
purity_df_real <- readRDS('datasets/purity_df_real.rds')

simulated_purity_score_old <- readRDS('datasets/simulated_purity_score_old')
real_purity_score_old <- readRDS('datasets/real_purity_score_old')



simulated_monocle <- vector(mode = 'list', length = 10)
names(simulated_monocle) <- c('bc_2','bc_3','bt_2','bt_3','l_2','l_3','c_2','c_3','d_2','d_3')
simulated_dbcti <- simulated_monocle
for (i in names(simulated_monocle)) {
  simulated_monocle[[i]] <- eval(parse(text = paste0(i,'_monocle')))
}
for (i in names(simulated_dbcti)) {
  simulated_dbcti[[i]] <- eval(parse(text = paste0(i,'_dbcti')))
}


calculate_purity_monocle <- function(x, ref){
  x <- as.factor(x@phenoData@data[["State"]])
  ref <- as.factor(ref$dataset$progressions$from)
  purity <- NMF::purity(x, ref)
  return(purity)
}
calculate_purity_dbcti <- function(x, ref){
  x <- as.factor(x@distribution_estimation$cluster_index)
  ref <- as.factor(ref$dataset$progressions$from)
  purity <- NMF::purity(x, ref)
  return(purity)
}

simulated_dbcti_purity <- simulated_dbcti
simulated_monocle_purity <- simulated_monocle

for (i in names(simulated_dbcti)) {
  ref <- eval(parse(text = i))
  simulated_dbcti_purity[[i]] <- calculate_purity_dbcti(simulated_dbcti[[i]], ref)
  simulated_monocle_purity[[i]] <- calculate_purity_monocle(simulated_monocle[[i]], ref)
}


#simulated
simulated_dbcti_purity
simulated_monocle_purity
purity_df_simulated
simulated_purity_score_old

#real
purity_df_real
real_purity_score_old


##############################################combien score
simulated_purity_score_old <- simulated_purity_score_old[, c(1,4,5,8,9,12,13,16,17,20)]
real_purity_score_old <- real_purity_score_old[, c(1,4,5,8,9,12,17,20,21,24)]


purity_df_simulated <- purity_df_simulated %>% add_column(model_monocle = NA, model_dbcti = NA)
purity_df_real <- purity_df_real %>% add_column(model_monocle = NA, model_dbcti = NA)
purity_df_real$model_monocle[1] <- real_purity_score_old[1,1]
purity_df_real$model_monocle[2] <- real_purity_score_old[1,7]
purity_df_real$model_monocle[3] <- real_purity_score_old[1,3]
purity_df_real$model_monocle[4] <- real_purity_score_old[1,9]
purity_df_real$model_monocle[5] <- real_purity_score_old[1,5]

purity_df_real$model_dbcti[1] <- real_purity_score_old[1,2]
purity_df_real$model_dbcti[2] <- real_purity_score_old[1,8]
purity_df_real$model_dbcti[3] <- real_purity_score_old[1,4]
purity_df_real$model_dbcti[4] <- real_purity_score_old[1,10]
purity_df_real$model_dbcti[5] <- real_purity_score_old[1,6]

#simulated
purity_df_simulated$model_monocle[1] <- simulated_purity_score_old[1,9]
purity_df_simulated$model_monocle[4] <- simulated_purity_score_old[1,5]
purity_df_simulated$model_monocle[7] <- simulated_purity_score_old[1,1]
purity_df_simulated$model_monocle[10] <- simulated_purity_score_old[1,3]
purity_df_simulated$model_monocle[13] <- simulated_purity_score_old[1,7]

purity_df_simulated$model_dbcti[1] <- simulated_purity_score_old[1,10]
purity_df_simulated$model_dbcti[4] <- simulated_purity_score_old[1,6]
purity_df_simulated$model_dbcti[7] <- simulated_purity_score_old[1,2]
purity_df_simulated$model_dbcti[10] <- simulated_purity_score_old[1,4]
purity_df_simulated$model_dbcti[13] <- simulated_purity_score_old[1,8]


#############
purity_df_simulated$model_monocle[c(2,3)] <- c(simulated_monocle_purity[[1]], simulated_monocle_purity[[2]])
purity_df_simulated$model_monocle[c(5,6)] <- c(simulated_monocle_purity[[3]], simulated_monocle_purity[[4]])
purity_df_simulated$model_monocle[c(8,9)] <- c(simulated_monocle_purity[[5]], simulated_monocle_purity[[6]])
purity_df_simulated$model_monocle[c(11,12)] <- c(simulated_monocle_purity[[7]], simulated_monocle_purity[[8]])
purity_df_simulated$model_monocle[c(14,15)] <- c(simulated_monocle_purity[[9]], simulated_monocle_purity[[10]])

purity_df_simulated$model_dbcti[c(2,3)] <- c(simulated_dbcti_purity[[1]], simulated_dbcti_purity[[2]])
purity_df_simulated$model_dbcti[c(5,6)] <- c(simulated_dbcti_purity[[3]], simulated_dbcti_purity[[4]])
purity_df_simulated$model_dbcti[c(8,9)] <- c(simulated_dbcti_purity[[5]], simulated_dbcti_purity[[6]])
purity_df_simulated$model_dbcti[c(11,12)] <- c(simulated_dbcti_purity[[7]], simulated_dbcti_purity[[8]])
purity_df_simulated$model_dbcti[c(14,15)] <- c(simulated_dbcti_purity[[9]], simulated_dbcti_purity[[10]])

saveRDS(purity_df_simulated, 'datasets/purity_df_simulated.rds')
saveRDS(purity_df_real, 'datasets/purity_df_real.rds')
purity_df_simulated <- purity_df_simulated %>% 
  add_column(topology = unlist(lapply(str_split(purity_df_simulated$.id, '_'), function(x) x[1]))) 

purity_df_simulated$model_monocle <- as.numeric(purity_df_simulated$model_monocle)
purity_df_simulated$model_dbcti <- as.numeric(purity_df_simulated$model_dbcti)
purity_df_real$model_monocle <- as.numeric(purity_df_real$model_monocle)
purity_df_real$model_dbcti <- as.numeric(purity_df_real$model_dbcti)


purity_df_simulated_summary <- purity_df_simulated %>% group_by(topology) %>% 
  summarise(slingshot = mean(model_slingshot), paga = mean(model_paga),tscan = mean(model_tscan),monocle2 = mean(model_monocle),dbcti = mean(model_dbcti))
purity_df_real_summary <- purity_df_real %>% group_by(.id) %>% 
  summarise(slingshot = mean(model_slingshot), paga = mean(model_paga),tscan = mean(model_tscan),monocle2 = mean(model_monocle),dbcti = mean(model_dbcti))

purity_df_real_summary$.id <- c('Camp_data', 'Kowalczyk_data', 'Leng_data', 'Nestorowa_data', 'Yan_data')


saveRDS(purity_df_simulated_summary, 'datasets/purity_df_simulated_summary.rds')
saveRDS(purity_df_real_summary, 'datasets/purity_df_real_summary.rds')

purity_simulated_ggplot <- purity_df_simulated_summary %>% pivot_longer(cols = !topology,names_to = 'Tool', values_to = 'Score')
purity_real_ggplot <- purity_df_real_summary %>% pivot_longer(cols = !.id,names_to = 'Tool', values_to = 'Score')


purity_simulated_score_boxplot <- ggplot(purity_simulated_ggplot, aes(x = Tool, y = Score, fill = Tool), color = 'black')+
  geom_boxplot() + 
  theme_classic() +
  ylab('Score') +
  ggtitle(label = 'Purity score boxplot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), legend.position = 'none', plot.title = element_text(face = 'bold', size=22))
pdf(file = 'purity_simulated_score_boxplot.pdf',width = 8,height = 6)
plot(purity_simulated_score_boxplot)
dev.off()

purity_real_score_boxplot <- ggplot(purity_real_ggplot, aes(x = Tool, y = Score, fill = Tool), color = 'black')+
  geom_boxplot() + 
  theme_classic() +
  ylab('Score') +
  ggtitle(label = 'Purity score boxplot') +
  theme(axis.text.x = element_text(face = 'bold', size = 22), axis.title.y = element_text(face = 'bold', size = 22), axis.title.x = element_blank(),title = element_text(face = 'bold'),
        axis.text.y = element_text(face = 'bold', size = 22), legend.position = 'none', plot.title = element_text(face = 'bold', size=22))

pdf(file = 'purity_real_score_boxplot.pdf',width = 8,height = 6)
plot(purity_real_score_boxplot)
dev.off()

