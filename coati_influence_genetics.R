#seeing whether influence scores correlate with the genetics matrix - i.e. are individuals who are more influential on one another closer relatives?

library(sna)  
library(asnipe) 
library(abind) 
library(fields) 
library(ecodist)
library(vegan)


#directories 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/ch1_cleancode/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level2/' 

#read in library of functions
setwd(code_dir) 
source('coati_function_library_V1.R') 

#load in data
load('C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_relatedness.Rdata') 
load('C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/presedente_coati_ids.Rdata')

#removing males and Wildflower 
coati_ids_cut <- coati_ids[-c(5, 8, 9, 10, 18, 21),] 

#order of individuals in the genetics matrices (so same as subgrouping matrices)
pres_inds_neworder <- c("B12.L001", "B17.L001", "B16.L001", "B09.L001", "B03.L001", "B01.L001", "B02.L001", "B04.L001", "B14.L001", "B10.L001", "B07.L001", "B05.L001", "B06.L001", "B13.L001", "B11.L001", "B15.L001") 
pres_neworder_indx <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14) #order for coati_ids when adult males are not in the group 


#load csv from Jack's influence metrics
influence <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/presidente_speed_influence_matrix_2024_09_04.csv', header = T)
#remove first column
influence <- influence[,-1]
inf_matrix <- as.matrix(influence) 
rownames(inf_matrix) <- coati_ids_cut$name
colnames(inf_matrix) <- coati_ids_cut$name

n_inds = 16
#visualize_network_matrix_presedente(inf_matrix, coati_ids_cut[1:16,])

#reorder to same as genetics/subgrouping matrices
inf_matrix_raw <- inf_matrix[pres_neworder_indx, pres_neworder_indx]
rownames(inf_matrix_raw) <- coati_ids_cut$name[pres_neworder_indx]
colnames(inf_matrix_raw) <- coati_ids_cut$name[pres_neworder_indx]

#saving plot 
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'pres_influence_raw.png')) 
visualize_network_matrix_presedente(inf_matrix_raw, coati_ids_cut[pres_neworder_indx,])
dev.off()

#as the matrix is not symetrical, getting the mean of each dyad (for now - but could do something different later)
inf_matrix_mean <- (inf_matrix_raw + t(inf_matrix_raw)) / 2

#mean scores for each dyad matrix:
visualize_network_matrix_presedente(inf_matrix_mean, coati_ids_cut[pres_neworder_indx,])

set.seed(2)

#MRQAP to test whether influence scores correlate to relatedness
mod <- mrqap.dsp(inf_matrix_mean~pres_gen_matrix, directed="undirected", diagonal = F) 
mod #doesn't look like it with the mean values, but the larger subgroup do have higher influence scores compared to the smaller subgroup (Gillard, )


#saving plot
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'pres_influence_meanscores.png')) 
visualize_network_matrix_presedente(inf_matrix_mean, coati_ids_cut[pres_neworder_indx,]) 
dev.off() 

