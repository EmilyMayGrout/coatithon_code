#this script is for the MRQAP analysis for Trago group  


library(sna)   
library(asnipe)  
library(abind)  

#--------PARAMS------- 
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/" 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/trago_results/' 
gps_file <- "trago_xy_10min_level0.RData" 
id_file <- 'trago_coati_ids.RData' 

trago_gps_matrix <- 'trago_matrix_10min_proptimeinsamesubgroup.txt' #saved from plot3 in fission_fusion_galaxy  
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)  


#make age sex matrices in excel 



#read in library of functions  
setwd(code_dir)  
source('coati_function_library.R')  

#load data  
setwd(data_dir)  
load(id_file)  

trago_matrix <- read.table(trago_gps_matrix, header = T)  
#make subgroup membership matrix a matrix  
trago_matrix <- as.matrix(trago_matrix)  


# G8 - Tequila 
# G9 - Gin 
# G32 - Bailey 
# G33 - Rum 
# G34 - Whisky 
# G35 - Sake 
# G36 - Amarulla 
# G37 - Tiger 
# G38 - Limoncello 
# G39 - Cerveza 


#Amarulla, Limoncello, Whisky, Cerveza, Sake, Tiger, Bailey, Rum, Tequila 
trago_neworder_indx <- c(2,3,7:9,1,5,4,6) 
trago_inds_neworder <- c("G36", "G38", "G34", "G39", "G35", "G37", "G32", "G33", "G8")  

#first column is index numbers but I want this to be the IDs  
rownames(all_matrix) <- all_matrix$X  
#remove first column as this is incorrect for the subsetting to work  
all_matrix <- all_matrix[,-1]  

gen_matrix <- all_matrix[trago_inds_neworder, trago_inds_neworder]  
gen_matrix <- as.matrix(gen_matrix)  

diag(gen_matrix) <- NA  

n_inds <- 9  
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'trago_genetics_level0_noGin.png'))  
visualize_network_matrix_trago(gen_matrix, coati_ids[trago_neworder_indx,])  
dev.off()  





# #adding Gin to coati_ids 
# coati_ids_withGin <- rbind(coati_ids, list('Gin', 'na', 'Adult', 'Female', '#FF0000')) 
#  
# trago_neworder_indx <- c(2,3,7:9,1,5,4,6,10) 
# trago_inds_neworder <- c("G36", "G38", "G34", "G39", "G35", "G37", "G32", "G33", "G8", "G9")  
# gen_matrix <- all_matrix[trago_inds_neworder, trago_inds_neworder]  
# gen_matrix <- as.matrix(gen_matrix)  
#  
# diag(gen_matrix) <- NA  
#  
# n_inds <- 10 
# png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'trago_genetics_withGin_level0.png'))  
# visualize_network_matrix_trago(gen_matrix, coati_ids_withGin[trago_neworder_indx,])  
# dev.off()  
#  





















