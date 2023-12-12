#this script is to make the figure for Trago groups relatedness

#load in libraries
library(sna)   
library(asnipe)  
library(abind)  

#--------PARAMS------- 
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/" 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/trago_results/' 
gps_file <- "trago_xy_10min_level0.RData" 
id_file <- 'trago_coati_ids.RData' 

trago_gps_matrix <- 'trago_matrix_10min_proptimeinsamesubgroup.txt' #saved from plot3/FigureS1 in fission_fusion_trago 

#read in the genetics matrix for all groups
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)  

#I manually made these matrices in excel - of the dyad are the same age/sex, they get 1, if different, they get 0
#load these matrices in from the GoogleDrive
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/trago/trago_sex_matrix.csv', header = T) 
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/trago/trago_age_matrix.csv', header = T) 

rownames(age_matrix) <- age_matrix$X 
age_matrix <- age_matrix[,-1] 
age_matrix <- as.matrix(age_matrix) 
diag(age_matrix) <- NA
#remove Rum and Tequila rows and columns
age_matrix <- age_matrix[-c(8,9), -c(8,9)]

rownames(sex_matrix) <- sex_matrix$X 
sex_matrix <- sex_matrix[,-1] 
sex_matrix <- as.matrix(sex_matrix) 
diag(sex_matrix) <- NA 
sex_matrix <- sex_matrix[-c(8,9), -c(8,9)]

#read in library of functions  
setwd(code_dir)  
source('coati_function_library_V1.R')  

#load data  
setwd(data_dir)  
load(id_file)  

trago_matrix <- read.table(trago_gps_matrix, header = T)  
#make subgroup membership matrix a matrix  
trago_matrix <- as.matrix(trago_matrix)  

# G33 - Rum 
# G34 - Whisky 
# G35 - Sake 
# G36 - Amarulla 
# G37 - Tiger 
# G38 - Limoncello 
# G39 - Cerveza 

#Bailey, Tiger, Amarulla, Limoncello, Whisky, Cerveza, Sake
trago_neworder_indx <- c(4,1:3, 5:7) 
trago_inds_neworder <- c("G32","G37","G36", "G38", "G34", "G39", "G35") 

#first column is index numbers but I want this to be the IDs  
rownames(all_matrix) <- all_matrix$X  
#remove first column as this is incorrect for the subsetting to work  
all_matrix <- all_matrix[,-1]  

gen_matrix <- all_matrix[trago_inds_neworder, trago_inds_neworder]  
gen_matrix <- as.matrix(gen_matrix)  

diag(gen_matrix) <- NA  

#plot S3
n_inds <- 7 
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir ,'trago_genetics_level0_noGin.png'))  
par(mfrow=c(1,1), mar = c(6,6,1,1)) #bottom, left, top, and right
visualize_network_matrix_trago(gen_matrix, coati_ids[trago_neworder_indx,])  
dev.off()  

