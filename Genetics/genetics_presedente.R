#this script is for the MRQAP analysis for Presidente group


library(sna) 
library(asnipe)
library(abind)

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
id_file <- 'presedente_coati_ids.RData'
pres_gps_matrix <- 'presedente_matrix_10min_proptimeinsamesubgroup.txt' #saved from plot3 in fission_fusion_galaxy
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)

#I manually made these matrices in excel
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/presedente/sex_matrix.csv', header = T)
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/presedente/age_matrix.csv', header = T)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(id_file)

pres_matrix <- read.table(pres_gps_matrix, header = T)
#make subgroup membership matrix a matrix
pres_matrix <- as.matrix(pres_matrix)
diag(pres_matrix) <- NA

#order - without Wildflower AND males: 
#new_order <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14)
#ardern, merkel, moscoso, cleopatra, gandhi, obama, castro, mujica, may, khan, peron, zelenskyy, meir, truss, gillard, torrijos


#Presidente ID codes to subset correct individuals from genetics matrix
# B12 - Ardern
# B17 - Merkel
# B16 - Moscoso
# B09 - Cleopatra
# B03 - Gandhi
# B01 - Obama
# B02 - Castro
# B04 - Mujica
# B14 - May
# B10 - Khan
# B07 - Peron
# B05 - Zelenskyy
# B06 - Meir
# B13 - Truss
# B11 - Gillard
# B15 - Torrijos

pres_inds_neworder <- c("B12.L001", "B17.L001", "B16.L001", "B09.L001", "B03.L001", "B01.L001", "B02.L001", "B04.L001", "B14.L001", "B10.L001", "B07.L001", "B05.L001", "B06.L001", "B13.L001", "B11.L001", "B15.L001")

t <- coati_ids[-c(5, 8, 9, 10, 18, 21),]
t$n <- 1:nrow(t)
 
pres_neworder_indx <- c(1,9,10,3,4,12,2,11,7,5,13,16,8,15,6,14) #order for coati_ids when adult males are not in the group

#first column is index numbers but I want this to be the IDs
rownames(all_matrix) <- all_matrix$X
#remove first column as this is incorrect for the subsetting to work
all_matrix <- all_matrix[,-1]

gen_matrix <- all_matrix[pres_inds_neworder, pres_inds_neworder]
gen_matrix <- as.matrix(gen_matrix)
diag(gen_matrix) <- NA

#removing males and Wildflower from coati_ids so I can visualise the relatedness
coati_ids_cut <- coati_ids[-c(5, 8, 9, 10, 18, 21),]


n_inds <- 16
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'pres_genetics.png'))
visualize_network_matrix(gen_matrix, coati_ids_cut[pres_neworder_indx,])
dev.off()


#make age and sex matrices actual matrices:
# order of sex/age matrix (same as coati_ids order): ardern, merkel, moscoso, cleopatra, gandhi, obama, castro, mujica, may, khan, peron, zelenskyy, meir, truss, gillard, torrijos
#pres_neworder_indx <- c(1,13,14,3,4,16,2,15,11,6,17,22,12,20,7,19) #order when males in the coati_ids

rownames(age_matrix) <- age_matrix$X
age_matrix <- age_matrix[,-1]
age_matrix <- as.matrix(age_matrix)
diag(age_matrix) <- NA
rownames(sex_matrix) <- sex_matrix$X
sex_matrix <- sex_matrix[,-1]
sex_matrix <- as.matrix(sex_matrix)
diag(sex_matrix) <- NA


visualize_network_matrix(pres_matrix, coati_ids_cut[pres_neworder_indx,])
visualize_network_matrix(age_matrix, coati_ids_cut[pres_neworder_indx,])
visualize_network_matrix(sex_matrix, coati_ids_cut[pres_neworder_indx,])
visualize_network_matrix(gen_matrix, coati_ids_cut[pres_neworder_indx,])

#MRQAP with Double-Semi-Partialing (DSP)
#trying out different interactions, but want to see whether age/sex/genetics influence the subgroup membership patterns
t1 <- mrqap.dsp(pres_matrix~gen_matrix+sex_matrix, directed="undirected")
t2 <- mrqap.dsp(pres_matrix~gen_matrix+age_matrix, directed="undirected")

t3 <- mrqap.dsp(pres_matrix~age_matrix+sex_matrix+gen_matrix, directed="undirected", diagonal = F)









