#this script is for the MRQAP analysis for Presidente group 

#TODO go through the age classes and check which ones need changing - make these changes for the coati_ids and the age matrices and rerun the analysis - but removing within group associations in the model and add it as a mantel test with subgorup membership. Then get the 3m matrices and run this mantel test
#then need to make the age class changes to the full group associations

library(sna)  
library(asnipe) 
library(abind) 
library(fields) 
library(ecodist)
library(vegan)

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/" 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/' 
id_file <- 'presedente_coati_ids_level1.RData'  #level one is where May and Cleopatra have been changed to subadults
pres_gps_matrix <- 'pres_matrix_10min_proptimeinsamesubgroup.txt' #saved from plot3 in fission_fusion_presedente
pres_gps_matrix_full <- 'presedente_matrix_10min_proptimeinfullgroup.txt' #saved from plot4 in fission_fusion_presedente
pres_gps_matrix_full3m <- 'presedente_matrix_10min_proptimeinfullgroup3m.txt' #saved from plot4 in 

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

#changing the juvenile colour to subadult colour for Cleopatra and May
#coati_ids[3,5] <- "#FFAA66"
#coati_ids[11,5] <- "#FFAA66"
#save(coati_ids, file="presedente_coati_ids_level1.RData")


pres_matrix <- read.table(pres_gps_matrix, header = T) 
#make subgroup membership matrix a matrix 
pres_matrix <- as.matrix(pres_matrix) 
diag(pres_matrix) <- NA 

pres_matrix_full <- read.table(pres_gps_matrix_full, header = T) 
#make subgroup membership matrix a matrix 
pres_matrix_full <- as.matrix(pres_matrix_full) 
diag(pres_matrix_full) <- NA 

pres_matrix_full3m <- read.table(pres_gps_matrix_full3m, header = T) 
#make subgroup membership matrix a matrix 
pres_matrix_full3m <- as.matrix(pres_matrix_full3m) 
diag(pres_matrix_full3m) <- NA 



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
visualize_network_matrix_presedente(gen_matrix, coati_ids_cut[pres_neworder_indx,]) 
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


setwd(plot_dir)

#png(height = 900, width = 1400, units = 'px', filename = paste0(plot_dir,'all_matrices.png'))
par(mfrow=c(2,3))
visualize_network_matrix_presedente(pres_matrix, coati_ids_cut[pres_neworder_indx,])
mtext("1) Proportion of time each dyad was in the same subgroup", cex = 1.2)
visualize_network_matrix_presedente(age_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("2) Age homophily", cex = 1.2)
visualize_network_matrix_presedente(sex_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("3) Sex homophily", cex = 1.2)
visualize_network_matrix_presedente(gen_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("4) Genetics - Triadic Maximum Likelihood method", cex = 1.2)
visualize_network_matrix_presedente(pres_matrix_full, coati_ids_cut[pres_neworder_indx,]) 
mtext('5) Within full group proportion of time within 10m', cex = 1.2)
visualize_network_matrix_presedente(pres_matrix_full3m, coati_ids_cut[pres_neworder_indx,]) 
mtext('6) Within full group proportion of time within 3m', cex = 1.2)

#dev.off()

#MRQAP with Double-Semi-Partialing (DSP) 
#trying out different interactions, but want to see whether age/sex/genetics influence the subgroup membership patterns 
t1 <- mrqap.dsp(pres_matrix~gen_matrix+sex_matrix, directed="undirected") 
t2 <- mrqap.dsp(pres_matrix~gen_matrix+age_matrix, directed="undirected") 

t3 <- mrqap.dsp(pres_matrix~age_matrix+sex_matrix+gen_matrix, directed="undirected", diagonal = F) 

t4 <- mrqap.dsp(pres_matrix~age_matrix+sex_matrix+gen_matrix+pres_matrix_full, directed="undirected", diagonal = F) 
t4.1 <- mrqap.dsp(pres_matrix~age_matrix+sex_matrix+gen_matrix+pres_matrix_full3m, directed="undirected", diagonal = F) 
#within group associations and genetics explains subgroup membership

#is within group associations affected by genetics 
t5 <- mrqap.dsp(pres_matrix_full ~ gen_matrix+age_matrix+sex_matrix, directed="undirected", diagonal = F)
#genetics, sex and age don't affect within group associations
#also running this test with the 3m within group associations 
t5.1 <- mrqap.dsp(pres_matrix_full ~ gen_matrix+age_matrix+sex_matrix, directed="undirected", diagonal = F)

#Ben comment: distance to individuals shouldn't be in the model as it subgroup membership could be dependent on it, so use results for t3 and then do a mantel test to look at the distance relationships

t6 <- mantel(xdis = pres_matrix, ydis = pres_matrix_full, method = "pearson", permutations = 999)
t7 <- mantel(xdis = pres_matrix, ydis = pres_matrix_full3m, method = "pearson", permutations = 999)

























#-----------------------------------------------------------------------