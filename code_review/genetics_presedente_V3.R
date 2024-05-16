#this script is for the MRQAP analysis for Presidente group 

#load in libraries
library(sna)  
library(asnipe) 
library(abind) 
library(fields) 
library(ecodist)
library(vegan)

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/" 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/' 
id_file <- 'presedente_coati_ids.RData'  #level one is where May and Cleopatra have been changed to subadults
pres_gps_matrix <- 'presedente_matrix_10min_proptimeinsamesubgroup_50m.txt' #saved from plot3c in fission_fusion_presedente

#read in the genetics matrix for all groups
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T) 

#I manually made these matrices in excel - of the dyad are the same age/sex, they get 1, if different, they get 0
#load from GoogleDrive 
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/presedente/pres_sex_matrix.csv', header = T) 
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/presedente/pres_age_matrix.csv', header = T) 

#read in library of functions 
setwd(code_dir) 
source('coati_function_library_V1.R') 

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

#ardern, merkel, moscoso, cleopatra, gandhi, obama, castro, mujica, may, khan, peron, zelenskyy, meir, truss, gillard, torrijos 
pres_inds_neworder <- c("B12.L001", "B17.L001", "B16.L001", "B09.L001", "B03.L001", "B01.L001", "B02.L001", "B04.L001", "B14.L001", "B10.L001", "B07.L001", "B05.L001", "B06.L001", "B13.L001", "B11.L001", "B15.L001") 

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

#figure 3d
n_inds <- 16 
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'pres_genetics.png')) 
visualize_network_matrix_presedente(gen_matrix, coati_ids_cut[pres_neworder_indx,]) 
dev.off() 


#make age and sex matrices actual matrices: 
# order of sex/age matrix (same as coati_ids order): ardern, merkel, moscoso, cleopatra, gandhi, obama, castro, mujica, may, khan, peron, zelenskyy, meir, truss, gillard, torrijos 

rownames(age_matrix) <- age_matrix$X 
age_matrix <- age_matrix[,-1] 
age_matrix <- as.matrix(age_matrix) 
diag(age_matrix) <- NA 
rownames(sex_matrix) <- sex_matrix$X 
sex_matrix <- sex_matrix[,-1] 
sex_matrix <- as.matrix(sex_matrix) 
diag(sex_matrix) <- NA 


setwd(plot_dir)


#rename Ardern in coati_ids_cut
coati_ids_cut$name[coati_ids_cut$name == "Ardera"] <- "Ardern"

#png(height = 900, width = 1400, units = 'px', filename = paste0(plot_dir,'all_matrices.png'))
par(mfrow=c(2,2))
visualize_network_matrix_presedente(pres_matrix, coati_ids_cut[pres_neworder_indx,])
mtext("1) Proportion of time each dyad was in the same subgroup", cex = 1.2)
visualize_network_matrix_presedente(age_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("2) Age homophily", cex = 1.2)
visualize_network_matrix_presedente(sex_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("3) Sex homophily", cex = 1.2)
visualize_network_matrix_presedente(gen_matrix, coati_ids_cut[pres_neworder_indx,]) 
mtext("4) Genetics - Triadic Maximum Likelihood method", cex = 1.2)

#dev.off()

#save these matrices to send to Tiffany
save(pres_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_subgrouping.Rdata')
save(age_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_age.Rdata')
save(sex_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_sex.Rdata')
pres_gen_matrix <- gen_matrix
save(pres_gen_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_relatedness.Rdata')
pres_ids <- coati_ids_cut[pres_neworder_indx,]
save(pres_ids, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/pres_ids.Rdata')



#MRQAP with Double-Semi-Partialing (DSP) 
set.seed(2)

#table 2
mod <- mrqap.dsp(pres_matrix~age_matrix+sex_matrix+gen_matrix, directed="undirected", diagonal = F) 


#-----------------------------------------------------------------------
#reordering the matrices to age/sex class to visualise th differences more clearly, this plot will not go in the paper as the order of the subgrouping matrix would not be as clear

#ardern, merkel, truss, moscoso, cleopatra, may, peron, gillard, torrijos, khan, gandhi, obama, mujica, castro, zelenskyy, meir
pres_inds_neworder_agesex <- c("B12.L001", "B17.L001", "B13.L001", "B16.L001", "B09.L001", "B14.L001",  "B07.L001", "B11.L001", "B15.L001","B10.L001", "B03.L001", "B01.L001","B04.L001",  "B02.L001",  "B05.L001", "B06.L001") 
pres_neworder_indx_agesex <- c(1,9,15,10,3,7,13,6,14,5,4,12,11,2,16,8)
gen_matrix_agesex <- all_matrix[pres_inds_neworder_agesex, pres_inds_neworder_agesex] 
gen_matrix_agesex  <- as.matrix(gen_matrix_agesex) 
diag(gen_matrix_agesex) <- NA 
n_inds <- 16 
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'pres_genetics_agesexorder.png')) 
visualize_network_matrix_presedente(gen_matrix_agesex, coati_ids_cut[pres_neworder_indx_agesex,]) 
dev.off() 
#-----------------------------------------------------------------------

#make plot with dyadic strength against relatedness

#figure S6
png(height = 1000, width = 1000, units = 'px', filename = paste0(plot_dir,'gen_pres_scatter_2.png'))
par(mar=c(8,8,6,3), mgp=c(5,1.2,0))
gen_vec <- gen_matrix[upper.tri(gen_matrix)]
pres_vec <- pres_matrix[upper.tri(pres_matrix)]
df <- data.frame(cbind(gen_vec, pres_vec))
plot(df$gen_vec, df$pres_vec, ylab = "Proportion of time together", xlab = "Relatedness (Triadic Maximum Likelihood)", pch = 19, col = "aquamarine3", cex = 3, cex.lab = 3, cex.axis = 2.5, yaxt = "n", ylim = c(0.3, 1))
axis(2, at = c(0,0.3,0.6,0.9), cex.axis = 3, las = 1)
abline(lm(pres_vec ~ gen_vec, data = data.frame(df)), col = "black", lwd = 2)
text(x = 0.4, y = 0.76,labels = expression(paste("r"^2, " = 0.089")), cex = 3)
text(x = 0.045, y = 0.99, labels = "(b)", cex = 3)
mtext("Presidente group", side = 3, line = 2, cex = 3)




#str(summary(lm(pres_vec ~ gen_vec, data = data.frame(df)))) #rsquared 0.0885

dev.off()






