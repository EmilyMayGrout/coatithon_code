#this script is for the MRQAP analysis for Galaxy group 

#load in libraries
library(sna)  
library(asnipe) 
library(abind) 
library(fields) 
library(ecodist)
library(vegan)

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/" 
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/' 
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/' 
id_file <- 'galaxy_coati_ids.RData' 
gal_gps_matrix <- 'gal_matrix_10min_proptimeinsamesubgroup_50m.txt' #saved from plot3a in fission_fusion_galaxy 

#read in the genetics matrix for all groups
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T) 

#I manually made these matrices in excel - of the dyad are the same age/sex, they get 1, if different, they get 0
#load from GoogleDrive
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/galaxy/gal_sex_matrix.csv', header = T) 
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/galaxy/gal_age_matrix.csv', header = T) 

#read in library of functions 
setwd(code_dir) 
source('coati_function_library_V1.R') 

#load data 
setwd(data_dir) 
load(id_file) 

gal_matrix <- read.table(gal_gps_matrix, header = T) 
#make subgroup membership matrix a matrix 
gal_matrix <- as.matrix(gal_matrix) 


#order of group in gal_gps matrix to the coati_ids order: c(5,1,11,4,10,2,3,6,7,8,9) 
#this order is gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto 

#Galaxy ID codes to subset correct individuals from genetics matrix 
# G1 -  Sol 
# G15B/C -  Estrella 
# G3 -  Luna 
# G10 - Gus 
# G11 - Saturno 
# G12 - Venus 
# G13 - Orbita 
# G14 - Pluto 
# G15C/B - Cometa 
# G16 - Lucero 
# G17 - Planeta 
# G18 - Quasar 

#order for gal_matrix is: gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto  

gal_inds_neworder <- c("G10", "G18", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14") 
gal_neworder_indx <- c(4,1,10,3,9,2,5,6,7,8) 

#first column is index numbers but I want this to be the IDs 
rownames(all_matrix) <- all_matrix$X 
#remove first column as this is incorrect for the subsetting to work 
all_matrix <- all_matrix[,-1] 

#gen_matrix <- all_matrix[gal_inds, gal_inds] 
gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder] 
gen_matrix <- as.matrix(gen_matrix) 


#----------------------------------------------------------------- 
#adding Sol into matrix (the only adult female we didn't collar because she avoided the traps but we caught her in the full group recapture to remove their collars)
#gus, quasar, sol, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto  
n_inds <- 12
gal_inds_neworder <- c("G10", "G18", "G1", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14") 

gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder] 
gen_matrix <- as.matrix(gen_matrix) 
#adding Sol to coati_ids
coati_ids_withSol <- rbind(coati_ids, list('Sol', 'na', 'Adult', 'Female', '#FF0000'))

gal_neworder_indx <-  c(5,1,12,11,4,10,2,3,6,7,8,9) 
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'gen_cometaG15C_estrellaG15B_Sol_level1.png')) 
visualize_network_matrix_galaxy(gen_matrix, coati_ids_withSol[gal_neworder_indx,]) 
dev.off() 

#---------------------------------------------------------------- 
#genetics matrix for all inds: 
# Cometa is 15B 
# Estrella is 15C 

n_inds <- 11 
#make order for gen_matrix the same as gal_matrix: 
#gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto  
gal_inds_neworder <- c("G10", "G18", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14") #same as above

gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder] 
gen_matrix <- as.matrix(gen_matrix) 
diag(gen_matrix) <- NA 
gal_neworder_indx <-  c(5,1,11,4,10,2,3,6,7,8,9) 

#Figure 3b
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir, 'gen_cometaG15C_estrellaG15B_level1.png'))
visualize_network_matrix_galaxy(gen_matrix, coati_ids[gal_neworder_indx,]) 
dev.off()

#make age and sex matrices actual matrices: 
# order of sex/age matrix (same as coati_ids order): Quasar,Estrella,Venus,Lucero,Gus,Orbita,Planeta,Saturno,Pluto,Luna,Cometa 

rownames(age_matrix) <- age_matrix$X 
age_matrix <- age_matrix[,-1] 
age_matrix <- as.matrix(age_matrix) 
diag(age_matrix) <- NA 
rownames(sex_matrix) <- sex_matrix$X 
sex_matrix <- sex_matrix[,-1] 
sex_matrix <- as.matrix(sex_matrix) 
diag(sex_matrix) <- NA 

#now need to put the age/sex classes into this order (so its the same for gen_matrix and gal_matrix) 
#order: gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto  
gal_neworder_indx <-  c(5,1,11,4,10,2,3,6,7,8,9) 
gal_inds_neworder_agesex <- c("Gus", "Quasar", "Cometa", "Lucero", "Luna", "Estrella","Venus", "Orbita", "Planeta", "Saturno", "Pluto") #same as above  
sex_matrix <- sex_matrix[gal_inds_neworder_agesex, gal_inds_neworder_agesex] 
age_matrix <- age_matrix[gal_inds_neworder_agesex, gal_inds_neworder_agesex] 

# NOW WE HAVE ALL 4 MATRICES WITH THE SAME ORDER OF GROUP MEMBERS!

#png(height = 900, width = 1400, units = 'px', filename = paste0(plot_dir,'all_matrices.png'))
#par(mfrow=c(2,2))
visualize_network_matrix_galaxy(gal_matrix, coati_ids[gal_neworder_indx,])
mtext("1) Proportion of time each dyad was in the same subgroup", cex = 1.2)
visualize_network_matrix_galaxy(age_matrix, coati_ids[gal_neworder_indx,]) 
mtext("2) Age homophily", cex = 1.2)
visualize_network_matrix_galaxy(sex_matrix, coati_ids[gal_neworder_indx,]) 
mtext("3) Sex homophily", cex = 1.2)
visualize_network_matrix_galaxy(gen_matrix, coati_ids[gal_neworder_indx,]) 
mtext("4) Genetics - Triadic Maximum Likelihood method", cex = 1.2)

#dev.off()

#save these matrices to send to Tiffany
save(gal_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/galaxy_subgrouping.Rdata')
save(age_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/galaxy_age.Rdata')
save(sex_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/galaxy_sex.Rdata')
galaxy_gen_matrix <- gen_matrix
save(galaxy_gen_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/galaxy_relatedness.Rdata')
galaxy_ids <- coati_ids[gal_neworder_indx,]
save(galaxy_ids, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/galaxy_ids.Rdata')


#-------------------------------------------------------------------------- 

#MRQAP with Double-Semi-Partialing (DSP) 
set.seed(2)

#table 2:
mod <- mrqap.dsp(gal_matrix~age_matrix+sex_matrix+gen_matrix, directed="undirected", diagonal = F) 
mod
#age and sex don't influence the subgroup membership but genetics does

#write.csv(gal_matrix, file = "C:/Users/egrout/Dropbox/coatithon/processed/genetics/galaxy/gal_matrix.csv")
#write.csv(gen_matrix, file = "C:/Users/egrout/Dropbox/coatithon/processed/genetics/galaxy/gen_matrix.csv")

#Figure S6 - scatter plot for genetics and prop time together
png(height = 1000, width = 1000, units = 'px', filename = paste0(plot_dir,'gen_gal_scatter_2.png'))
par(mar=c(8,8,6,3), mgp=c(5,1.2,0))
gen_vec <- gen_matrix[upper.tri(gen_matrix)]
gal_vec <- gal_matrix[upper.tri(gal_matrix)]
df <- data.frame(cbind(gen_vec, gal_vec))
plot(df$gen_vec, df$gal_vec, ylab = "Proportion of time together", xlab = "Relatedness (Triadic Maximum Likelihood)", pch = 19, col = "darkolivegreen3", cex = 3, cex.lab = 3, cex.axis = 2.5, yaxt = "n", ylim = c(0.3,1))
axis(2, at = c(0,0.3,0.6,0.9), cex.axis = 3, las = 1)
abline(lm(gal_vec ~ gen_vec, data = data.frame(df)), col = "black", lwd = 2)

#str(summary(lm(gal_vec ~ gen_vec, data = data.frame(df)))) #rsquared = 0.182
dev.off()

