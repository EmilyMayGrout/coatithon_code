#this script is for the MRQAP analysis for Galaxy group

library(sna) 
library(asnipe)
library(abind)

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
id_file <- 'galaxy_coati_ids.RData'
gal_gps_matrix <- 'gal_matrix_10min_proptimeinsamesubgroup.txt'
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)

#I manually made these matrices in excel
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/sex_matrix.csv', header = T)
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/age_matrix.csv', header = T)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

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
# G2 -  Estrella
# G3 -  Luna
# G10 - Gus
# G11 - Saturno
# G12 - Venus
# G13 - Orbita
# G14 - Pluto
# G15 - Cometa
# G16 - Lucero
# G17 - Planeta
# G18 - Quasar

#order for gal_matrix is: gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto 
#G10, G18, G15A, G16, G3, G2, G12, G13, G17, G11, G14

gal_inds <- c("G1", "G3", "G10", "G11", "G12", "G13", "G14", "G15A", "G16", "G17", "G18")
gal_inds_neworder <- c("G10", "G18", "G15A", "G16", "G3", "G12", "G13", "G17", "G11", "G14")
gal_neworder_indx <- c(4,1,10,3,9,2,5,6,7,8)

#missing G2 in genetics
#why is there a G15A, G15B, G15C? which should I use?

#first column is index numbers but I want this to be the IDs
rownames(all_matrix) <- all_matrix$X
#remove first column as this is incorrect for the subsetting to work
all_matrix <- all_matrix[,-1]

#gen_matrix <- all_matrix[gal_inds, gal_inds]
gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder]
gen_matrix <- as.matrix(gen_matrix)

#removing Estrella from coati_ids so I can visualise the relatedness
coati_ids_cut <- coati_ids[-2,]
diag(gen_matrix) <- NA

n_inds <- 10
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'gen_cometa_G15A_level1.png'))
visualize_network_matrix(gen_matrix, coati_ids_cut[gal_neworder_indx,])
dev.off()


#-----------------------------------------------------------------
#trying Estrella as G15A and Cometa as G15B
#gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto 
n_inds <- 11
gal_inds_neworder <- c("G10", "G18", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14")
gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder]
gen_matrix <- as.matrix(gen_matrix)
gal_neworder_indx <-  c(5,1,11,4,10,2,3,6,7,8,9)
png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'gen_cometaG15C_estrellaG15B_level1.png'))
visualize_network_matrix(gen_matrix, coati_ids[gal_neworder_indx,])
dev.off()

#----------------------------------------------------------------
#for now will do the genetics matrix for all inds:
# Cometa will be 15B
# Estrella will be 15C

n_inds <- 11
#make order for gen_matrix the same as gal_matrix:
#gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto 
gal_inds_neworder <- c("G10", "G18", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14") #same as above 

gen_matrix <- all_matrix[gal_inds_neworder, gal_inds_neworder]
gen_matrix <- as.matrix(gen_matrix)
diag(gen_matrix) <- NA
gal_neworder_indx <-  c(5,1,11,4,10,2,3,6,7,8,9)

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
#gus, quasar, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto 
gal_neworder_indx <-  c(5,1,11,4,10,2,3,6,7,8,9)
gal_inds_neworder_agesex <- c("Gus", "Quasar", "Cometa", "Lucero", "Luna", "Estrella","Venus", "Orbita", "Planeta", "Saturno", "Pluto") #same as above 
sex_matrix <- sex_matrix[gal_inds_neworder_agesex, gal_inds_neworder_agesex]
age_matrix <- age_matrix[gal_inds_neworder_agesex, gal_inds_neworder_agesex]


# NOW WE HAVE ALL 4 MATRICES WITH THE SAME ORDER OF GROUP MEMBERS!!
#  only thing that needs to be checked is the genetics matrix (due to not having Estrella and having 3x Cometa genetic values G15A/B/C)

#can check the arrays to make sure they're all in the right order:
#visualize_network_matrix(gal_matrix, coati_ids[gal_neworder_indx,])
#visualize_network_matrix(age_matrix, coati_ids[gal_neworder_indx,])
#visualize_network_matrix(sex_matrix, coati_ids[gal_neworder_indx,])
#visualize_network_matrix(gen_matrix, coati_ids[gal_neworder_indx,])

#--------------------------------------------------------------------------

#put all matrices into an array
#array <- abind(gal_matrix, gen_matrix, sex_matrix, age_matrix, along=3)

#MRQAP with Double-Semi-Partialing (DSP)
#trying out different interactions, but want to see whether age/sex/genetics influence the subgroup membership patterns
t1 <- mrqap.dsp(gal_matrix~gen_matrix+sex_matrix, directed="undirected")
t2 <- mrqap.dsp(gal_matrix~gen_matrix+age_matrix, directed="undirected")

t3 <- mrqap.dsp(gal_matrix~age_matrix+sex_matrix+gen_matrix, directed="undirected", diagonal = F)
#not sure how to interpret the results, but from what I can understand age and sex don't influence the subgroup membership but genetics does (if these genetics are correct which is unlikely)


#-----------------------------------------------------------------

