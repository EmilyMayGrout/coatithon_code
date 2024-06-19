#this script is getting the proportion of time Trago group members are together, and removing times when Tequila and Rum leave the group

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/trago_results/'
gps_file <- "trago_xy_10min_level0_all.RData"
id_file <- 'trago_coati_ids_all.RData'

#I manually made these matrices in excel - of the dyad are the same age/sex, they get 1, if different, they get 0
#load these matrices in from the GoogleDrive
sex_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/trago/trago_sex_matrix.csv', header = T) 
age_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/trago/trago_age_matrix.csv', header = T) 

#read in the genetics matrix for all groups
all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)  


#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(sna)   
library(asnipe)  
library(abind)  

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


subgroup_data <- get_subgroup_data(xs, ys, R=50)

ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    
    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

diag(ff_net) <- NA
new_order <- c(2,3,7:9,1,5,4,6)
ffnet_reorder <- ff_net[new_order, new_order]

visualize_network_matrix_trago(ffnet_reorder, coati_ids[new_order,])

#get timestamp when Rum and Tequila not in group

#Tequila 
#ts[180:1399]
#rum
#ts[320:1399]

xs[6, 180:1399] <- NA
ys[6, 180:1399] <- NA

xs[4, 320:1399] <- NA
ys[4, 320:1399] <- NA

subgroup_data <- get_subgroup_data(xs, ys, R=50)
ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)
for(i in 1:n_inds){
  for(j in 1:n_inds){
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}
diag(ff_net) <- NA
new_order <- c(1,2,3,7,8,9,4,5,6)
ffnet_reorder <- ff_net[new_order, new_order]
visualize_network_matrix_trago(ffnet_reorder, coati_ids[new_order,])

#--------------------------------------------------------------------


rownames(age_matrix) <- age_matrix$X 
age_matrix <- age_matrix[,-1] 
age_matrix <- as.matrix(age_matrix) 
diag(age_matrix) <- NA

rownames(sex_matrix) <- sex_matrix$X 
sex_matrix <- sex_matrix[,-1] 
sex_matrix <- as.matrix(sex_matrix) 
diag(sex_matrix) <- NA 
 
trago_matrix <- ffnet_reorder

# G8 - Tequila 
# G32 - Bailey 
# G33 - Rum 
# G34 - Whisky 
# G35 - Sake 
# G36 - Amarulla 
# G37 - Tiger 
# G38 - Limoncello 
# G39 - Cerveza 

#Amarulla, Limoncello, Whisky, Cerveza, Sake, Tiger, Bailey, Rum, Tequila
trago_neworder_indx <- c(2,3,7,8,9,1,5,4,6)
trago_inds_neworder <- c("G36", "G38","G34","G39", "G35","G37","G32","G33","G8") 


#first column is index numbers but I want this to be the IDs  
rownames(all_matrix) <- all_matrix$X  
#remove first column as this is incorrect for the subsetting to work  
all_matrix <- all_matrix[,-1]  

gen_matrix <- all_matrix[trago_inds_neworder, trago_inds_neworder]  
gen_matrix <- as.matrix(gen_matrix)  

diag(gen_matrix) <- NA  

#plot S3
n_inds <- 9
par(mfrow=c(1,1), mar = c(6,6,1,1)) #bottom, left, top, and right
visualize_network_matrix_trago(gen_matrix, coati_ids[trago_neworder_indx,])  
visualize_network_matrix_trago(sex_matrix, coati_ids[trago_neworder_indx,])  
visualize_network_matrix_trago(age_matrix, coati_ids[trago_neworder_indx,])  
visualize_network_matrix_trago(trago_matrix, coati_ids[trago_neworder_indx,])  

#save these matrices to send to Tiffany
save(trago_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/trago_subgrouping.Rdata')
save(age_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/trago_age.Rdata')
save(sex_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/trago_sex.Rdata')
trago_gen_matrix <- gen_matrix
save(trago_gen_matrix, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/trago_relatedness.Rdata')
trago_ids <- coati_ids[trago_neworder_indx,]
save(trago_ids, file = 'C:/Users/egrout/Dropbox/coatithon/processed/microbiome/trago_ids.Rdata')



