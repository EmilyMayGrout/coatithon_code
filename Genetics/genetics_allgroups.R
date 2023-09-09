#this script is bringing the genetics for all coatis into the same matrices


code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/' 
#read in library of functions 
setwd(code_dir) 
source('coati_function_library.R') 

all_matrix <- read.csv('C:/Users/egrout/Dropbox/coatithon/processed/genetics/CoatiTrioMLmatrix.csv', header = T)  
#first column is index numbers but I want this to be the IDs  
rownames(all_matrix) <- all_matrix$X  
#remove first column as this is incorrect for the subsetting to work  
all_matrix <- all_matrix[,-1]  


#load in ID files
load("C:/Users/egrout/Dropbox/coatithon/processed/2022/trago/trago_coati_ids.RData")
#adding Gin
trago_withGin <- rbind(coati_ids, list('Gin', 'na', 'Adult', 'Female', '#FF0000'))

load("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/galaxy_coati_ids.RData")
#adding Sol to coati_ids
galaxy_withSOl <- rbind(coati_ids, list('Sol', 'na', 'Adult', 'Female', '#FF0000'))

load("C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/presedente_coati_ids_level1.RData")
pres_ids <- coati_ids

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

trago_neworder_indx <- c(2,3,7:9,1,5,4,6,10) 
trago_inds_neworder <- c("G36", "G38", "G34", "G39", "G35", "G37", "G32", "G33", "G8", "G9") 
trago_withGin_ordered <- trago_withGin[trago_neworder_indx,]

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

#gus, quasar, sol, cometa, lucero, luna, estrella, venus, orbita, planeta, saturno, pluto  
gal_inds_neworder <- c("G10", "G18", "G1", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14") 
gal_neworder_indx <-  c(5,1,12,11,4,10,2,3,6,7,8,9) 
galaxy_withSOl_ordered <- galaxy_withSOl[gal_neworder_indx,]
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
# B08 - Mandela
# B18 - Lula
# B19 - Kenyatta
# B21 - Gendry
# B20 - Sam
# B22 - Wildflower

#ardern, merkel, moscoso, cleopatra, gandhi, obama, castro, mujica, may, khan, peron, zelenskyy, meir, truss, gillard, torrijos, mandela, lula, kenyatta, gendry, sam, wildflower
pres_inds_neworder <- c("B12.L001", "B17.L001", "B16.L001", "B09.L001", "B03.L001", "B01.L001", "B02.L001", "B04.L001", "B14.L001", "B10.L001", "B07.L001", "B05.L001", "B06.L001", "B13.L001", "B11.L001", "B15.L001", "B08.L001", "B18.L001", "B19.L001", "B21.L001", "B20.L001", "B22.L001") 
pres_neworder_indx <- c(1,13,14,3,4,16,2,15,11,6,17,22,12,20,7,19, 10,9,8,5,18,21)
pres_ids_ordered <- pres_ids[pres_neworder_indx,]


#first need to rbind the coati IDs to make one table for all inds
all_ids <- rbind(trago_withGin_ordered, galaxy_withSOl_ordered, pres_ids_ordered)

#then need to match this order with the genetic codes to subset the data

all_neworder <- c("G36", "G38", "G34", "G39", "G35", "G37", "G32", "G33", "G8", "G9","G10", "G18", "G1", "G15C", "G16", "G3", "G15B","G12", "G13", "G17", "G11", "G14", "B12.L001", "B17.L001", "B16.L001", "B09.L001", "B03.L001", "B01.L001", "B02.L001", "B04.L001", "B14.L001", "B10.L001", "B07.L001", "B05.L001", "B06.L001", "B13.L001", "B11.L001", "B15.L001", "B08.L001", "B18.L001", "B19.L001", "B21.L001", "B20.L001", "B22.L001") 

gen_matrix <- all_matrix[all_neworder, all_neworder]  
gen_matrix <- as.matrix(gen_matrix) 
diag(gen_matrix) <- NA  

#need to change the rep value to -.018 in the visualize_network_matrix_galaxy function for the agesex dots to be in correct place
n_inds <- 44 
png(height = 1200, width = 1300, units = 'px', filename = 'C:/Users/egrout/Dropbox/coatithon/results/all_genetics.png')  
visualize_network_matrix_galaxy(gen_matrix, all_ids)  
dev.off()  


