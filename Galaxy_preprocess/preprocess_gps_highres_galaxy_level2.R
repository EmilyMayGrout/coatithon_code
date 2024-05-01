#This code is removing errors which were manually found by going through all events detected with the "analyse_ff_event" function in the identify_splits_and_merges code to see which events were caused by GPS error - the time indexes were extracted and these will be replaced with NA's here

group <- 'galaxy' #subdirectory where the group data is stored
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
codedir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'

#SOURCE FUNCTIONS
setwd(codedir)
source('coati_function_library.R')

#LOAD DATA
#read in coati ids
setwd(groupdir)
load(file=paste0(group,'_coati_ids.RData'))

#read in timestamp data - notice the level number, if 2, then its not got GPS speed errors
load(file=paste0(group,'_xy_highres_level1.RData'))

# 1 - "Quasar" 
# 2 - "Estrella" 
# 3 - "Venus"  
# 4 - "Lucero"  
# 5 - "Gus"   
# 6 - "Orbita" 
# 7 - "Planeta" 
# 8 - "Saturno" 
# 9 - "Pluto"  
# 10 - "Luna"  
# 11 - "Cometa" 


#Luna
xs[10,111:120] <- NA
ys[10,111:120] <- NA

xs[10, 53500:53600] <- NA
ys[10, 53500:53600] <- NA

#Orbita
xs[6,180:200] <- NA
ys[6,180:200] <- NA

xs[6,250:350] <- NA
ys[6,250:350] <- NA

#Pluto
xs[9,2120:2140] <- NA
ys[9,2120:2140] <- NA

#Estrella
xs[2,10800:10850] <- NA
ys[2,10800:10850] <- NA

#Cometa
xs[11, 26600:26800] <- NA
ys[11, 26600:26800] <- NA

xs[11,50450:50460] <- NA
ys[11,50450:50460] <- NA

xs[11, 86280:86370] <- NA
ys[11, 86280:86370] <- NA

#Quasar
xs[1,50430:50480] <- NA
ys[1,50430:50480] <- NA

xs[1,59085:59150] <- NA
ys[1,59085:59150] <- NA

#Venus
xs[3,65720:65850] <- NA
ys[3,65720:65850] <- NA

#Saturno
xs[8,163650:163700] <- NA
ys[8,163650:163700] <- NA




save(list=c('xs','ys','ts'), file = paste0(groupdir,'galaxy_xy_highres_level2.RData'))





