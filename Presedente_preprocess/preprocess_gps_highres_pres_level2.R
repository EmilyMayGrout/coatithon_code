#This code is removing errors which were manually found by going through all events detected with the "analyse_ff_event" function in the identify_splits_and_merges code to see which events were caused by GPS error - the time indexes were extracted and these will be replaced with NA's here

group <- 'presedente' #subdirectory where the group data is stored
groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
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

# 1 - "Ardera"    
# 2 - "Castro"     
# 3 - "Cleopatra"  
# 4 - "Gandhi"     
# 5 - "Gendry"     
# 6 - "Khan"       
# 7 - "Gillard"   
# 8 - "Kenyatta"   
# 9 - "Lula"       
# 10 - "Mandela"   
# 11 - "May"       
# 12 - "Meir"      
# 13 - "Merkel"     
# 14 - "Moscoso"   
# 15 - "Mujica"     
# 16 - "Obama"     
# 17 - "Peron"    
# 18 - "Sam"      
# 19 - "Torrijos"  
# 20 - "Truss"    
# 21 - "Wildflower"
# 22 - "Zelenskyy" 

#get the time durations where we want to replace with NA's - these are in the google sheet called Finding correct ff from level1

#Khan
xs[6,1200:1500] <- NA
ys[6,1200:1500] <- NA

xs[6,62500:62600] <- NA
ys[6,62500:62600] <- NA

#Ardern
xs[1,12250:12350] <- NA
ys[1,12250:12350] <- NA

xs[1,76210:76280] <- NA
ys[1,76210:76280] <- NA

#Torrijos
xs[19,17200:17290] <- NA
ys[19,17200:17290] <- NA

#Kenyatta
xs[8,46400:46560] <- NA
ys[8,46400:46560] <- NA

#Merkel
xs[13,59250:59300] <- NA
ys[13,59250:59300] <- NA

xs[13,72050:72080] <- NA
ys[13,72050:72080] <- NA

#Mandela
xs[10,67700:67800] <- NA
ys[10,67700:67800] <- NA

xs[10,116980:117480] <- NA
ys[10,116980:117480] <- NA

#Meir
xs[12,120640:120750] <- NA
ys[12,120640:120750] <- NA

#Gillard
xs[7, 130820:130950] <- NA
ys[7, 130820:130950] <- NA


save(list=c('xs','ys','ts'), file = paste0(groupdir,'presedente_xy_highres_level2.RData'))

