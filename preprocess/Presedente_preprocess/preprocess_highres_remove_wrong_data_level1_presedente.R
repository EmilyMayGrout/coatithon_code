#FOR PRESEDENTE GROUP

#this script is removing the times in the highres period where the collar was not on the animal

#remove Peron times when collar was off

#output of this code will make a new gps_file called "presedente_xy_10min_level1.RData"
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/'
gps_file <- "presedente_xy_highres_level0.RData" 
id_file <- 'coati_ids.RData'

outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#Peron collar fell off on 21.01.23 at 17:16 UTC
#Peron recaptured on 25.01.23 and released at 22:40 UTC

#ts[32401] is at 22.01.23 at 11:00 (remove data from here as the high res finishes at 14:00 UTC)
# ts[75601] is at 26.01.23 at 11:00

#add NA's for Peron
xs[17,32401:75601] <- NA
ys[17,32401:75601] <- NA


save(list=c('xs','ys','ts'), file = paste0(outdir,'presedente_xy_highres_level1.RData'))
