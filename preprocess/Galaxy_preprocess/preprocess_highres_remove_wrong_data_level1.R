#FOR GALAXY GROUP

#Removing time when Venus was in the tree stationary
#She likely found the group after the 1Hz period on the 28.12.21. So we should add NA's from the start of behaving weird (26.12.21) to 29.12.21 11:00 (when she was moving with the group)


#output of this code will make a new gps_file called "presedente_xy_highres_level1.RData"
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/'
gps_file <- "galaxy_xy_highres_level0.RData" 
id_file <- 'coati_ids.RData'

outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#ts[54001] is 2021-12-29 11:00:00

#need to keep the 24th and 25th, as Venus is with the group then
#ts[21601] is  "2021-12-26 11:00:00 UTC"

#so going to add NA's for Venus from ts[21601] to ts[54001] 
xs[3,21601:54001] <- NA
ys[3,21601:54001] <- NA


save(list=c('xs','ys','ts'), file = paste0(outdir,'galaxy_xy_highres_level1.RData'))



