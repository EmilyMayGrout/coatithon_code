# this code is to use the output matrix from the preprocess_gps_lowres and remove the data from times when either:
#1. the ind was not in the group because of capture/handling effect (e.g. Venus)
#2. the ind's collar came off
#3. the gps points were goofy (e.g. Gillard on 29.01.23 from Presedente going in a straight line)

#output of this code will make a new gps_file called "presedente_xy_10min_level1.RData"
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/'
gps_file <- "presedente_xy_10min_level0.RData" 
id_file <- 'coati_ids.RData'

outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"

#load data
setwd(data_dir)
load(gps_file)
load(id_file)


#Peron collar fell off on 21.01.23 at 17:16 UTC
#Peron recaptured on 25.01.23 and released at 22:40 UTC

#double checking that row 17 is Peron, which it is!
plot(xs[2,200:400], ys[2,200:400], type = "l", col = "blue")
lines(xs[17,200:400], ys[17,200:400], type = "l", col = "red")
lines(xs[3,200:400], ys[3,200:400], type = "l", col = "black")

#want to remove from ts[195] to ts[540]
#add NA's for Peron
xs[17,195:540] <- NA
ys[17,195:540] <- NA

#Moscoso collar fell on 28.01.23 ar 17:20 UTC
#Moscoso released 28.01.23 at 22:30 UTC

#remove from ts[741] to ts[772]
xs[14,741:772] <- NA
ys[14,741:772] <- NA


save(list=c('xs','ys','ts'), file = paste0(outdir,'presedente_xy_10min_level1.RData'))
