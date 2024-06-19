# this code is to use the output matrix from the preprocess_gps_lowres and remove the data from times when either:
#1. the ind was not in the group because of capture/handling effect (e.g. Venus)
#2. the ind's collar came off
#3. the gps points were goofy (e.g. Gillard on 29.01.23 from Presedente going in a straight line)

#output of this code will make a new gps_file called "galaxy_xy_10min_level1.RData"

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/'
gps_file <- "galaxy_xy_10min_level0.RData" 
id_file <- 'coati_ids.RData'

outdir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"


#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#want to add NA's for the data we want to exclude from the analysis
#in Galaxy group this was when Venus responded badly to the collaring and stayed in a tree for 2 days 

#finding time when Venus left the tree

#Venus is xs[3,]

#checked by plotting everyone for the first day ish to see where Venus (black) was
plot(xs[1,100:200], ys[1,100:200], type = "l", col = "blue")
lines(xs[2,100:200], ys[2,100:200], type = "l", col = "red")
lines(xs[3,100:200], ys[3,100:200], type = "l", col = "black")
lines(xs[4,100:200], ys[4,100:200], type = "l", col = "orange")
lines(xs[5,100:200], ys[5,100:200], type = "l", col = "green")
lines(xs[6,100:200], ys[6,100:200], type = "l", col = "purple")
lines(xs[7,100:200], ys[7,100:200], type = "l", col = "pink")
lines(xs[8,100:200], ys[8,100:200], type = "l", col = "grey")
lines(xs[9,100:200], ys[9,100:200], type = "l", col = "cyan")
lines(xs[10,100:200], ys[10,100:200], type = "l", col = "pink4")
lines(xs[11,100:200], ys[11,100:200], type = "l", col = "yellow4")

#when did she start moving with the group?
#from looking at the movevis videos, I know she wasn't moving on the 27.12.21 which is the point at ts[250], started moving on the 28.12.21 at 12:35 ts[322], and was with a subgroup with Orbita [6] on the 29.12.21 ts[391]. She likely found the group on the 28.12.21 

lines(xs[6,356:360], ys[6,356:360], type = "l", col = "red")
plot(xs[3,322:360], ys[3,322:360], type = "l", col = "green")

#[355:360] Venus was close to Orbita by 20m
plot(xs[6,200:350], ys[6,200:350], type = "l", col = "yellow3")
lines(xs[6,335:360], ys[6,335:360], type = "l", col = "red")
lines(xs[3,300:360], ys[3,300:360], col = "green")

#met at ts[342] ish which is 28.12.21 15:50 UTC, looks like the group found Venus and she went with them after

#so going to add NA's for Venus from ts[1] to ts[341]
xs[3,1:341] <- NA
ys[3,1:341] <- NA

save(list=c('xs','ys','ts'), file = paste0(outdir,'galaxy_xy_10min_level1.RData'))









