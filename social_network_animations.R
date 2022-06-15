#this script is to visualise the social networks 

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/animations/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

date <- '2021-12-24'

#-------SETUP-------

library(fields)
library(viridis)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#narrow down to that date
t_idxs <- which(as.Date(ts)==date)

#-----MAIN------

xs_day <- xs[,t_idxs]
ys_day <- ys[,t_idxs]
ts_day <- ts[,t_idxs]

n_inds <- nrow(xs)
n_times <- ncol(xs)

#getting window axis to fit entire data
max_y <- max(ys_day, na.rm = T)
min_y <- min(ys_day, na.rm = T)
max_x <- max(xs_day, na.rm = T)
min_x <- min(xs_day, na.rm = T)

#making animation
setwd(plot_dir)
subgroup_data <- get_subgroup_data(xs_day, ys_day, R=50)
i=100

for (t in 1:ncol(xs_day)){
  
  png(height = 1080, width = 1080, units = 'px', filename = paste0(plot_dir,t, '.png'))
  plot(xs_day[,t], ys_day[,t], xlim = c(min_x, max_x), ylim = c(min_y, max_y))
  
  subgroup_data$ind_subgroup_membership[,i]
  for(i in 1:n_inds){
    for(j in 1:n_inds){
      
      #getting subgroup id for individual i and j
      sub_id_i <- subgroup_data$ind_subgroup_membership[i,t]
      sub_id_j <- subgroup_data$ind_subgroup_membership[j,t]
      
      #make a line between i and j
      if(!is.na(sub_id_i) & !is.na(sub_id_j)){
        if(sub_id_i == sub_id_j){
          xi <- xs_day[i,t]
          yi <- ys_day[i,t]
          xj <- xs_day[j,t]
          yj <- ys_day[j,t]
          
          lines(c(xi, xj), c(yi, yj))
        }
      }

    }
  }
  dev.off()

}



