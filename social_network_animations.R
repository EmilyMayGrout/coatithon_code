#this script is to visualise the social networks 

data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/animations/white_background_nolegend/'
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/animations/white_background_withlegend/'
#plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/animations/black_background_withlegend/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/animations/bigger_labels/'
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

dates <- seq.Date(from = as.Date('2021-12-24'), to = as.Date('2022-01-07'), by = 1)

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

for(d in 1:length(dates)){
  
  date <- as.character(dates[d])
  print(date)
  #narrow down to that date
  t_idxs <- which(as.Date(ts)==date)
  
  #-----MAIN------
  
  xs_day <- xs[,t_idxs]
  ys_day <- ys[,t_idxs]
  ts_day <- ts[t_idxs]
  
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #getting window axis to fit entire data
  max_y <- max(ys_day, na.rm = T)
  min_y <- min(ys_day, na.rm = T)
  max_x <- max(xs_day, na.rm = T)
  min_x <- min(xs_day, na.rm = T)
  
  #making animation
  setwd(plot_dir)
  subgroup_data <- get_subgroup_data(xs=xs_day, ys=ys_day, R=50)
  
  id_colors <- rainbow(11)
  
  dir.create(date)
  for (t in 1:ncol(xs_day)){
    
    png(height = 1080, width = 1080, units = 'px', filename = paste0(plot_dir, date, '/', t, '.png'), bg = 'black')
    plot(NULL, xlim = c(min_x, max_x), ylim = c(min_y, max_y))
    
    legend('bottomleft',legend=coati_ids$name, pch = 21, col = id_colors, pt.bg = coati_ids$color, text.col='white' , cex = 2.5, pt.cex = 3, pt.lwd = 3)
    
    subgroup_data$ind_subgroup_membership[,t]
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
            
            lines(c(xi, xj), c(yi, yj), col = 'white')
          }
        }
  
      }
      range_x <- max_x - min_x
      range_y <- max_y - min_y
      lines(c(min_x + 0.9*range_x, min_x + 0.9*range_x-50), c(min_y + range_y*0.05, min_y + range_y*0.05), col = 'white', lwd = 2)
      text(x = min_x + 0.9*range_x - 25, y = min_y + range_y*0.1, labels = '50 m' ,col = 'white', cex = 3)
      points(xs_day[,t], ys_day[,t], pch = 21, bg = coati_ids$color, cex = 3, col = id_colors, lwd = 3 )
    }
    dev.off()
  
  }
  
}
  
  date <- as.character(dates[d])
  print(date)
  #narrow down to that date
  t_idxs <- which(as.Date(ts)==date)
  
  #-----MAIN------
  
  xs_day <- xs[,t_idxs]
  ys_day <- ys[,t_idxs]
  ts_day <- ts[t_idxs]
  
  n_inds <- nrow(xs)
  n_times <- ncol(xs)
  
  #getting window axis to fit entire data
  max_y <- max(ys_day, na.rm = T)
  min_y <- min(ys_day, na.rm = T)
  max_x <- max(xs_day, na.rm = T)
  min_x <- min(xs_day, na.rm = T)
  
  #making animation
  setwd(plot_dir)
  subgroup_data <- get_subgroup_data(xs=xs_day, ys=ys_day, R=50)
  
  id_colors <- rainbow(11)
  
dir.create(date)
for (t in 1:ncol(xs_day)){
  
  png(height = 450, width = 450, units = 'px', filename = paste0(plot_dir, date, '/', t, '.png'), bg = 'white')
  plot(NULL, xlim = c(min_x, max_x), ylim = c(min_y, max_y), xlab = "", ylab = "")
  
  #legend('bottomleft',legend=coati_ids$name, pch = 21, col = id_colors, pt.bg = coati_ids$color, text.col='black' , cex = 2.5, pt.cex = 3, pt.lwd = 3)
  
  subgroup_data$ind_subgroup_membership[,t]
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
          
          lines(c(xi, xj), c(yi, yj), col = 'black')
        }
      }
    }
    range_x <- max_x - min_x
    range_y <- max_y - min_y
    lines(c(min_x + 0.9*range_x, min_x + 0.9*range_x-50), c(min_y + range_y*0.05, min_y + range_y*0.05), col = 'black', lwd = 2)
    text(x = min_x + 0.9*range_x - 25, y = min_y + range_y*0.1, labels = '50 m' ,col = 'black', cex = 1.5)
    points(xs_day[,t], ys_day[,t], pch = 21, bg = coati_ids$color, cex = 3, col = id_colors, lwd = 3 )
   }
  dev.off()
  
 }
