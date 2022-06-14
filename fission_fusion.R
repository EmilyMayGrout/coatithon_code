#fission fusion analysis script
#first define groups for each time step (each 10 mins)

#--------priority list----------
#how the radius size influences the distribution of the number of sub groups over time
#who tends to be on their own - plot each individuals subgroup sizes  
#when they split, the distribution of sub group sizes 

#network of within group spatial associations and is that consistent over time
#network of fraction of time each pair is in the same sub group and is it consistent over time
#are these two things related? does the within group network correlate with the fission-fusion subgroup network


#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
in_file <- "galaxy_xy_10min_level0.RData"

#list of Rs
Rs <- c(10, 20,30,40, 50, 100)

#-------SETUP-------

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(in_file)

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists.png'))
par(mfrow=c(6,1), mar = c(4,5,1,1))
for (i in 1:length(Rs)){

  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(1,11), cex.lab = 2, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7))

}
dev.off()
 