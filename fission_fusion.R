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
gps_file <- "galaxy_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10, 20,30,40, 50, 100)

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

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


#------plot 1: number of sub groups when the radius is changed -----------
png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists.png'))
par(mfrow=c(6,1), mar = c(6,5,1,1))
for (i in 1:length(Rs)){

  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7))

}
dev.off()
 
#------plot 2: number of individuals in each sub group when radius is 50m -----------

png(height = 540, width = 270, units = 'px', filename = paste0(plot_dir,'subgroup_size_hists.png'))
subgroup_data <- get_subgroup_data(xs, ys, R=50)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)


par(mfrow=c(2,1), mar = c(6,5,1,1))
hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))

dev.off()
#making tables to look at numbers of individuals in each split
subgroup_df <- as.data.frame(t(subgroup_counts))
table(paste(subgroup_df[s4,1], subgroup_df[s4,2], subgroup_df[s4,3],subgroup_df[s4,4], sep= '_'))
table(paste(subgroup_df[s3,1], subgroup_df[s3,2], subgroup_df[s3,3], sep= '_'))
table(paste(subgroup_df[s2,1], subgroup_df[s2,2], sep= '_'))


#------------which individuals tend to be in the same subgroup--------------

subgroup_data <- get_subgroup_data(xs, ys, R=50)


ff_net <- matrix(NA, nrow = n_inds, ncol = n_inds)

#going through each dyad and calculating fraction of time they are in the same subgroup (out of all time both are tracked)
for(i in 1:n_inds){
  for(j in 1:n_inds){

    #getting subgroup id for individual i and j
    sub_ids_i <- subgroup_data$ind_subgroup_membership[i,]
    sub_ids_j <- subgroup_data$ind_subgroup_membership[j,]
    
    #computing edge weight (fraction of time in same subgroup)
    ff_net[i,j] <- mean(sub_ids_i == sub_ids_j, na.rm=T)
  }
}

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network.png'))
diag(ff_net) <- NA
image.plot(ff_net, col = viridis(256), zlim=c(0.3,1), xaxt= 'n', yaxt = 'n')

axis(1, at = seq(0,1,length.out= n_inds), labels = coati_ids$name, las = 2)
axis(2, at = seq(0,1,length.out= n_inds), labels = coati_ids$name, las = 2)

points(rep(-.08,n_inds),seq(0,1,length.out=n_inds),col=coati_ids$color, xpd = T, pch = 19)
points(seq(0,1,length.out=n_inds),rep(-.08,n_inds),col=coati_ids$color, xpd = T, pch = 19)
dev.off()












