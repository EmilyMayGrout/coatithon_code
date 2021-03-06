#fission fusion analysis script
#first define groups for each time step (each 10 mins)

#--------priority list----------
#how the radius size influences the distribution of the number of sub groups over time
#who tends to be on their own - plot each individuals subgroup sizes  
#when they split, the distribution of sub group sizes 

#network of within group spatial associations and is that consistent over time
#network of fraction of time each pair is in the same sub group and is it consistent over time
#are these two things related? does the within group network correlate with the fission-fusion subgroup network
#who displaces who - how do the groups split - who moves more

###to do next
#write up lab note book - for the coatithon - with plots and interpretations - go through code again
#make animations of fission fusion with ffmpeg
#once written up, send to Meg and Ben
#start doing some writing to frame the paper
#run the stuff that makes sense on Trago

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

diag(ff_net) <- NA
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]


png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network.png'))

visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
dev.off()


#--------------within full group individual associations----------------
#---------------to compare with the sub group memberships---------------

R = 50
subgroup_data <- get_subgroup_data(xs, ys, R)

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)
#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)
#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs[,full_group_index]
subset_y <- ys[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 10)


new_order <- c(1,11,4,10,2,3,6,7,8,9,5)

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'withingroup_network_withgus.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids[new_order,])
dev.off()


#removing gus from the full group to get proximity data within group



xs_nogus <- xs[-c(5), ]
ys_nogus <- ys[-c(5), ]
coati_ids_nogus <- coati_ids[-c(5),]

R = 50
subgroup_data <- get_subgroup_data(xs_nogus, ys_nogus, R)

n_inds <- nrow(xs_nogus)
n_times <- ncol(xs_nogus)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs_nogus))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)
#find index for when there is only 1 subgroup
s1 <- which(subgroup_data$n_subgroups == 1)

full_group_index <- intersect(all_tracked_idxs, s1)
#subset the x's and y's for the moments the full group has a gps and is together
subset_x <- xs_nogus[,full_group_index]
subset_y <- ys_nogus[, full_group_index]

#get proximity network
within_group_data <- get_proximity_data(subset_x, subset_y, 10)


new_order <- c(1,10,4,9,2,3,5,6,7,8)

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'withingroup_network_withoutgus.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids_nogus[new_order,])
dev.off()

#there's an additional 42 data points but the proximity values don't change much
#but gives more data points by not including him which is good
#other option to get more data is looking at all data (not removing points when not all individuals are tracked)







