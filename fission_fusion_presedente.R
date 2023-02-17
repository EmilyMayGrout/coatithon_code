#fission fusion analysis script with Presedente group
#first define groups for each time step (each 10 mins)

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/'
gps_file <- "presedente_xy_10min_level0.RData"
id_file <- 'coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#removing Wildflower
xs <- xs[c(1:20,22),]
ys <- ys[c(1:20,22),]
coati_ids <- coati_ids[-21,]

#removing males and wildflower
xs <- xs[c(1:4,6,7,11:17,19,20,22),]
ys <- ys[c(1:4,6,7,11:17,19,20,22),]
coati_ids <- coati_ids[-c(5,8,9,10,18,21),]





#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#this is very low, likely because the chance that ALL collars get a GPS fix at the same time is less when there are more collars.
png(height = 400, width = 480, units = 'px', filename = paste0(plot_dir,'hist_n_tracked.png'))
hist(n_tracked, breaks = 15 )
dev.off()






#------plot 1: number of sub groups when the radius is changed (graph put in dropbox results folder) -----------
png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists.png'))
par(mfrow=c(6,1), mar = c(6,5,1,1))

for (i in 1:length(Rs)){
  
  R <- Rs[i]
  subgroup_data <- get_subgroup_data(xs, ys, R)
  xlab <- ''
  if(i == length(Rs)){
    xlab <- 'Number of subgroups'
  }
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(.5,22,1), cex.lab = 2, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7))
  
}
dev.off()



#------plot 2: number of individuals in each sub group when radius is 50m -----------

png(height = 540, width = 270, units = 'px', filename = paste0(plot_dir,'subgroup_size_hists_50m.png'))
subgroup_data <- get_subgroup_data(xs, ys, R=50)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)


par(mfrow=c(2,1), mar = c(6,5,1,1))
hist(subgroup_counts[,s2], breaks=seq(0.5,22,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
hist(subgroup_counts[,s3], breaks=seq(0.5,22,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))

dev.off()



#-------plot 2.2 subgroups at 50m omitting times when there are solitary individuals in their own group

# remove instances when the sub-group is 1 individual as this is not really a sub-group

#first make values equal to 1 and 10 equal to NA for the 2 sub groups and change 9 to NA in 3 subgroups
subgroup_data <- get_subgroup_data(xs, ys, R=50)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)

par(mfrow=c(2,1), mar = c(6,5,2,1))
#for 2 subgroups
subgroup_counts_NA <- subgroup_counts
#subgroup_counts_NA[subgroup_counts_NA == "1"] <- NA
#subgroup_counts_NA[subgroup_counts_NA == "10"] <- NA
hist(subgroup_counts_NA[,s2], breaks=seq(0.5,22,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
#for 3 subgroups
subgroup_counts_NA[subgroup_counts_NA == "9"] <- NA
hist(subgroup_counts_NA[,s3], breaks=seq(0.5,22,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))



#----------plot 3: which individuals tend to be in the same subgroup--------------

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
new_order <- c(1:16)
#order 1 - in age/sex class order: new_order <- c(21,1,13,20,7,14,17,6,10,19,3,4,11,12,16,2,15,22,5,8,9,18)

#order 2 - with subgroups rearranged: new_order <- c(17,1,13,14,3,4,16,2,15,11,6,22,12,20,7,19,9,10,5,8,18,21)

#order 3 - without Wildflower: new_order <- c(17,1,13,14,3,4,16,2,15,11,6,21,12,20,7,19,9,10,5,8,18)

#order 3 - without Wildflower AND males: new_order <- c(13,1,9,10,3,4,12,2,11,7,5,16,8,15,6,14)


par(mgp=c(3, 1, 0), mar=c(11,10,4,2))

ffnet_reorder <- ff_net[new_order, new_order]
png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_network_onlyGroup.png'))

visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])


dev.off()


#--------------plot 4: within full group individual associations----------------
#---------------to compare with the sub group memberships---------------
#this is the probability of individuals being in the same sub-group using absolute dyadic distances 

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


new_order <- c(13,1,9,10,3,4,12,2,11,7,5,16,8,15,6,14)

png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'withingroup_network_onlyGroup.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids[new_order,])
dev.off()



#-------------------------------------------------------------------
#Plotting the mean sub-group size for each hour for all individuals

#get mean group size for each timepoint - so need to get the total number of individuals divided by the number of groups

sub_counts <- as.data.frame(t(subgroup_data$subgroup_counts))
colnames(sub_counts) <- c("sub1", "sub2", "sub3", "sub4", "sub5")
n_groups <- as.data.frame(subgroup_data$n_subgroups)
#combine to make dataframe with the number of individuals in each sub group and the number of subgroups
n_subs <- cbind(sub_counts, n_groups)
colnames(n_subs)[colnames(n_subs) == 'subgroup_data$n_subgroups'] <- 'n_groups'
#adding time to n_subs
n_subs <- cbind(n_subs, ts)

#get the hour
n_subs$hour <- hour(n_subs$ts)

#get total number of individuals
n_subs$n_inds <- n_tracked

#calculate the mean group size
n_subs$mean_group_size <- NA

n_subs$mean_group_size <- (n_subs$n_inds/n_subs$n_groups)

#change hour to Panama time
n_subs$panama_time <- n_subs$hour-5
n_subs$date <- as.Date(n_subs$ts)

#save n_subs as an RData object
#save(n_subs, file = "C:/Users/egrout/Dropbox/stats_Franzi/data/n_subs.Rdata")

library("vioplot")
#now plotting mean group size for each hour of the day
png(height = 500, width = 800, units = 'px', filename = paste0(plot_dir, "mean_group_size_violin.png"))
vioplot(n_subs$mean_group_size ~ n_subs$panama_time,  xlab = "panama time", ylab = "mean subgroup size", col = "cyan3")
dev.off()



#subsetting to smaller date range
n_subs_subset <- n_subs %>%  filter(date > as.Date("2021-12-25") & date < as.Date("2022-01-10"))
#looking at mean group size per day 
g1 <- ggplot(data=n_subs_subset, aes(x=`hour`, y=mean_group_size, group = `hour`))+geom_boxplot()+ theme_classic() + facet_wrap(~date, ncol = 4)
#png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir, "mean_group_size.png"))
g1 + stat_summary(fun = median,
                  geom = "line",
                  aes(group = 1),
                  col = "red")
#dev.off()

#look at the number of groups
g2 <- ggplot(data=n_subs_subset, aes(x=`hour`, y=n_groups, group = `hour`))+geom_boxplot()+ theme_classic() + facet_wrap(~date, ncol = 4)

g2 + stat_summary(fun = median,
                  geom = "line",
                  aes(group = 1),
                  col = "red")


#plot number of subgroups over time

plot(ts, subgroup_data$n_subgroups, type = 'l')






