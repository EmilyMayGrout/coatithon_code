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
#make animations of fission fusion with ffmpeg
#start doing some writing to frame the paper
#run the stuff that makes sense on Trago

#further ideas:
#is one of the sub-groups displaced or do they naturally drift apart from loss of cohesion?
#in the grouping matrix, it could hide individuals who are closely associated because of their association with a different individual,
#how can we see whether this is occuring?
#is there multi-level sub grouping patterns? - groups in groups?
#can we predict a fission by looking at the stability of the group?


#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/'
gps_file <- "galaxy_xy_10min_level0.RData"
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

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)


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
  hist(subgroup_data$n_subgroups[all_tracked_idxs],main = paste(R, "m"), xlab = xlab, col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2, cex.main = 2, cex.axis=2, freq = FALSE, ylim=c(0,.7))

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
hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))

dev.off()


#------plot 2.1: subgroup size at 100m----------------------

png(height = 540, width = 270, units = 'px', filename = paste0(plot_dir,'subgroup_size_hists_100m.png'))
subgroup_data <- get_subgroup_data(xs, ys, R=100)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)


par(mfrow=c(2,1), mar = c(6,5,1,1))
hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))

dev.off()

#not much difference between the 50m and 100m radius cut of for sub group structure 

#making tables to look at numbers of individuals in each split
subgroup_df <- as.data.frame(t(subgroup_counts))
table(paste(subgroup_df[s4,1], subgroup_df[s4,2], subgroup_df[s4,3],subgroup_df[s4,4], sep= '_'))
table(paste(subgroup_df[s3,1], subgroup_df[s3,2], subgroup_df[s3,3], sep= '_'))
table(paste(subgroup_df[s2,1], subgroup_df[s2,2], sep= '_'))


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
subgroup_counts_NA[subgroup_counts_NA == "1"] <- NA
subgroup_counts_NA[subgroup_counts_NA == "10"] <- NA
hist(subgroup_counts_NA[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
#for 3 subgroups
subgroup_counts_NA[subgroup_counts_NA == "9"] <- NA
hist(subgroup_counts_NA[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))


#--------plot 2.3 subgroups at 100m omitting times when there are solitary individuals in their own group

subgroup_data <- get_subgroup_data(xs, ys, R=100)
subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)
s4 <- which(n_subgroups == 4)

par(mfrow=c(2,1), mar = c(6,5,2,1))
#for 2 subgroups
subgroup_counts_NA <- subgroup_counts
subgroup_counts_NA[subgroup_counts_NA == "1"] <- NA
subgroup_counts_NA[subgroup_counts_NA == "10"] <- NA
hist(subgroup_counts_NA[,s2], breaks=seq(0.5,11,1), xlab = '', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))
#for 3 subgroups
subgroup_counts_NA[subgroup_counts_NA == "9"] <- NA
hist(subgroup_counts_NA[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 1.5, cex.main = 1.5, cex.axis=1.5, freq = FALSE, ylim=c(0,.6))




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
new_order <- c(1,11,4,10,2, 3,6,7,8,9,5)
ffnet_reorder <- ff_net[new_order, new_order]


png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'subgroup_network.png'))

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


new_order <- c(1,11,4,10,2,3,6,7,8,9,5)

png(height = 400, width = 400, units = 'px', filename = paste0(plot_dir,'withingroup_network_withgus.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids[new_order,])
dev.off()


#to visualise the absolute dyadic distances for any time 
#png(height = 800, width = 800, units = 'px', filename = paste0(plot_dir,'withingroup_network_withgus_distovertime_2.png'))
#par(mfrow=c(10,11), mar = c(1,1,1,1))
#for (i in 145:253){
#  visualize_network_matrix(within_group_data$dist_over_time[,,i], coati_ids[new_order,])
#the darker the colour, the closer the distance (opposite)
#  }
#dev.off()



#---------------plot 5: full group associations without Gus 
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





#plotting number of subgroups over time when the radius is 50m
plot(subgroup_data$n_subgroups, xlim = c(0, 1633))



#making the subgroup data into a dataframe so I can make social networks 
R = 50
subgroup_data <- get_subgroup_data(xs, ys, R)
t <- as.data.frame(rbind(subgroup_data$ind_subgroup_membership[,1:1633]))
#put id into dataframe
t <- cbind(coati_ids$name, t)
#changing column names into a list from 1 to 1633 - one day would be good to have as time
column_names <- c(1:1633)
colnames(t)[-1] <- column_names
#change ID's into numbers
t[t == 'Quasar'] <- '1'
t[t == 'Estrella'] <- '2'
t[t == 'Venus'] <- '3'
t[t == 'Lucero'] <- '4'
t[t == 'Gus'] <- '5'
t[t == 'Orbita'] <- '6'
t[t == 'Planeta'] <- '7'
t[t == 'Saturno'] <- '8'
t[t == 'Pluto'] <- '9'
t[t == 'Luna'] <- '10'
t[t == 'Cometa'] <- '11'

library(tidyr)
#pivot the data frame into a long format
test <- t %>% pivot_longer(-c(`coati_ids$name`), names_to='time', values_to='sub-group')
test$time <- as.numeric(test$time)

#time df
time_df <- data.frame(ts, 1:length(ts))
names(time_df)[2] <- "time"
test <- left_join(test, time_df)


library(dplyr)
test$subgroup_mod <- case_when(
        test$`coati_ids$name` == 1  ~ test$`sub-group` -0.25,
        test$`coati_ids$name` == 2  ~ test$`sub-group` -0.20,
        test$`coati_ids$name` == 3  ~ test$`sub-group` -0.15,
        test$`coati_ids$name` == 4  ~ test$`sub-group` -0.10,
        test$`coati_ids$name` == 5  ~ test$`sub-group` -0.05,
        test$`coati_ids$name` == 6  ~ test$`sub-group`+ 0.0,
        test$`coati_ids$name` == 7  ~ test$`sub-group`+ 0.05,
        test$`coati_ids$name` == 8  ~ test$`sub-group`+ 0.10,
        test$`coati_ids$name` == 9  ~ test$`sub-group`+ 0.15,
        test$`coati_ids$name` == 10  ~ test$`sub-group`+ 0.20,
        test$`coati_ids$name` == 11  ~ test$`sub-group`+ 0.25

)


plot(test$time, test$subgroup_mod)


#remove rows with NA's
test <- test[complete.cases(test), ]

library(lubridate)
library(hms)
test$hours <- as_hms(test$ts)




# plot subgroup changes
#test %>% 
  #filter(time < 100) %>% 
  ggplot(data = test, aes(x = ts, 
             y = subgroup_mod, 
             color = `coati_ids$name`, 
             group = `coati_ids$name`)) +
  geom_point() +
  geom_line(aes(group = `coati_ids$name`)) +
  theme_classic()

  
  
  
#now adding grey areas for night time

#need to make separate dataframe for day and night times then add it in to geom_rect, calling the different dataframe for each geom
firstday <- as.POSIXct('2021-12-24 11:00', tz = 'UTC')
lastday <-  as.POSIXct('2022-01-13 23:00', tz = 'UTC')
xmax <- seq.POSIXt(from = firstday, to = lastday,  by = 'day')
firstnight <- as.POSIXct('2021-12-24 23:00', tz = 'UTC')
lastnight <-  as.POSIXct('2022-01-13 23:00', tz = 'UTC')
xmin <- seq.POSIXt(from = firstnight, to = lastnight,  by = 'day')
ymin = 0
ymax = 6
daynight <- data.frame(1:21,xmax, xmin, ymax, ymin)
colnames(daynight)[colnames(daynight) == 'X1.17'] <- 'rect_id'

#rename coati_id column
colnames(test)[colnames(test) == 'coati_ids$name'] <- 'id'
colnames(test)[colnames(test) == 'ts'] <- 'Time'


 #aes(xmin = test$ts[which(test$hours == as_hms("23:00:00"))], 
  #              xmax = test$ts[which(test$hours == as_hms("11:00:00"))],
   #             ymin = 0, 
    #            ymax = 5)


t[t == 'Quasar'] <- '1'
t[t == 'Estrella'] <- '2'
t[t == 'Venus'] <- '3'
t[t == 'Lucero'] <- '4'
t[t == 'Gus'] <- '5'
t[t == 'Orbita'] <- '6'
t[t == 'Planeta'] <- '7'
t[t == 'Saturno'] <- '8'
t[t == 'Pluto'] <- '9'
t[t == 'Luna'] <- '10'
t[t == 'Cometa'] <- '11'


library(ggthemes)
#final plot:

#png(height = 800, width = 1600, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m.png'))


ggplot(data = test, aes(x = Time, 
             y = subgroup_mod, 
             color = id, 
             group = id)) +
  scale_color_discrete(name="Coati ID", labels=c("Quasar", "Luna", "Cometa", 
                                             "Estrella", "Venus", "Lucero", 
                                             "Gus", "Orbita", "Planeta",
                                             "Saturno", "Pluto")) +
  geom_rect(data = daynight, 
            aes(xmin = xmax, xmax = xmin, ymin = ymin, ymax = ymax), 
            inherit.aes = FALSE, fill = "white") +
  
  geom_point(data = test, aes(x = Time, 
                              y = subgroup_mod, 
                              color = id, 
                              group = id), size = 1.8) +
  geom_line(data = test, aes(x = Time, 
                             y = subgroup_mod, 
                             color = id, 
                             group = id)) +
  scale_y_continuous("Sub-group number", limits = c( 0, 6 ), breaks = 0:5) +
  theme(panel.background = element_rect(fill = 'snow2'), 
        panel.grid.major = element_line(color = 'snow2'),  
        panel.grid.minor = element_line(color = 'snow2', size = 2)) 
  

#dev.off()


#next make plot for specific days - zoomed in

#first plot for 6th Jan
