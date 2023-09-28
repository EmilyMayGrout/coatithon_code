
#fission fusion analysis script for Galaxy group 

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/code_review/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
gps_file <- "galaxy_xy_10min_level1.RData" #level0 is when Venus is not removed
id_file <- 'galaxy_coati_ids.RData'

#list of Rs
Rs <- c(10,20,30,40,50,100)
R <- 50

#-------SETUP-------

library(fields)
library(viridis)
library(tidyverse)
library(lubridate)
library(hms)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vioplot)

#read in library of functions
setwd(code_dir)
source('coati_function_library_V1.R')

#load data
setwd(data_dir)
load(gps_file)
load(id_file)

#-----FUNCTIONS----
mode <- function(x) {
  return(as.numeric(names(which.max(table(x)))))
}

#-----MAIN------

n_inds <- nrow(xs)
n_times <- ncol(xs)


#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#----------------------------------------------------------------------

# Figure 2a,b,c: Characterizing the subgroup patterns when group was split into 2 or 3 subgroups
png(height = 400, width = 1280, units = 'px', filename = paste0(plot_dir,'50m_charecterisations.png'))

par(mfrow=c(1,3), mar = c(6,5,2,1)) #(bottom, left, top, right)
subgroup_data <- get_subgroup_data(xs, ys, R)
hist(subgroup_data$n_subgroups[all_tracked_idxs],main = "50m radii", xlab =  'Number of subgroups', col = "darkolivegreen3", breaks = seq(.5,11,1), cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 6), border = "transparent")

subgroup_counts <- subgroup_data$subgroup_counts[,all_tracked_idxs]
n_subgroups <- subgroup_data$n_subgroups[all_tracked_idxs]

s2 <- which(n_subgroups == 2)
s3 <- which(n_subgroups == 3)

hist(subgroup_counts[,s2], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '2 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), border = "transparent")
hist(subgroup_counts[,s3], breaks=seq(0.5,11,1), xlab = 'Subgroup size', main = '3 subgroups', col = "darkolivegreen4", cex.lab = 2, cex.main = 3, cex.axis=2, freq = FALSE, ylim=c(0,.6), xlim = c(0, 11), border = "transparent")

dev.off()

#-----------------------------------------------------------

#Figure S2a: Number of sub groups when the radius is changed (graph put in dropbox results folder)

png(height = 1080, width = 480, units = 'px', filename = paste0(plot_dir,'n_subgroups_hists_level1.png'))
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

#--------------------------------------------------------------

#Figure 3a: which individuals tend to be in the same subgroup

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
new_order <- c(5,1,11,4,10,2, 3,6,7,8,9)
ffnet_reorder <- ff_net[new_order, new_order]

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_level1.png'))

visualize_network_matrix(ffnet_reorder, coati_ids[new_order,])
dev.off()

#--------------------------------------------------------------------------

#Figure 6a: Probability of individuals being in the same sub-group using absolute dyadic distances 
#removing Gus (adult male) from the full group to get proximity data within group

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

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'withingroup_network_withoutgus_level1.png'))
visualize_network_matrix(within_group_data$proximity_net, coati_ids_nogus[new_order,])
dev.off()

#there's an additional 42 data points but the proximity values don't change much
#this gives more data points by not including Gus


#-------------------------------------------------------------------

#Figure S3a: Plotting the mean sub-group size for each hour for all individuals

#get mean group size for each time point - so need to get the total number of individuals divided by the number of groups

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

#now plotting mean group size for each hour of the day
png(height = 500, width = 700, units = 'px', filename = paste0(plot_dir, "mean_group_size_violin_level1_red.png"))
vioplot(n_subs$mean_group_size ~ n_subs$panama_time,  cex.axis = 1.5, cex.lab = 1, xlab = "panama time", ylab = "mean subgroup size", col = "darkolivegreen3", colMed = "black")
dev.off()


#---------------------------------------------------------------------

#Figure S4ab: make plot for the number of GPS points recorded for each individual

png(height = 600, width = 1400, units = 'px', filename = paste0(plot_dir, "number_tracked.png"))
par(mfrow=c(1,2), mar = c(8,8,2,1)) #c(bottom, left, top, right)
sum_tracked <- colSums(!is.na(xs))
hist(sum_tracked[ !(sum_tracked==0) ], ylim = c(0, 800), main = "", xlab = "Number of individuals tracked", ylab = "Frequency", col = "darkolivegreen4", cex.lab = 2, cex.axis = 2)
each_sum <- data.frame(sum = rowSums(!is.na(xs)))
barplot(each_sum$sum, names.arg = coati_ids$name, las=2, col = "darkolivegreen3", ylab = "Number of GPS points",  cex.lab = 2, cex.axis = 2, cex.names=2, mgp=c(5,1,0))
dev.off()

#---------------------------------------------------------------------

#Figure 4a: histogram for consistency 

#need to run this again as Gus was removed for figure 6a
n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, 50)


#run through time by time and check if each time is a split
splits <- c()
for(t in 1:(n_times-1)){
  
  #get subgroup membership now and later
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #if we have a time of completely nas, pass
  if(sum(!is.na(subgroups_now))==0 | sum(!is.na(subgroups_later))==0){
    next
  }
  
  #transfer NAs from now to later and later to now
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #get number of subgroups now and in next time step (later)
  n_subgroups_now <- length(unique(subgroups_now[!is.na(subgroups_now)]))
  n_subgroups_later <- length(unique(subgroups_later[!is.na(subgroups_later)]))
  
  #get number of singleton groups now and later
  singletons_now <- sum(table(subgroups_now)==1)
  singletons_later <- sum(table(subgroups_later)==1)
  
  #determine if this time step is a split
  #if we have one group that goes to more than one, and there are no singletons subsequently, then it's a split
  if((n_subgroups_now == 1 )
     & n_subgroups_later > n_subgroups_now
     & singletons_later==0
  ){
    splits <- c(splits, t)
  }
  
  #if we have more than one group, but rest are singletons, and number of singletons doesn't change (so we don't have just one loner moving off), then it's a split
  if(n_subgroups_now > 1 
     & ((singletons_now + 1) == n_subgroups_now )
     & n_subgroups_later > n_subgroups_now
     & singletons_now == singletons_later
  ){
    splits <- c(splits, t)
  }
  
  
}

#make a data frame of splits, with original group and subgroups
splits_df <- data.frame(t = splits, orig_group=NA, sub1=NA, sub2=NA, sub3=NA, sub4=NA, sub5=NA)
for(i in 1:nrow(splits_df)){
  t <- splits_df$t[i]
  subgroups_now <- subgroup_data$ind_subgroup_membership[,t]
  subgroups_later <- subgroup_data$ind_subgroup_membership[,t+1]
  
  #transfer NAs from now to later and later to now so errenous splits are not detected
  subgroups_now[which(is.na(subgroups_later))] <- NA
  subgroups_later[which(is.na(subgroups_now))] <- NA
  
  #original group = largest subgroup (no singletons)
  orig_subgroup_number <- mode(subgroups_now)
  orig_subgroup_members <- which(subgroups_now == orig_subgroup_number)
  
  #store original group membership in data frame
  splits_df$orig_group[i] <- list(orig_subgroup_members)
  
  #find the groups where the original members went
  group_ids_later <- unique(subgroups_later[orig_subgroup_members])
  
  for(j in 1:length(group_ids_later)){
    group_id <- group_ids_later[j]
    inds_in_group <- which(subgroups_later==group_id)
    orig_inds_in_group <- intersect(inds_in_group, orig_subgroup_members) #only count the original group members 
    
    #really hacky shit to get R to put lists into a data frame :(
    if(j==1){
      splits_df$sub1[i] <- list(orig_inds_in_group)
    }
    if(j==2){
      splits_df$sub2[i] <- list(orig_inds_in_group)
    }
    if(j==3){
      splits_df$sub3[i] <- list(orig_inds_in_group)
    }
    if(j==4){
      splits_df$sub4[i] <- list(orig_inds_in_group)
    }
    if(j==5){
      splits_df$sub5[i] <- list(orig_inds_in_group)
    }
  }
}

#number in each subgroup
splits_df$n_orig <- sapply(splits_df$orig_group, function(x){return(sum(!is.na(x)))})
splits_df$n_sub1 <- sapply(splits_df$sub1, function(x){return(sum(!is.na(x)))})
splits_df$n_sub2 <- sapply(splits_df$sub2, function(x){return(sum(!is.na(x)))})
splits_df$n_sub3 <- sapply(splits_df$sub3, function(x){return(sum(!is.na(x)))})

#should save this dataframe as its needed for the merge_analysis_galaxy 
save(splits_df, file = "C:/Users/egrout/Dropbox/coatithon_notgithub/Galaxy_fission_fusion/splits_df.Rdata")  

#DONE WITH DATAFRAME

p_dyad_together <- get_p_dyad_together(splits_df_local = splits_df, n_inds_local = n_inds)

#Compute consistency based on consistency metric (see coati_function_library)
consistency_data <- get_consistency(p_dyad_together)

#randomize splits and recompute consistency metric
n_rands <- 1000
rando_consistencies <- rep(NA, n_rands)
for (i in 1:n_rands){
  
  rando <- randomise_splits(splits_df)
  rando_dyad <- get_p_dyad_together(splits_df_local = rando, n_inds_local = n_inds)
  consist <- get_consistency(rando_dyad)
  rando_consistencies[i] <- consist
  
}

png(height = 600, width = 900, units = 'px', filename = paste0(plot_dir,'consistency_splits_hist.png'))
par(mfrow=c(1,1), mar = c(6,6,3,3))#(bottom, left, top, right)
hist(rando_consistencies, breaks=seq(0,0.5,0.005), main = "", xlab = "Consistency of random sub-group allocations", col = "slategray4",  cex.lab = 2.5, cex.axis=2.5)
abline(v=consistency_data, col = 'orange2', lwd=4)
dev.off()


#this graph shows that the mean consistency of individuals splitting with others is more likely than randomly assigning individuals into sub-groups


#----------------------------------------------------------------

#Figure 3b: consistency matrix

#should look into repeatability of binary data - to see if there is a better metric for getting consistency values from binary data

#symmetrize matrix - copy entries from p_dyad_together[i,j] to p_dyad_together[j,i]
for(i in 1:(n_inds-1)){
  for(j in (i+1):n_inds){
    p_dyad_together[j,i] <- p_dyad_together[i,j]
  }
}


diag(p_dyad_together) <- NA
new_order <- c(5,1,11,4,10,2,3,6,7,8,9)
p_dyad_together_reorder <- p_dyad_together[new_order, new_order]

png(height = 600, width = 650, units = 'px', filename = paste0(plot_dir,'subgroup_network_splits.png'))
par(mfrow=c(1,1), mar = c(1,2,1,1))#(bottom, left, top, right)
visualize_network_matrix(p_dyad_together_reorder, coati_ids[new_order,])
dev.off()

#---------------------------------------------------------------------

#Figure 5: Visualisation of the temporal fission-fusion dynamics for one day in Galaxy group. 

#making the subgroup data into a dataframe so I can make social networks 
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


#pivot the data frame into a long format
test <- t %>% pivot_longer(-c(`coati_ids$name`), names_to='time', values_to='sub-group')
test$time <- as.numeric(test$time)

#time df
time_df <- data.frame(ts, 1:length(ts))
names(time_df)[2] <- "time"
test <- left_join(test, time_df)

test$subgroup_mod <- case_when(
  test$`coati_ids$name` == 1  ~ test$`sub-group` -0.32,
  test$`coati_ids$name` == 2  ~ test$`sub-group` -0.24,
  test$`coati_ids$name` == 3  ~ test$`sub-group` -0.18,
  test$`coati_ids$name` == 4  ~ test$`sub-group` -0.12,
  test$`coati_ids$name` == 5  ~ test$`sub-group` -0.06,
  test$`coati_ids$name` == 6  ~ test$`sub-group`+ 0.0,
  test$`coati_ids$name` == 7  ~ test$`sub-group`+ 0.06,
  test$`coati_ids$name` == 8  ~ test$`sub-group`+ 0.12,
  test$`coati_ids$name` == 9  ~ test$`sub-group`+ 0.18,
  test$`coati_ids$name` == 10  ~ test$`sub-group`+ 0.24,
  test$`coati_ids$name` == 11  ~ test$`sub-group`+ 0.33
  
)

#remove rows with NA's
test <- test[complete.cases(test), ]

#rename coati_id column
colnames(test)[colnames(test) == 'coati_ids$name'] <- 'id'
colnames(test)[colnames(test) == 'ts'] <- 'Time'

test$hours <- as_hms(test$Time)

#change time to Panama time
test$Panama_time <- with_tz(test$Time, tzone = "America/Panama")

#need to make separate dataframe for day and night times then add it in to geom_rect, calling the different dataframe for each geom
firstday <- as.POSIXct('2021-12-24 06:00', tz = 'America/Panama')
lastday <-  as.POSIXct('2022-01-13 18:00', tz = 'America/Panama')
xmax <- seq.POSIXt(from = firstday, to = lastday,  by = 'day')
firstnight <- as.POSIXct('2021-12-24 18:00', tz = 'America/Panama')
lastnight <-  as.POSIXct('2022-01-13 18:00', tz = 'America/Panama')
xmin <- seq.POSIXt(from = firstnight, to = lastnight,  by = 'day')
ymin = 0
ymax = 6
daynight <- data.frame(1:21,xmax, xmin, ymax, ymin)
colnames(daynight)[colnames(daynight) == 'X1.17'] <- 'rect_id'

test$hour <- as_hms(test$Panama_time)

png(height = 600, width = 1200, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_2_1day.png'))

ggplot(data = test, aes(x = Panama_time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_color_discrete(name="Coati ID", 
                       labels=c("Quasar", "Luna", "Cometa", "Estrella", "Venus", "Lucero", "Gus", "Orbita", "Planeta", "Saturno", "Pluto")) +
  geom_rect(data = daynight, 
            aes(xmin = xmax, xmax = xmin, ymin = ymin, ymax = ymax), 
            inherit.aes = FALSE, fill = "white") +
  
  geom_point(data = test, aes(x = Panama_time, 
                              y = subgroup_mod, 
                              color = id, 
                              group = id), size = 1.9) +
  geom_line(data = test, aes(x = Panama_time, 
                             y = subgroup_mod, 
                             color = id, 
                             group = id)) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous("Sub-group number", limits = c(-0.1, 6),expand=c(0,-0.2), breaks = 0:5) +
  theme(panel.background = element_rect(fill = 'lightsteelblue3'), #changed colour to snow2 for the recursion markdown
        panel.grid.major = element_line(color = 'lightsteelblue3'),  
        panel.grid.minor = element_line(color = 'lightsteelblue3', size = 2))+
  scale_x_datetime(limits=c(as.POSIXct("2022-01-05 11:00:00"), as.POSIXct("2022-01-06 01:00:00"), tz = "America/Panama"), position = "top", date_breaks="6 hour", expand=c(0,0)) +
  xlab("Panama time") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.title = element_text(size=25),
        legend.text = element_text(size=25), legend.key=element_rect(fill="white"))

dev.off()





