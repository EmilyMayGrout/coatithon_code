#this script is getting the speed of the subgroup and the full group
#need to rerun this script with level2 speeds

#---PARAMS----
R <- 50
dt <- 30 #time interval between points (=10 for coati low res). 10 is meters per minute.
min_tracked <- 7 #minimum number of individuals tracked to include in analysis (=7)

#---DIRECTORIES----
#which group - galaxy or presedente
group <- 'presedente'

code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'

if(group == 'galaxy'){
  groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
  plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level2/'

  } else if(group == 'presedente'){
  groupdir <- "C:/Users/egrout/Dropbox/coatithon/processed/split_analysis_processed/"
  plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level2/'
  }


#-----LIBRARIES-----

library(fields)
library(viridis)
library(dplyr)
library(hms) 
library(ggplot2)
library(vioplot)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)
library(doBy)
library(FSA)
library(lubridate)
library(hms)

#read in library of functions
setwd(code_dir)
source('coati_function_library.R')


#----LOAD DATA----
#read in coati ids
setwd(groupdir)

#read in timestamp data
load(file=paste0(group,'_xy_highres_level2.RData'))
load(file=paste0(group,'_coati_ids.RData'))

#-----MAIN------
#finding indexes of the adult males to remove from this analysis - as they're not really group members....
males <- which(coati_ids$age == "Adult" & coati_ids$sex == "Male")
#removing these males from the matrices
if(group == "presedente"){
  #remove Wildflower
  wf <- which(coati_ids$name == "Wildflower")
  md <- which(coati_ids$name == "Mandela")
  xs <- xs[-c(males, wf, md),]
  ys <- ys[-c(males, wf, md),]
  coati_ids <- coati_ids[-c(males, wf, md),]
  
}else if(group == "galaxy"){
  xs <- xs[-c(males),]
  ys <- ys[-c(males),]
  coati_ids <- coati_ids[-c(males),]
  
}

#downsample the xs and ys for one point every 30 seconds as the 1Hz speeds will be much higher than real speeds due to GPS error

xs <- xs[,seq(1, ncol(xs), 30)]
ys <- ys[,seq(1, ncol(ys), 30)]
ts <- ts[seq(1, length(ts), 30)]

n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#get the subgroup data when radius is 50m
subgroup_data <- get_subgroup_data(xs, ys, R)

#------------------------

#make a dataframe with time and individual
speed_df <- data.frame(t = rep(1:(n_times-1), n_inds), 
                       UTC_time = rep(ts[1:(n_times-1)], n_inds) ,
                       ind = rep(1:n_inds, each = (n_times-1)), 
                       n_tracked = rep(n_tracked[1:(n_times-1)], n_inds), 
                       subgroup_size = NA, 
                       split = NA, 
                       speed = NA, 
                       context = NA, 
                       name = rep(coati_ids$name, each = n_times-1))

#loop through to calculate subgroup size, whether there is a split, and ind speed
#i=1
for (i in 1:nrow(speed_df)){
  
  #get the current individual and time for that row
  time <- speed_df$t[i]
  ind <- speed_df$ind[i]
  
  #get current subgroups for all individuals at that time
  sub_data <- subgroup_data$ind_subgroup_membership[,time]
  
  #get current subgroup for the focal individual at that time
  ind_subgroup <- subgroup_data$ind_subgroup_membership[ind,time]
  
  #if that individual is not tracked, leave the NA and go to next row
  if(is.na(ind_subgroup)){
    next
  }
  
  #get subgroup size for the focal individual
  subgroup_size <- sum(sub_data == ind_subgroup, na.rm=T)
  
  #store subgroup size in df
  speed_df$subgroup_size[i] <- subgroup_size
  
  #determine whether group is "split" i.e. has at least 2 subgroups with at least 2 members each
  sub_counts <- subgroup_data$subgroup_counts[,time]
  n_subgroups_with_at_least_2_members <- sum(sub_counts >= 2, na.rm=T)
  split <- n_subgroups_with_at_least_2_members >= 2
  
  #store the split in the dataframe
  speed_df$split[i] <- split
  
  #calculate speed with the xs and ys
  ind_xs_current <- xs[ind, time]
  ind_ys_current <- ys[ind, time]
  
  next_time <- time + 1
  ind_xs_next <- xs[ind, next_time]
  ind_ys_next <- ys[ind, next_time]
  
  dx <- ind_xs_next - ind_xs_current
  dy <- ind_ys_next - ind_ys_current
  
  #getting speed
  speed <- (sqrt((dx)^2 + (dy)^2))/dt
  
  speed_df$speed[i] <- speed
  
  #adding context: alone, split or group together
  if(subgroup_size == 1){
    speed_df$context[i] <-  "alone"
  }
  
  if(split == TRUE & subgroup_size > 1){
    speed_df$context[i] <- "split"
  }
  
  if(split == FALSE & subgroup_size > 1){
    speed_df$context[i] <- "together"
  }
  
  
}

speed_days <- speed_df

#remove NA's (when individual was not tracked) or when too few individuals were tracked
speed_days <- speed_days %>%  filter(n_tracked > min_tracked)
speed_days <- speed_days %>% filter(!is.na(speed))

#remove times when alone
#speed_days <- speed_days[!(speed_days$context == "alone"),]

hist(speed_days$speed, breaks = 1000)

vioplot(speed_days$speed ~ speed_days$context)

# function for number of observations 
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


if(group == "galaxy"){
  colors <- c("darkolivegreen2","olivedrab3", "olivedrab4")
  label <- "Galaxy group"
  
} else if (group == "presedente"){
  colors <- c("coral1","coral3", "coral4")
  label <- "Presidente group"
}



#if dt is 30, the speed is meters/second
#should rerun this with the level2 data to remove errenious GPS points
ggplot(speed_days, aes(x = context, y = speed, fill = context)) + 
  geom_violin() +
  theme_classic()+
  scale_fill_manual(values= colors)+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        legend.title = element_text(size=24),
        legend.text = element_text(size=20), 
        strip.text = element_text(size = 20)) +
  labs(x = label, y = "Speed (m/s)")+
  ylim(0, 2)+
  guides(fill="none")+
  stat_summary(fun.data = give.n, geom = "text", cex = 6, position = position_nudge(x=0.2, y = 1.4))+
  facet_wrap(~name)
  
ggsave(paste0(plot_dir, "speeds_highres_byind.png"), width = 25, height = 15)


mean(speed_days$speed[speed_days$context == "together"])
sd(speed_days$speed[speed_days$context == "together"])

mean(speed_days$speed[speed_days$context == "split"])
sd(speed_days$speed[speed_days$context == "split"])

mean(speed_days$speed[speed_days$context == "alone"])
sd(speed_days$speed[speed_days$context == "alone"])







