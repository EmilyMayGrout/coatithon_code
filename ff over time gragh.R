## this script is for making the fission fusion animation over time

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

#rename coati_id column
colnames(test)[colnames(test) == 'coati_ids$name'] <- 'id'
colnames(test)[colnames(test) == 'ts'] <- 'Time'

library(lubridate)
library(hms)
test$hours <- as_hms(test$Time)


# plot subgroup changes
#this is for the 6th of Jan, need to find a better way of subsetting

ggplot(data = test, aes(x = Time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_x_datetime(limits=c(as.POSIXct("2022-01-06 11:00:00"), as.POSIXct("2022-01-06 23:00:00"))) +
  scale_color_discrete(name="Coati ID", 
                       labels=c("Quasar", "Luna", "Cometa", "Estrella", "Venus", "Lucero", "Gus", "Orbita", "Planeta", "Saturno", "Pluto")) +
  geom_point() +
  geom_line(aes(group = id)) +
  theme_classic()

#---------------------------------------------------------------------
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

library(ggthemes)
#final plot:

#png(height = 800, width = 1600, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m.png'))


ggplot(data = test, aes(x = Time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_color_discrete(name="Coati ID", 
                       labels=c("Quasar", "Luna", "Cometa", "Estrella", "Venus", "Lucero", "Gus", "Orbita", "Planeta", "Saturno", "Pluto")) +
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

