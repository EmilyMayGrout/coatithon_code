## this script is for making the fission fusion animation over time for galaxy group

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/'
gps_file <- "galaxy_xy_10min_level1.RData"
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

#filter test df to 05.01.22 for Alie to look at
start_time <- as.POSIXct("2022-01-05 11:00:00")
end_time <- as.POSIXct("2022-01-05 23:00:00")

# Subset the dataframe to include only rows within the date range
subset_df <- test[test$ts >= start_time & test$ts <= end_time, ]

write.csv(subset_df, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/050122_gal_subgrouping.csv", row.names = FALSE)


#filter test df to 05.01.22 for Alie to look at
start_time <- as.POSIXct("2022-01-06 11:00:00")
end_time <- as.POSIXct("2022-01-06 23:00:00")

# Subset the dataframe to include only rows within the date range
subset_df <- test[test$ts >= start_time & test$ts <= end_time, ]

write.csv(subset_df, file = "C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/060122_gal_subgrouping.csv", row.names = FALSE)


#read.csv("C:/Users/egrout/Dropbox/coatithon/processed/2022/galaxy/050122_gal_subgrouping.csv")

library(dplyr)
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

#png(height = 400, width = 800, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_Jan6th.png'))
#g2 <- 
ggplot(data = test, aes(x = Time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_x_datetime(limits=c(as.POSIXct("2022-01-06 11:00:00"), as.POSIXct("2022-01-06 23:00:00"))) +
  scale_color_manual(name="Coati ID", 
                       labels=c("Quasar", "Luna", "Cometa", "Estrella", "Venus", "Lucero", "Gus", "Orbita", "Planeta", "Saturno", "Pluto"), values = c("#fcfdbf", "#fecf92", "#fe9f6d", "#f7705c", "#de4968", "#b73779", "#8c2981", "#641a80", "#3b0f70", "#140e36", "#000004")) +
  geom_point() +
  geom_line(aes(group = id)) +
  theme_classic()

#saveRDS(g2, file = paste0(plot_dir,"sub_groupings_over_time_50m_Jan6th.rds"))
#readRDS(file = "C:/Users/egrout/Dropbox/coatithon/results/sub_groupings_over_time_50m_Jan6th.rds")
dev.off()

#-----------------------------------------------------------------------------
#day 29.12.21 to 01.01.22 were when the 2 subgroups were not together

g3 <- ggplot(data = test, aes(x = Time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_x_datetime(limits=c(as.POSIXct("2021-12-29 11:00:00"), as.POSIXct("2022-01-01 23:00:00"))) +
  scale_color_manual(name="Coati ID", 
                     labels=c("Quasar", "Luna", "Cometa", "Estrella", "Venus", "Lucero", "Gus", "Orbita", "Planeta", "Saturno", "Pluto"), values = c("#fcfdbf", "#fecf92", "#fe9f6d", "#f7705c", "#de4968", "#b73779", "#8c2981", "#641a80", "#3b0f70", "#140e36", "#000004")) +
  geom_point() +
  geom_line(aes(group = id)) +
  theme_classic()

saveRDS(g3, file = paste0(plot_dir,"sub_groupings_over_time_50m_Dec29th_Jan1st.rds"))
readRDS(file = "C:/Users/egrout/Dropbox/coatithon/results/galaxy_results/level1/sub_groupings_over_time_50m_Dec29th_Jan1st.rds")

png(height = 500, width = 800, units = 'px', filename = paste0(plot_dir,'subgroup_4days.png'))
g3
dev.off()




#---------------------------------------------------------------------
#now adding grey areas for night time

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

library(ggthemes)
#final plot:


png(height = 800, width = 2000, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_2.png'))


#g1 <- 
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
  scale_y_continuous("Sub-group number", limits = c(0, 6),expand=c(0,0), breaks = 0:5) +
  theme(panel.background = element_rect(fill = 'lightsteelblue3'), #changed colour to snow2 for the recursion markdown
        panel.grid.major = element_line(color = 'lightsteelblue3'),  
        panel.grid.minor = element_line(color = 'lightsteelblue3', size = 2)) +
  scale_x_datetime(limits=c(as.POSIXct("2022-01-05 10:00:00"), as.POSIXct("2022-01-06 00:00:00"), tz = "America/Panama"), position = "top", date_breaks="1 day", expand=c(0,0)) +
  xlab("Panama time") +
  theme(axis.text.x=element_text(size=25),
        axis.text.y=element_text(size=30),
        axis.title=element_text(size=30),
        legend.title = element_text(size=35),
        legend.text = element_text(size=35), legend.key=element_rect(fill="white"))
#saveRDS(g1, file = paste0(plot_dir,"sub_groupings_over_time_50m_2.rds"))

#readRDS(file = "C:/Users/egrout/Dropbox/coatithon/results/sub_groupings_over_time_50m_2.rds")

dev.off()


#hist(test$`sub-group`, col = "azure3")







#---------------------------------------------------------------
#for 1 day with the same these as above graph


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
