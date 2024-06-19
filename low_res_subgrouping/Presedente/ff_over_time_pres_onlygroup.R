## this script is for making the fission fusion animation over time for presedente group

##WITHOUT MALES AND WILDFLOWER

#--------PARAMS-------
data_dir <- "C:/Users/egrout/Dropbox/coatithon/processed/2023/presedente/"
code_dir <- 'C:/Users/egrout/Dropbox/coatithon/coatithon_code/'
plot_dir <- 'C:/Users/egrout/Dropbox/coatithon/results/presedente_results/level1/'
gps_file <- "presedente_xy_10min_level1.RData"
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


#REMOVING MALES AND WILDFLOWER -- ids: 5,8,9,10,18,21

xs <- xs[-c(5,8,9,10,18,21),]
ys <- ys[-c(5,8,9,10,18,21),]

coati_ids <- coati_ids[-c(5,8,9,10,18,21),]

#-----MAIN------
n_inds <- nrow(xs)
n_times <- ncol(xs)

#number of individuals tracked at each time point
n_tracked <- colSums(!is.na(xs))

#indexes to time points where all individuals were tracked
all_tracked_idxs <- which(n_tracked==n_inds)

#making the subgroup data into a dataframe so I can make social networks 
R = 50
subgroup_data <- get_subgroup_data(xs, ys, R) #this is correct

t <- as.data.frame(rbind(subgroup_data$ind_subgroup_membership[,1:1165]))


#put id into dataframe
t <- cbind(coati_ids$name, t)
#changing column names into a list from 1 to 1165 - one day would be good to have as time
column_names <- c(1:1165)
colnames(t)[-1] <- column_names

#change ID's into numbers
t[t == "Ardera"] <- '1'
t[t == "Castro" ] <- '2'
t[t == "Cleopatra" ] <- '3'
t[t ==  "Gandhi" ] <- '4'
t[t == "Ghengis Khan"] <- '5'
t[t == "Gillard"  ] <- '6'
t[t == "May" ] <- '7'
t[t == "Meir" ] <- '8'
t[t == "Merkel" ] <- '9'
t[t == "Moscoso" ] <- '10'
t[t == "Mujica" ] <- '11'
t[t == "Obama" ] <- '12'
t[t == "Peron" ] <- '13'
t[t == "Torrijos" ] <- '14'
t[t == "Truss" ] <- '15'
t[t == "Zelenskyy" ] <- '16'


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
  test$`coati_ids$name` == 1  ~ test$`sub-group` -0.24,
  test$`coati_ids$name` == 2  ~ test$`sub-group` -0.21,
  test$`coati_ids$name` == 3  ~ test$`sub-group` -0.18,
  test$`coati_ids$name` == 4  ~ test$`sub-group` -0.15,
  test$`coati_ids$name` == 5  ~ test$`sub-group` -0.12,
  test$`coati_ids$name` == 6  ~ test$`sub-group` -0.06,
  test$`coati_ids$name` == 7  ~ test$`sub-group` -0.03,
  test$`coati_ids$name` == 8  ~ test$`sub-group` -0.0,
  test$`coati_ids$name` == 9  ~ test$`sub-group` -0.03,
  test$`coati_ids$name` == 10  ~ test$`sub-group`-0.06,
  test$`coati_ids$name` == 11  ~ test$`sub-group`+ 0.09,
  test$`coati_ids$name` == 12  ~ test$`sub-group`+ 0.12,
  test$`coati_ids$name` == 13  ~ test$`sub-group`+ 0.15,
  test$`coati_ids$name` == 14  ~ test$`sub-group`+ 0.18,
  test$`coati_ids$name` == 15  ~ test$`sub-group`+ 0.21,
  test$`coati_ids$name` == 16  ~ test$`sub-group`+ 0.24
)

#remove rows with NA's
test <- test[complete.cases(test), ]

#rename coati_id column
colnames(test)[colnames(test) == 'coati_ids$name'] <- 'id'
colnames(test)[colnames(test) == 'ts'] <- 'Time'


library(lubridate)
library(hms)
library(RColorBrewer)
test$hours <- as_hms(test$Time)


#choosing a colour palette

#display.brewer.all(colorblindFriendly = TRUE)

#getting the colour IDs for ggplot
#brewer.pal(n = 11, name = "PiYG")
#"#8E0152" "#C51B7D" "#DE77AE" "#F1B6DA" "#FDE0EF" "#F7F7F7" "#E6F5D0" "#B8E186" "#7FBC41" "#4D9221" "#276419"
#brewer.pal(n = 11, name = "RdBu")
#"#67001F" "#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#F7F7F7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC" "#053061"
#brewer.pal(n = 12, name = "Paired")
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"



#because ggplot makes id numbers in order as though they are a factor, I need to rearrange
as.factor(test$id)
#Levels: 1 10 11 12 13 14 15 16 2 3 4 5 6 7 8 9 <- this the order the ids need to be in

coati_ids$new_order <- 1:nrow(coati_ids)
#"Ardern", "Moscoso", "Mujica","Obama","Peron","Torrijos","Truss","Zelenskyy","Castro","Cleopatra","Gandhi","Ghengis Khan","Gillard","May", "Meir","Merkel"





#plot sub groups over time

png(height = 1000, width = 1300, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_nomales_29th.png'))

ggplot(data = test, aes(x = Time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_x_datetime(limits=c(as.POSIXct("2023-01-29 17:00:00"), as.POSIXct("2023-01-29 23:30:00"))) +
  scale_color_manual(name="Coati ID", 
                     labels= c("Ardern", "Moscoso", "Mujica","Obama","Peron","Torrijos","Truss","Zelenskyy","Castro","Cleopatra","Gandhi","Ghengis Khan","Gillard","May", "Meir","Merkel"), 
                     values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "grey", "#B15928", "#B2182B", "#D6604D", "#F4A582", "cyan4")) +
  geom_point(size = 4) +
  geom_line(aes(group = id)) +
  labs(y = "sub group ID", x = "UTC time")+
  ylim(0, 3)+ 
  geom_hline(linetype = 'dotted', yintercept = 1.5)+ 
  theme( panel.background =element_blank(),
                     legend.key = element_rect(fill = "white"),
                    text=element_text(size=20), #change font size of all text
                    axis.text=element_text(size=20), #change font size of axis text
                    axis.title=element_text(size=20), #change font size of axis titles
                    plot.title=element_text(size=20), #change font size of plot title
                    legend.text=element_text(size=20), #change font size of legend text
                    legend.title=element_text(size=20))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) #change font size of legend title  

dev.off()

#---------------------------------------------------------------------------------------------------------
#now adding grey areas for night time

#change time to Panama time
test$Panama_time <- with_tz(test$Time, tzone = "America/Panama")

#need to make separate dataframe for day and night times then add it in to geom_rect, calling the different dataframe for each geom
firstday <- as.POSIXct('2023-01-19 06:10', tz = 'America/Panama')
lastday <-  as.POSIXct('2023-02-02 18:10', tz = 'America/Panama')
xmax <- seq.POSIXt(from = firstday, to = lastday,  by = 'day')
firstnight <- as.POSIXct('2023-01-19 18:10', tz = 'America/Panama')
lastnight <-  as.POSIXct('2023-02-02 18:10', tz = 'America/Panama')
xmin <- seq.POSIXt(from = firstnight, to = lastnight,  by = 'day')
ymin = 0
ymax = 5

daynight <- data.frame(1:15,xmax, xmin, ymax, ymin)
colnames(daynight)[colnames(daynight) == 'X1.15'] <- 'rect_id'

library(ggthemes)
#final plot:
png(height = 600, width = 2000, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_nomales.png'))

ggplot(data = test, aes(x = Panama_time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_color_discrete(name="Coati ID", 
                       labels=c("Ardern", "Moscoso", "Mujica","Obama","Peron","Torrijos","Truss","Zelenskyy","Castro","Cleopatra","Gandhi","Ghengis Khan","Gillard","May", "Meir","Merkel")) +
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
  scale_y_continuous("sub group ID", limits = c(0, 5),expand=c(0,0), breaks = 0:8) +
  theme(panel.background = element_rect(fill = 'azure3'), #changed colour to snow2 for the recursion markdown
        panel.grid.major = element_line(color = 'azure3'),  
        panel.grid.minor = element_line(color = 'azure3', size = 2)) +
  scale_x_datetime(limits=c(as.POSIXct("2023-01-19 10:00:00"), as.POSIXct("2023-02-02 02:00:00"), tz = "America/Panama"), position = "top", date_breaks="1 day", expand=c(0,0)) +
  xlab("Panama time") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.title = element_text(size=25),
        legend.text = element_text(size=25), legend.key=element_rect(fill="white"))

dev.off()

#----------------------------------------------------------------------------------------------------------------

#just looking at 1 day

png(height = 600, width = 2000, units = 'px', filename = paste0(plot_dir,'sub_groupings_over_time_50m_nomales_29th.png'))

ggplot(data = test, aes(x = Panama_time, 
                        y = subgroup_mod, 
                        color = id, 
                        group = id)) +
  scale_color_discrete(name="Coati ID", 
                       labels=c("Ardern", "Moscoso", "Mujica","Obama","Peron","Torrijos","Truss","Zelenskyy","Castro","Cleopatra","Gandhi","Ghengis Khan","Gillard","May", "Meir","Merkel")) +
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
  scale_y_continuous("sub group ID", limits = c(0, 5),expand=c(0,0), breaks = 0:5) +
  theme(panel.background = element_rect(fill = 'azure3'), #changed colour to snow2 for the recursion markdown
        panel.grid.major = element_line(color = 'azure3'),  
        panel.grid.minor = element_line(color = 'azure3', size = 2)) +
  scale_x_datetime(limits=c(as.POSIXct("2023-01-29 10:00:00"), as.POSIXct("2023-01-30 02:00:00"), tz = "America/Panama"), position = "top", date_breaks="3 hour", expand=c(0,0)) +
  xlab("Panama time") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.title = element_text(size=25),
        legend.text = element_text(size=25), legend.key=element_rect(fill="white"))

dev.off()





