#this script is getting behavioural classes from the tracks using tau and eta as well as speed and turning angle and the FPT


library(move)
library(lubridate)
library(MASS)
library(smoove)
library(plyr)
library(doParallel)

options(digits=15)


plot_dir <- "C:/Users/egrout/Dropbox/stats_Franzi/plots/"
data_dir <- "C:/Users/egrout/Dropbox/stats_Franzi/data/"

#this data is from lowres embc script in the animove folder in dropbox
load(file = paste0(data_dir,"resample_Galaxy.Rdata"))

coati_id=read.csv(paste0(data_dir, "coati_id.csv"), header = F)


#convert resample_Galaxy to a movestack
#load(file = paste0(data_dir,"resample_Galaxy.Rdata"))
resamp_gal <- moveStack(resample_Galaxy)
resamp_gal <- spTransform(resamp_gal, CRS("+init=epsg:32617"))
#--------------------------------------------------------
# make a dataframe with ID, time, speed, age and sex
#first collapse the list
Galaxy_df <- as.data.frame(resamp_gal)
names(Galaxy_df)
Galaxy_df <- Galaxy_df[,c(35, 36, 27, 28, 29, 30, 31, 32)]
#read in the coati ids with the age and sex class
#coati_id <- read.csv(paste0(data_dir, "coati_id.csv"), header = F)
#add the age and sex to each individual using the match function
Galaxy_df$age <- coati_id$V3[match(Galaxy_df$trackId, coati_id$V1)]
Galaxy_df$sex <- coati_id$V4[match(Galaxy_df$trackId, coati_id$V1)]
#define the age and sex classes as factors
Galaxy_df$age_logical <- factor(Galaxy_df$age, levels = c("Juvenile", "Sub-adult", "Adult"))
Galaxy_df$age <- factor(Galaxy_df$age)
Galaxy_df$sex <- as.factor(Galaxy_df$sex)
#add column for hour
Galaxy_df$hours <- factor(hour(Galaxy_df$timestamps))

##Making the coordinates a complex number for convinenince purposes
Galaxy_df$Z=Galaxy_df$coords.x1.1+1i*Galaxy_df$coords.x2.1

colnames(Galaxy_df)[7]="utm.easting"
colnames(Galaxy_df)[8]="utm.northing"

Galaxy_df$local_time=with_tz(Galaxy_df$timestamps, tzone = "America/Panama")
Galaxy_df$hours_local=hour(Galaxy_df$local_time)
Galaxy_df2=Galaxy_df[which(Galaxy_df$hours_local>5 & Galaxy_df$hours_local<19),]
Galaxy_df2=Galaxy_df2[-which(date(Galaxy_df2$local_time)==as_date("2021-12-25")),]
Galaxy_df2=Galaxy_df2[-which(date(Galaxy_df2$local_time)==as_date("2021-12-26")),]
Galaxy_df2$Date=date(Galaxy_df2$local_time)
Galaxy_df2$interpolated="No"
unique_dates=unique(Galaxy_df2$Date)

Galaxy_df2_split=split(Galaxy_df2,Galaxy_df2$trackId)
interpolated_data=c()
count=1
for(i in 1:length(Galaxy_df2_split)){
  for(j in 2:length(unique_dates)){
      
    temp=Galaxy_df2_split[[i]]
    
    sleep_loc=temp$Z[which(temp$Date==unique_dates[j-1])][length(which(temp$Date==unique_dates[j-1]))]
    wake_loc=temp$Z[which(temp$Date==unique_dates[j])][1]
    
    new_loc=mean(c(sleep_loc,wake_loc))
    
    sleep_time=max(temp$local_time[which(temp$Date==unique_dates[j-1])], na.rm = TRUE)
    wake_time=min(temp$local_time[which(temp$Date==unique_dates[j])], na.rm = TRUE)
    
    newtimes= seq( from = sleep_time, to = wake_time, by = "10 min")
    newtimes=newtimes[-c(1,length(newtimes))]
    drift = rnegbin(n=length(newtimes), mu = 5, theta = 2)
    directions= runif(n=length(newtimes), min = -pi, max = pi)
    offsets= complex(modulus  = drift, argument = directions)
    tempDF= data.frame(matrix(ncol = ncol(Galaxy_df2), nrow = length(newtimes)))
    colnames(tempDF)=colnames(Galaxy_df2)
    
    New_points=new_loc+offsets
    tempDF$local_time=newtimes
    tempDF$Z=New_points
    tempDF$interpolated = "Yes"
    tempDF$trackId=names(Galaxy_df2_split[i])
    tempDF$utm.easting=Re(New_points)
    tempDF$utm.northing=Im(New_points)
    interpolated_data[count]=list(tempDF)
    count=count+1
  }
}

interpolated_data=do.call(rbind,interpolated_data)
Galaxy_df2$trackId=as.character(Galaxy_df2$trackId)
interpolated_data$trackId=as.character(interpolated_data$trackId)

Galaxy_df3=rbind(Galaxy_df2,interpolated_data)
Galaxy_df3$trackId=as.factor(Galaxy_df3$trackId)

Galaxy_df3=Galaxy_df3[
  with(Galaxy_df3, order(as.character(Galaxy_df3$trackId), Galaxy_df3$local_time)),
]

splitdata=split(Galaxy_df3,Galaxy_df3$trackId)

cl <- makeCluster(detectCores())
registerDoParallel(cl)



# #for loop to get tau and eta for each track - this takes a long time to run so have saved it and reloaded it after forloop
# tracks=c()
# try(for( i in 1:length(splitdata)){
# 
#   track=as.data.frame(splitdata[[i]])
#   track <- track[order(track$local_time),] 
#   track=track[complete.cases(track$local_time),]
#   
#   t=track$local_time
# 
#   simSweep1 <- try(sweepRACVM(Z=track$Z, T=track$local_time, windowsize = 300, windowstep = 10, model = "UCVM", progress=FALSE, time.unit = "mins",.parallel = TRUE))
# 
#   
#   CP1.all <- try(findCandidateChangePoints(windowsweep = simSweep1, clusterwidth = 70))
#   
# 
#   simCP1.table=try(getCPtable(CPs = CP1.all, modelset = "UCVM", tidy = FALSE, iterate = FALSE))
#  
#   simPhases1 <- estimatePhases(simCP1.table)
#   phaseTable1 <- summarizePhases(simPhases1)
#   
#   
# ##Loop to append the eta and tau estimates back to original track data
#   for ( j in 1:nrow(phaseTable1)){
#     for ( k in 1:nrow(track)){
#       
#       if((track$local_time[k]>=phaseTable1$start[j]) & (track$local_time[k]<phaseTable1$end[j])) {track$phase[k]= phaseTable1$phase[j]}
#       
#       if((track$local_time[k]>=phaseTable1$start[j]) & (track$local_time[k]<phaseTable1$end[j])) {track$tau[k]= phaseTable1$tau[j]}
#       
#       if((track$local_time[k]>=phaseTable1$start[j]) & (track$local_time[k]<phaseTable1$end[j])) {track$eta[k]= phaseTable1$eta[j]}
#       
#     }
#     
#   }
#   tracks[i]=list(track)
# })
# tracks2=plyr::ldply(tracks, data.frame)
# tracks2 = tracks2[-which(tracks2$interpolated=="Yes"),]
# 
# splitftp=split(tracks2,tracks2$trackId)
# 
# 
# tracks2$Daylabel=as.numeric(as.factor(as.character(tracks2$Date)))
# tracks2$loopID=paste(tracks2$trackId, tracks2$Daylabel, sep = "_")
# 
# #need to ask Shauhin if I can do this or not....(as the above code takes a long time to run)
# 
# write.csv(tracks2, "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/behav_seg/eta_tau_df.csv")
# 




tracks2 <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/behav_seg/eta_tau_df.csv")

#when loaded need to remove first row and add 
tracks2 <- tracks2[, -1]
tracks2$local_time <- with_tz(tracks2$local_time, tzone = "America/Panama")
tracks2$timestamps <- with_tz(tracks2$timestamps)
tracks2$trackId <- as.factor(tracks2$trackId)
tracks2$age <- as.factor(tracks2$age)
tracks2$sex <- as.factor(tracks2$sex)
tracks2$age_logical <- as.factor(tracks2$age_logical)
tracks2$Date <- date(tracks2$local_time)
tracks2$Daylabel <- as.numeric(tracks2$Daylabel)



library(adehabitatLT)

splitftp=split(tracks2,as.factor(tracks2$loopID))

ftpdata=c()
for ( i in 1:length(splitftp)) {
  ftpday=as.data.frame(splitftp[[i]])
  #colnames(ftpday)=names(data)
  #if(nrow(ftpday)==0){next}
  
  L1.traj <- as.ltraj(xy = ftpday[,c("utm.easting", "utm.northing")], date = ftpday$local_time, id=as.character(ftpday$trackId), typeII=TRUE)
  radii <- (20) #radius to use for fpt
  radii2 <- (50) #radius to use for fpt
  
  fewdaysfpt <- fpt(L1.traj, radii)
  fewdaysfpt2 <- fpt(L1.traj, radii2)
  
  fewdaysfpt1=as.data.frame(fewdaysfpt[])
  fewdaysfpt1b=as.data.frame(fewdaysfpt2[])
  
  ftpday$fpt20=fewdaysfpt1$r1
  ftpday$fpt50=fewdaysfpt1b$r1
  
  ftpdata[i]=list(ftpday)
}
ftpdata=plyr::ldply(ftpdata, data.frame)

ftpdata$logETA=log(ftpdata$eta)
ftpdata$logTAU=log(ftpdata$tau)
ftpdata$log1pTau=log(1+ftpdata$tau)

ftpdata$FPT_above_20=NA	
ftpdata$FPT_above_50=NA

##above 1800 means spending more time in the same place
ftpdata$FPT_above_20[which(ftpdata$fpt20>1800)]=1
ftpdata$FPT_above_20[which(ftpdata$fpt20<=1800)]=0
ftpdata$FPT_above_50[which(ftpdata$fpt50>1800)]=1
ftpdata$FPT_above_50[which(ftpdata$fpt50<=1800)]=0

hist(ftpdata$eta, breaks = 100)
hist(ftpdata$tau, breaks = 100)
hist(ftpdata$fpt20, breaks = 100)
hist(ftpdata$fpt50, breaks = 100)

##ztransform data before clustering
ftpdata$eta_scaled=(ftpdata$eta-mean(ftpdata$eta, na.rm = TRUE))/sd(ftpdata$eta, na.rm = TRUE)
ftpdata$logETA_scaled=(ftpdata$logETA-mean(ftpdata$logETA, na.rm = TRUE))/sd(ftpdata$logETA, na.rm = TRUE)

ftpdata$log1pTau_scaled=(ftpdata$log1pTau-mean(ftpdata$log1pTau, na.rm = TRUE))/sd(ftpdata$log1pTau, na.rm = TRUE)
ftpdata$fpt20_scaled=(ftpdata$fpt20-mean(ftpdata$fpt20, na.rm = TRUE))/sd(ftpdata$fpt20, na.rm = TRUE)
ftpdata$fpt50_scaled=(ftpdata$fpt50-mean(ftpdata$fpt50, na.rm = TRUE))/sd(ftpdata$fpt50, na.rm = TRUE)
ftpdata$turnAngle_scaled=(ftpdata$turnAngle-mean(ftpdata$turnAngle, na.rm = TRUE))/sd(ftpdata$turnAngle, na.rm = TRUE)
ftpdata$speed_ms_scaled=(ftpdata$speed_ms-mean(ftpdata$speed_ms, na.rm = TRUE))/sd(ftpdata$speed_ms, na.rm = TRUE)

library(dbscan)
ftpdata2=ftpdata[complete.cases(ftpdata$fpt20),]
ftpdata2=ftpdata2[complete.cases(ftpdata2$fpt50),]

##Will need to experiment wiht eps and minPts values 
##minPts is the minimum number of points allowed to be in each cluster group

###scaled versions
cluster1 = dbscan(ftpdata2[,c(28,29,31,32)], eps= .9, minPts = 10) ### Clustering with eta, tau, and both FTPs 
cluster2 = dbscan(ftpdata2[,c(28,31,32)], eps= .9, minPts = 10) ## Clustering with eta, tau, and FPT20
cluster3 = dbscan(ftpdata2[,c(29,31,32)], eps= 9, minPts = 10) #Clustering with eta, tau, and FPT50
cluster4 = dbscan(ftpdata2[,c(35,36,28,29)], eps= .6, minPts = 10) ### Clustering with speed, turnangle, and both FTPs 
cluster5 = dbscan(ftpdata2[,c(35,36,28)], eps= .6, minPts = 10) ### Clustering with speed, turnangle, and  FTP20 
cluster6 = dbscan(ftpdata2[,c(35,36,29)], eps= .6, minPts = 10) ### Clustering with speed, turnangle, and FTP50 
cluster7 = dbscan(ftpdata2[,c(35,36)], eps= .6, minPts = 10) ### Clustering with speed, turnangle 
cluster8 = dbscan(ftpdata2[,c(31,32)], eps= .9, minPts = 10) ### Clustering with eta, tau 
cluster9 = dbscan(ftpdata2[,c(28,29)], eps= .9, minPts = 10) ### Clustering with only fpt
cluster10 = dbscan(ftpdata2[,c(28,29,31,32,35,36)], eps= 1, minPts = 10) ### Clustering with only fpt
cluster11 = dbscan(ftpdata2[,c(33,34,35,36,28,29)], eps= 1, minPts = 10) ### Clustering with speed, turnangle, and both FTPs 
cluster12 = dbscan(ftpdata2[,c(36,28,29)], eps= .6, minPts = 10) ### Clustering with speed, turnangle, and both FTPs 




#### Attach the predicted cluster assignements to the original data

#scaled
prediction1 = predict(cluster1, data = ftpdata2[,c(28,29,31,32)], newdata = ftpdata2[,c(28,29,31,32)])
prediction2 = predict(cluster2, data = ftpdata2[,c(28,31,32)], newdata = ftpdata2[,c(28,31,32)])
prediction3 = predict(cluster3, data = ftpdata2[,c(29,31,32)], newdata = ftpdata2[,c(29,31,32)])
prediction4 = predict(cluster4, data = ftpdata2[,c(35,36,28,29)], newdata = ftpdata2[,c(35,36,28,29)])
prediction5 = predict(cluster5, data = ftpdata2[,c(35,36,28)], newdata = ftpdata2[,c(35,36,28)])
prediction6 = predict(cluster6, data = ftpdata2[,c(35,36,29)], newdata = ftpdata2[,c(35,36,29)])
prediction7 = predict(cluster7, data = ftpdata2[,c(35,36)], newdata = ftpdata2[,c(35,36)])
prediction8 = predict(cluster8, data = ftpdata2[,c(31,32)], newdata = ftpdata2[,c(31,32)])
prediction9 = predict(cluster9, data = ftpdata2[,c(25:26)], newdata = ftpdata2[,c(25:26)])
prediction10 = predict(cluster10, data = ftpdata2[,c(28,29,31,32,35,36)], newdata = ftpdata2[,c(28,29,31,32,35,36)])
prediction11 = predict(cluster11, data = ftpdata2[,c(33,34,35,36,28,29)], newdata = ftpdata2[,c(33,34,35,36,28,29)])
prediction12 = predict(cluster11, data = ftpdata2[,c(36,28,29)], newdata = ftpdata2[,c(36,28,29)])

unique(prediction1)
unique(prediction2)
unique(prediction3)
unique(prediction4)
unique(prediction5)
unique(prediction6)
unique(prediction7)
unique(prediction8)
unique(prediction9)
unique(prediction10)
unique(prediction11)
unique(prediction12)

ftpdata2$prediction1=prediction1
ftpdata2$prediction2=prediction2
ftpdata2$prediction3=prediction3
ftpdata2$prediction4=prediction4
ftpdata2$prediction5=prediction5
ftpdata2$prediction6=prediction6
ftpdata2$prediction7=prediction7
ftpdata2$prediction8=prediction8
ftpdata2$prediction9=prediction9
ftpdata2$prediction10=prediction10
ftpdata2$prediction11=prediction11
ftpdata2$prediction12=prediction12

clusterDF=ftpdata2

write.csv(clusterDF,"clusterDF.csv")

library(ggplot2)
library(plotly)

speed_clusters=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction4), y= speed_ms_scaled, color = as.factor(prediction4)))+ theme_classic()

turning_clusters=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction4), y= turnAngle_scaled,  color = as.factor(prediction4)))+ theme_classic()

fpt20_clusters=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction4), y= FPT_above_20, color = as.factor(prediction4)))+ theme_classic()

fpt50_clusters=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction4), y= FPT_above_50, color = as.factor(prediction4)))+ theme_classic()

speed_clusters2=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction12), y= speed_ms_scaled, color = as.factor(prediction12)))+ theme_classic()

fpt20_clusters2=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction12), y= FPT_above_20, color = as.factor(prediction12)))+ theme_classic()

fpt50_clusters2=ggplot(data=clusterDF)+ geom_boxplot(aes(x=as.factor(prediction12), y= FPT_above_50, color = as.factor(prediction12)))+ theme_classic()




ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction1)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .9, minPts = 10, eta + tau + fpt20 + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction2)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .9, minPts = 10, eta + tau + fpt20")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction3)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .9, minPts = 10, eta + tau + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction4)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .6, minPts = 10, speed + turning + fpt20 + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction5)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .6, minPts = 10, speed + turning + fpt20")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction6)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .6, minPts = 10, speed + turning + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction7)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .6, minPts = 10, speed + turning")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction8)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = .9, minPts = 10, eta + tau")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction9)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = .9, minPts = 10, fpt20 + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction10)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = 1, minPts = 10, speed, turning, eta, tau, fpt20 + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction11)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = 1, minPts = 10, speed, turning, eta, tau, fpt20 + fpt50")

ggplot(data=clusterDF)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction12)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = .6, minPts = 10, speed, fpt20 + fpt50")


###models 1,4 and 10 are winning from visual inspection 

Estrella = clusterDF[which(as.character(clusterDF$trackId)=="Estrella"),]
library(ggpubr)

Estrella_model1=ggplot(data=Estrella)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction1)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .9, minPts = 10, eta + tau + fpt20 + fpt50")

Estrella_model4= ggplot(data=Estrella)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction4)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId) + ggtitle("eps = .6, minPts = 10, speed + turning + fpt20 + fpt50")

Estrella_model10 = ggplot(data=Estrella)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction10)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = 1, minPts = 10, speed, turning, eta, tau, fpt20 + fpt50")

Estrella_model11 = ggplot(data=Estrella)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction11)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = 1, minPts = 10, speed, turning, eta, tau, fpt20 + fpt50 fpt20bin + fpt50bin")

Estrella_model12 = ggplot(data=Estrella)+ geom_point(aes(x=utm.easting, y= utm.northing, group = 1, color = as.factor(prediction12)))+ theme_classic()+ coord_equal()+ facet_wrap(~trackId)  + ggtitle("eps = .6, minPts = 10, speed, fpt20 + fpt50")

ggarrange(Estrella_model4, speed_clusters, common.legend = TRUE, nrow=2)
ggarrange(Estrella_model4,  turning_clusters, common.legend = TRUE, nrow=2)
ggarrange(Estrella_model4, fpt20_clusters, common.legend = TRUE, nrow=2)
ggarrange(Estrella_model4, fpt50_clusters, common.legend = TRUE, nrow=2)

ggarrange(Estrella_model4, Estrella_model12, common.legend = TRUE)


plot_ly(x=clusterDF$turnAngle_scaled, y=clusterDF$speed_ms_scaled, z=clusterDF$FPT_above_20, type="scatter3d", mode="markers", color=as.factor(clusterDF$prediction4))
