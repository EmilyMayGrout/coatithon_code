#this script is calculating the vedba for Galaxy group

#devtools::install_github("edwinth/thatssorandom")
##devtools::install_github("EliGurarie/smoove") 

library(move)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

#try dbscan with speed and turning angle and compare with eta/tau dbscan
#add vedba 
#parameter tuning for epsilon values and the min number of points 
#remove outliers (check for them)



#loading the acc data
#cred <- movebankLogin() # specially useful when sharing 
#Galaxy <- getMovebankData(study="Coati Galaxy Group Gamboa", login=cred)

#this reads in the acc for all individuals
#Galaxy_acc <- getMovebankNonLocationData(study="Coati Galaxy Group Gamboa" , sensorID="Acceleration", login=cred)

#write.csv(Galaxy_acc, "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/all_Galaxy_acc.csv")

Galaxy_acc2 <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/all_Galaxy_acc.csv")

Galaxy_acc2 <- Galaxy_acc2[,-1]
Galaxy_acc2$timestamp <- as.POSIXct(Galaxy_acc2$timestamp)

#--------------------------------------------------------

acc <- Galaxy_acc2
acc$individual_local_identifier <- as.factor(acc$individual_local_identifier)
acc_new <- split(acc, acc$individual_local_identifier)

#want to get pitch, roll, yaw and vedba from ACC
results = c()
count = 1
for(j in 1:length(acc_new)){
  temp <- acc_new[[j]]
  for(i in 1:nrow(temp)){
    if(i %% 1000 == 0){
      print(i)
    }
    burst <- temp$eobs_accelerations_raw[i]
    burstsplit <- strsplit(burst, split = ' ')
    xyz <- burstsplit[[1]]
    xindices <- seq(from = 1, to = length(xyz), by = 3)
    x <- xyz[xindices]
    x <- as.numeric(x)
    yindices <- seq(from = 2, to = length(xyz), by = 3)
    y <- xyz[yindices]
    y <- as.numeric(y)
    zindices <- seq(from = 3, to = length(xyz), by = 3)
    z <- xyz[zindices]
    z <- as.numeric(z)
    dx <- x - mean(x)
    dy <- y - mean(y)
    dz <- z - mean(z)
    vedba <- sqrt((dx^2)+(dy^2)+(dz^2))
    mean_vedba <- mean(vedba)
    temp$mean_vedba[i] <- mean_vedba
    #can add pitch, roll and yaw columns when I've worked out the issue with yaw
    temp$pitch[i] <- atan2(-mean(x), sqrt(mean(y)*mean(y) + mean(z)*mean(z)))
    temp$roll[i] <- atan2(mean(y), mean(z))
    temp$yaw[i] <- atan(mean(z)/sqrt(mean(x)*mean(x) + mean(z)*mean(z)))
    #saves each individuals data into a list called results
  }
  results[count] = list(temp)
  count = count + 1
}
results_merged <- do.call(rbind, results)
all_vedba <- results_merged

#write.csv(results_merged, "C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/all_Galaxy_vedba.csv", row.names = F)

acc_Galaxy <- read.csv("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/all_Galaxy_vedba.csv")

acc_Galaxy$individual_local_identifier <- as.factor(acc_Galaxy$individual_local_identifier)
acc_Galaxy$timestamp <- as.POSIXct(acc_Galaxy$timestamp)


#plot to each individual
acc_split <- split(acc_Galaxy, acc_Galaxy$individual_local_identifier)

#open the raw ACC from the decoder and merge the column with the ACC or ACCN so I know where the IMU acc is

# filenames <- list.files("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/decoder_acc_id/", pattern="*.txt", full.names=TRUE)
# ldf <- lapply(filenames, read.csv)
# name <- list.files("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/decoder_acc_id/", pattern="*.txt", full.names=FALSE)
# names(ldf) <- name


#open the decoder ACC and make into one dataframe
filenames <- dir("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/decoder_acc_id/", full.names = TRUE, pattern="*.txt") 
filenames2 <- dir("C:/Users/egrout/Dropbox/coatithon/rawdata/2022/galaxy/acc/decoder_acc_id/", full.names = FALSE, pattern="*.txt") 

decoder_acc <- read.table(filenames[1], sep=',')
decoder_acc$filename <-filenames2[1]
for( i in 2:length(filenames)){
  print(i)
  curr_data <- read.table(filenames[i], sep=',')
  curr_data$filename <- filenames2[i]
  decoder_acc <- rbind(decoder_acc, curr_data)
  
}




#coati_id=read.csv("C:/Users/egrout/Dropbox/stats_Franzi/data/coati_id.csv", header = F)








#histogram each individual

ggplot(data = acc_Galaxy, aes(x = log(mean_vedba)))+
geom_histogram(color="black", fill="white", bins = 50)+facet_wrap(~individual_local_identifier)+theme_classic()

#dbscan vedba
library(dbscan)

plot(acc_split$Cometa$pitch, acc_split$Cometa$roll)


cluster = dbscan(acc_Galaxy[,c(22, 23)], eps= .9, minPts = 50)
prediction = predict(cluster, data = acc_Galaxy[,c(22, 23)], newdata = acc_Galaxy[,c(22, 23)])
unique(prediction)
acc_Galaxy$prediction = prediction
clusterDF = acc_Galaxy

kmeanClust <- kmeans(acc_Galaxy$mean_vedba, centers=3)
acc_Galaxy$kmeans <- as.factor(kmeanClust$cluster)

ggplot(acc_Galaxy, aes(x = log(mean_vedba), fill = kmeans, color = kmeans)) + geom_histogram(alpha=0.5, bins = 50)+ theme_classic()+facet_wrap(~individual_local_identifier)




#for loading the resampled 10 min gps
#load(file = "C:/Users/egrout/Dropbox/stats_Franzi/data/resample_Galaxy.Rdata")


