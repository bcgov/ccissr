library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/suitability_v13_12.csv")
suit$X<-NULL
names(suit)

#read in update 
SWBun_update<-read.csv("tables/regional_updates/SWBun_updatedSS.csv") 

#remove SWBun old SS
SWBun<-subset(suit, bgc=="SWBun")#45 obs
suit<-subset(suit, bgc!="SWBun")

#bind in new units x/m/h & ratings to main suitability table 
suit<-rbind(suit, SWBun_update)

#write as new version
write.csv(suit, "tables/versioned/suitability_v13_13.csv")


