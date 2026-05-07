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


#update naming to SWBun1
suit<-read.csv("tables/versioned/Suitability_v13_27.csv")
suit$X<-NULL
names(suit)

suit <- suit %>%
  mutate(bgc = str_replace_all(bgc, "SWBun", "SWBun1"))%>%
  mutate(ss_nospace = str_replace_all(ss_nospace, "SWBun", "SWBun1"))
write.csv(suit, "tables/versioned/suitability_v13_28.csv")

edat<-read.csv("tables/versioned/Edatopic_v13_15.csv") 
edat$X<-NULL
names(edat)

uns<-subset(edat,BGC=="SWBuns")
edat<-subset(edat,BGC!="SWBuns")#pull these out so not str replaced 

edat <- edat %>%
  mutate(BGC = str_replace_all(BGC, "SWBun", "SWBun1"))%>%
  mutate(SS_NoSpace = str_replace_all(SS_NoSpace, "SWBun", "SWBun1"))
edat<-rbind(edat, uns) #add back in 
write.csv(edat, "tables/versioned/Edatopic_v13_16.csv") 
