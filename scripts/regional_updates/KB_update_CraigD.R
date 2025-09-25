library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_21.csv")
suit$X<-NULL
names(suit)

#read in update 
SE_update<-read.csv("tables/regional_updates/needsreviewSE2.csv") 
names(SE_update)

#calculate change stats  
SE_update<-mutate(SE_update, suit_change= if_else(newsuit==suitability, 0, 1))
nchange<-sum(SE_update$suit_change, na.rm = T)
#percent change 
pchange<-nchange/nrow(SE_update)*100 #ratings were changed/updated
#add to Deb M.'s ratings 
#16.3 % of reviewed ratings were updated (118/724)
ntot<-724+nrow(SE_update)
nchange<-118+nchange
pchange<-nchange/ntot*100 #ratings were changed/updated

#select relevant columns and rename for join
SE_update<-select(SE_update, -Zone, -Note, -Comments, -suit_change)

bgcs<-unique(SE_update$bgc)
suitSE<-subset(suit, bgc %in% bgcs)%>%subset(., is.na(mod)|mod=='inputed')%>%subset(., spp!="X")

#remove prev ratings from main dataframe
suit<-anti_join(suit, suitSE)

#add updated ratings back in 
suit<-rbind(suit, SE_update)

#write as new version
write.csv(suit, "tables/versioned/suitability_v13_22.csv")

