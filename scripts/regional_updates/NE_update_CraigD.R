library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_20.csv")
suit$X<-NULL
names(suit)

#read in update 
NE_update<-read.csv("tables/regional_updates/needsreviewNE.csv") 
names(NE_update)

#calculate change stats  
NE_update<-mutate(NE_update, suit_change= if_else(newsuit==suitability, 0, 1))
nchange<-sum(NE_update$suit_change, na.rm = T)
#percent change 
pchange<-nchange/nrow(NE_update)*100 #ratings were changed/updated

#select relevant columns and rename for join
#need to save footnotes somewhere else to be input into CCISS. Asked Kiri 9/10/25
NE_update<-select(NE_update, -Zone, -Note, -Comments, -Footnotes, -suit_change)

bgcs<-unique(NE_update$bgc)
suitNE<-subset(suit, bgc %in% bgcs)%>%subset(., is.na(mod)|mod=='inputed')%>%subset(., spp!="X")

#remove prev ratings from main dataframe
suit<-anti_join(suit, suitNE)

#add updated ratings back in 
suit<-rbind(suit, NE_update)

#write as new version
write.csv(suit, "tables/versioned/suitability_v13_21.csv")

