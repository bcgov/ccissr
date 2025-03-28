library(tidyverse)
#read in current table
suit<-read.csv("tables/versioned/suitability_v13_14.csv")
suit$X<-NULL
names(suit)

#read in update table
TO_update<-read.csv("tables/regional_updates/Thompson_Okanagan_Mar2025.csv") 
names(TO_update)

#pull out ratings of bgcs in TO update from main to compare  
suitx<-subset(suit, bgc %in% TO_update$bgc)

#make list of site series Mike has added suitability for (not currently in main table) 
new<-anti_join(TO_update, suitx,  by="ss_nospace") #added ratings 

#calculate change stats  
suitx<-rename(suitx, newsuit1=newsuit, mod1=mod)%>%select(-outrange)
suitx<-left_join(suitx, TO_update, by=c('bgc', 'ss_nospace', 'spp'))

#determine which ratings are not filled in (forested units only) that were previously 
check<-subset(suitx, is.na(newsuit))%>%separate(ss_nospace, into=c("zone", "ss"), sep='/', remove=F)%>%
  filter(. , !grepl("Fm|Gg|Gs|Ro|Rt|Ws", ss))%>%select(-ss, -zone)%>%
  subset( bgc!="ESSFdh2"& ss_nospace!="ESSFdc3/04"& spp!="X")
write.csv(check, 'checkTOmissing.csv')

#ESSFdc3 does not have an /04


#calculate change stats  
TO_update<-mutate(TO_update, feas_change= if_else(newfeas==FEAS2025, 0, 1))
nchange<-sum(TO_update$feas_change, na.rm = T)

#subset changed ratings 
updated_ratings<-subset(coast_update, feas_change>0)

#how many of these were initially blank?
blanks<-subset(updated_ratings, newfeas=="")
nchange=nchange-nrow(blanks)

#percent change 
pchange<-nchange/nrow(coast_update)*100 #ratings were changed/updated
padd<-nrow(blanks)/nrow(coast_update)*100 #ratings were added to site series where there were none

#bind updates to full feas table 
coast_update$newfeas<-NULL
coast_update$feas_change<-NULL
coast_update<-rename(coast_update, newfeas=FEASFEB2025)%>%relocate(newfeas,.after = spp)

names(coast_update)
names(feas)

#suit<-subset(suit, !bgc %in% TO_update$bgc)


#bind back together
feas<-rbind(feas, coast_update)

#write as new version
write.csv(feas, "feas_tables/versioned/Feasibility_v13_6.csv")
