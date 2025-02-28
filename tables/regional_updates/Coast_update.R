library(tidyverse)
#read in current table
feas<-read.csv("tables/versioned/Feasibility_v13_5.csv")
feas$X<-NULL
names(feas)

#read in updated table
coast_update<-read.csv("tables/regional_updates/Coast_Feb2025.csv") 
#coast_updatex<-subset(coast_update, ss_nospace %in% feas$ss_nospace) #already there -update suit ratings
#new_ss<-anti_join(coast_update, coast_updatex) #add new 

#delete bgcs from coast update in main 
feas<-subset(feas, !bgc %in% coast_update$bgc)
#delete CWHms1 & 2-> replaced with CWHms 4 & 5
feas<-filter(feas, bgc!="CWHms1" & bgc!="CWHms2")

check<-as.data.frame(unique(feas$bgc)) #look good 


#calculate change stats  
coast_update<-mutate(coast_update, feas_change= if_else(newfeas==FEASFEB2025, 0, 1))
nchange<-sum(coast_update$feas_change, na.rm = T)

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

#bind back together
feas<-rbind(feas, coast_update)

#write as new version
write.csv(feas, "feas_tables/versioned/Feasibility_v13_6.csv")
