library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_17.csv")
suit$X<-NULL
names(suit)

#read in update 
SkeenaNE_update<-read.csv("tables/regional_updates/Skeena&NE_review.csv") 
names(SkeenaNE_update)

#select relevant columns and rename for join
SkeenaNE_update<-select(SkeenaNE_update, BGC, SS_NoSpace, expert, spp, expert_updated, expert_initials)%>%
  rename(bgc=BGC, ss_nospace=SS_NoSpace)

#join 
suit2<-left_join(suit, SkeenaNE_update)

#how many of these are new? 
suit3<-subset(suit2, !is.na(expert_initials))
new<-anti_join(SkeenaNE_update, suit3)

#fill in not changed for Erica 
suit2<-mutate(suit2,expert_updated=if_else(expert_initials=="EBL"&is.na(expert_updated), expert, expert_updated))

#calculate change stats  
suit2<-mutate(suit2, suit_change= if_else(newsuit==expert_updated, 0, 1))
nchange<-sum(suit2$suit_change, na.rm = T)

#add in other 2 new SBS units from Erica
new2<-read.csv("tables/regional_updates/Skeena_SBSmz_vz_Feb2025.csv")

#percent change 
pchange<-nchange/(nrow(SkeenaNE_update)+ nrow(new2))*100 #ratings were changed/updated
padd<-(nrow(new)+ nrow(new2))/(nrow(SkeenaNE_update) + nrow(new2))*100 #ratings were added to site series where there were none

#replace with updated ratings & expert initials
suit<-mutate(suit2, newsuit= if_else(!is.na(expert_initials), expert_updated, newsuit))
suit<-mutate(suit, mod= if_else(!is.na(expert_initials), expert_initials, mod))

#remove additional cols
suit<-select(suit, -expert, -expert_updated, -expert_initials, -suit_change)
names(suit)

#add back new ratings 
new<-rename(new, mod=expert_initials, suitability=expert, newsuit=expert_updated)
new$sppsplit<-NA
new$outrange<-NA
new<-relocate(new, sppsplit, .before = suitability)
new2<-rename(new2, suitability=feasible, newsuit=newfeas)%>%select(-comments)
names(new)
names(new2)
new<-rbind(new, new2)

suit<-rbind(suit, new)
#writenew#write asnew#writenew#write as new version
write.csv(suit, "tables/versioned/suitability_v13_18.csv")
