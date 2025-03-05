library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Feasibility_v13_7.csv")
suit$X<-NULL
names(suit)
suit<-rename(suit, suitability=feasible, newsuit=newfeas)#rename to suitability

#read in updated table
omineca_update<-read.csv("tables/regional_updates/omineca_cariboo_Feb2025.csv") 
names(omineca_update)
unique(omineca_update$expert_initials)
#subset to Daniel Sklar
omineca_update<-subset(omineca_update,  expert_initials=="DS")

#select relevant columns and rename for join
omineca_update<-select(omineca_update, BGC, SS_NoSpace, expert, spp, expert_updated, expert_initials)%>%rename(bgc=BGC, ss_nospace=SS_NoSpace)

#join 
suit2<-left_join(suit, omineca_update)

#calculate change stats  
suit2<-mutate(suit2, suit_change= if_else(newsuit==expert_updated, 0, 1))
nchange<-sum(suit2$suit_change, na.rm = T)

#subset changed ratings 
updated_ratings<-subset(suit2, suit_change>0)

#how many of these were initially blank?
blanks<-subset(updated_ratings, newsuit=="")
#nchange=nchange-nrow(blanks)

#percent change 
pchange<-nchange/nrow(omineca_update)*100 #ratings were changed/updated
#padd<-nrow(blanks)/nrow(omineca_update)*100 #ratings were added to site series where there were none


#replace with updated ratings & expert initials
suit3<-mutate(suit2, newsuit= if_else(!is.na(expert_initials), expert_updated, newsuit))
suit3<-mutate(suit3, mod= if_else(!is.na(expert_initials), expert_initials, mod))

#remove additional cols
suit3<-select(suit3, -expert, -expert_updated, -expert_initials, -suit_change)
names(suit3)

#write as new version
write.csv(suit3, "tables/versioned/suitability_v13_8.csv")
