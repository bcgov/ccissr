library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/suitability_v13_8.csv")
suit$X<-NULL
names(suit)

#read in update 
cariboo_update<-read.csv("tables/regional_updates/omineca_cariboo_Feb2025.csv") 
unique(cariboo_update$expert_initials)
#subset to Kristi Iversen 
cariboo_update<-subset(cariboo_update,  expert_initials=="KI"| expert_initials=="KEI")

#select relevant columns and rename for join
cariboo_update<-select(cariboo_update, BGC, SS_NoSpace, expert, spp, expert_updated, expert_initials)%>%rename(bgc=BGC, ss_nospace=SS_NoSpace)

#join 
suit2<-left_join(suit, cariboo_update)

#calculate change stats  
suit2<-mutate(suit2, suit_change= if_else(newsuit==expert_updated, 0, 1))
nchange<-sum(suit2$suit_change, na.rm = T)

#subset changed ratings 
updated_ratings<-subset(suit2, suit_change>0)

#how many of these were initially blank?
blanks<-subset(updated_ratings, newsuit=="")
#nchange=nchange-nrow(blanks)

#percent change 
pchange<-nchange/nrow(cariboo_update)*100 #ratings were changed/updated
#padd<-nrow(blanks)/nrow(cariboo_update)*100 #ratings were added to site series where there were none


#replace with updated ratings & expert initials
suit3<-mutate(suit2, newsuit= if_else(!is.na(expert_initials), expert_updated, newsuit))
suit3<-mutate(suit3, mod= if_else(!is.na(expert_initials), expert_initials, mod))

#remove additional cols
suit3<-select(suit3, -expert, -expert_updated, -expert_initials, -suit_change)
names(suit3)

#write as new version
write.csv(suit3, "tables/versioned/suitability_v13_9.csv")


# add back in new units - 3/11/25 
suit<-read.csv("tables/versioned/suitability_v13_11.csv")
suit$X<-NULL

#read in updated table
cariboo_update<-read.csv("tables/regional_updates/omineca_cariboo_Feb2025.csv") 
#subset to Kristi Iversen 
cariboo_update<-subset(cariboo_update,  expert_initials=="KI"| expert_initials=="KEI")

#select relevant columns and rename for join
cariboo_update<-select(cariboo_update, BGC, SS_NoSpace, expert, spp, expert_updated, expert_initials)%>%rename(bgc=BGC, ss_nospace=SS_NoSpace)

#join 
suit3<-left_join(suit, cariboo_update)

#figure out which were missed
cariboo_check<-subset(suit3, expert_initials=="KI"| expert_initials=="KEI")
diff<-anti_join(cariboo_update, cariboo_check)

#remove duplicate cols
suit3<-select(suit3, -expert, -expert_updated, -expert_initials)
names(suit3)
names(diff)

#rename or add cols to match 
diff<-rename(diff, mod=expert_initials, suitability=expert, newsuit=expert_updated)
diff$outrange<-"FALSE"
unique(diff$spp)
diff<-mutate(diff, sppsplit= case_when(spp=="Fd"~"Fdi", 
                                       spp=="Py"~"Pyi",
                                       spp=="Pl"~"Pli", 
                                       TRUE~ NA))
diff<-mutate(diff, sppsplit=if_else(is.na(sppsplit), spp, sppsplit))

#reorder cols
diff<-select(diff, bgc, ss_nospace, sppsplit, suitability, spp, newsuit, mod, outrange)

#bind back to main suitability table
suit3<-rbind(suit3, diff)

#write as new version
write.csv(suit3, "tables/versioned/suitability_v13_12.csv")


