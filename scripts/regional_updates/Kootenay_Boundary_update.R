library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_18.csv")
suit$X<-NULL
names(suit)

#read in update 
KB_update<-read.csv("tables/regional_updates/KB_review.csv") 
names(KB_update)

#clean up initials column
KB_update<-separate(KB_update, expert_initials, into = c("expert_initials", "Notes"), sep = "DJM")%>%
  mutate(Notes2=if_else(!is.na(expert_initials), expert_initials, Notes))%>% unite(Notes, c("Notes", "Notes2"), sep = "", na.rm = T)%>%
  mutate(expert_initials="DJM")

#clean up ratings column
KB_update<-separate(KB_update, expert_updated, into = c("expert_updated", "Notes2"), sep = "(?<=[0-9]) (?=[A-Za-z])")%>%
  unite(Notes, c("Notes", "Notes2"), sep = "", na.rm = T)
KB_update$expert_updated<-as.numeric(KB_update$expert_updated)

#fill in not changed 
KB_update<-mutate(KB_update,expert_updated=if_else(is.na(expert_updated), expert, expert_updated))

#select relevant columns and rename for join
KB_update<-select(KB_update, BGC, SS_NoSpace, expert, spp, expert_updated, expert_initials)%>%
  rename(bgc=BGC, ss_nospace=SS_NoSpace)

#join 
suit2<-left_join(suit, KB_update)

#calculate change stats  
suit2<-mutate(suit2, suit_change= if_else(newsuit==expert_updated, 0, 1))
nchange<-sum(suit2$suit_change, na.rm = T)

#percent change 
pchange<-nchange/nrow(KB_update)*100 #ratings were changed/updated

#replace with updated ratings & expert initials
suit<-mutate(suit2, newsuit= if_else(!is.na(expert_initials), expert_updated, newsuit))
suit<-mutate(suit, mod= if_else(!is.na(expert_initials), expert_initials, mod))

#remove additional cols
suit<-select(suit, -expert, -expert_updated, -expert_initials, -suit_change)
names(suit)

#write as new version
write.csv(suit, "tables/versioned/suitability_v13_19.csv")
