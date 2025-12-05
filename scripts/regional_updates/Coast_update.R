library(tidyverse)
#read in current table
feas<-read.csv("tables/versioned/Feasibility_v13_5.csv")
feas$X<-NULL
names(feas)
#first update Feb 2025----
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

#second update Haida Gwaii (pt 1) Aug 2025----
library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_22.csv")
suit$X<-NULL
names(suit)

#read in update 
Coast_update2<-read.csv("tables/regional_updates/Coast_review_missing_SCS_HAK_review.csv") 
Coast_update2<-subset(Coast_update2, !is.na(newsuit))
unique(Coast_update2$bgc)
Coast_update2$Zone<-NULL

#remove updated ratings from main dataset
suit<-mutate(suit, remove= case_when(bgc=="CWHvh3" &  !is.na(suitability)~1, 
                                     bgc=="MHwh"&  !is.na(suitability)~1, 
                                     bgc=="MHwhp"&  !is.na(suitability)~1, 
                                     TRUE~0)) %>%subset(., remove<1)%>%select(-remove)

#bind back in 
suit<-rbind(suit, Coast_update2)

#calculate change stats  
Coast_update2<-mutate(Coast_update2, change= if_else(suitability==newsuit, 0, 1))
nchange<-sum(Coast_update2$change, na.rm = T) #29
#how many of these were initially blank?
nadd<-sum(is.na(Coast_update2$change)) #54 added 

#add values to what is in TR 
#192/3027->221/3230->6.8
#91/3027->145/3230->4.5

#fill in other blanks in mod column from February review 
coast<-read.csv("tables/regional_updates/Coast_Feb2025.csv") 
suit<-mutate(suit, fill= case_when(is.na(mod) & bgc %in% coast$bgc~ 'SAS-HAK', TRUE~ NA))%>%mutate(mod=if_else(is.na(mod), fill, mod))
suit$fill<-NULL
  
#write out new version
suit<-write.csv(suit, "tables/versioned/Suitability_v13_23.csv")


#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_23.csv")
suit$X<-NULL
names(suit)

#take out the US and alberta stuff because it won't match plot data
nonBC<-filter(suit, grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))

suit<-anti_join(suit, nonBC)

#remove units that were in needs_rev Nov 2025
needs_rev<-subset(suit, is.na(mod)|mod=="inputed")
suit<-anti_join(suit, needs_rev)

#Bring in last reviewed units Dec 4, 2025
#Mostly coast but also Cariboo
updates<-read.csv("needs_review_Nov2025.csv")
updates<-mutate(updates, mod=if_else(Notes=='imputed'& mod=="", Notes, mod))
updates<-subset(updates, mod!='imputed')#if imputed values not reviewed, remove from main df

updates<-mutate(updates, newsuit= if_else(is.na(newsuit), suitability, newsuit))
updates<-mutate(updates, outrange= ifelse(is.na(outrange), "FALSE", outrange))
updates<-select(updates, -Notes, -Region)

updates<-mutate(updates, mod = na_if(mod, " "))
suit<-rbind(suit, updates, nonBC)
suit<-mutate(suit, newsuit= if_else(is.na(newsuit), suitability, newsuit))
#fill in missing spp codes 
suit<-mutate(suit, spp = na_if(spp, ""))
filter(suit, is.na(spp))
suit<-mutate(suit, spp= if_else(is.na(spp), sppsplit, spp))

#remove retired units from BEC13 CMAunp 
suit<-filter(suit, bgc!="CMAunp")

#write out version
write.csv(suit, "tables/versioned/Suitability_v13_24.csv")


