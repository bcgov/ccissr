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

#merge remaining 
suitx<-rename(suitx, newsuit1=newsuit, mod1=mod)%>%select(-outrange)
suitx<-left_join(suitx, TO_update, by=c('bgc', 'ss_nospace', 'spp'))

#determine which ratings are not filled in that were previously in the table 
check<-subset(suitx, is.na(newsuit))%>%separate(ss_nospace, into=c("zone", "ss"), sep='/', remove=F)%>%
  filter(. , !grepl("Fm|Gg|Gs|Ro|Rt|Ws", ss))%>%select(-ss, -zone)%>% #forested units only
  subset( bgc!="ESSFdh2"& ss_nospace!="ESSFdc3/04"& spp!="X") #ESSFdc3 does not have an /04
#write.csv(check, 'tables/regional_updates/checkTOmissing.csv')


#read in Mike's responses 4/9/25
check2<-read.csv("tables/regional_updates/checkTOmissing_MR.csv")
suitx<-left_join(suitx, check2)

#remove ones Mike says to delete 
suitx<-subset(suitx, Action!="delete"|is.na(Action))

#update ratings for some 
suitx<- mutate(suitx, 
  action_numeric = suppressWarnings(as.numeric(Action)),  # non-numeric → NA
  newsuit = if_else(is.na(newsuit) & !is.na(action_numeric), action_numeric, newsuit)) %>%
  select(-action_numeric) 

#calculate change stats  
suitx<-mutate(suitx, feas_change= if_else(newsuit==newsuit1, 0, 1))%>%
  mutate(feas_change=if_else(is.na(feas_change), 0, feas_change))

#percent change 
nchange<-sum(suitx$feas_change == 1, na.rm = TRUE)
nsame<-sum(suitx$feas_change == 0, na.rm = TRUE)
total<-nrow(suitx)-nblank
nrow(new)

round(nchange/total, 2)*100 
round(nrow(new)/total, 2)*100

blanks<-subset(suitx, is.na(newsuit)& is.na(Action))
sort(unique(blanks$ss_nospace)) #check remaining blanks-> mostly non forested units 

#clean up 
suitxy<-select(suitx,bgc, ss_nospace, spp, newsuit, mod, outrange)
TO_update2<-rbind(suitxy, new)%>%mutate(outrange=if_else(spp=="Lw", TRUE, FALSE))

#Mike provided ratings for 01, 02 etc. ESSFdh2
TO_update3<-filter(TO_update2, bgc=='ESSFdh2' & is.na(mod))
TO_update2<-anti_join(TO_update2, TO_update3)  
TO_update2<-subset(TO_update2, ss_nospace!="ESSFdc3/04")

#keep previous non-forested unit ratings 
blanks<-select(blanks, bgc, ss_nospace, spp, newsuit1, mod1)
TO_update4<-left_join(TO_update2, blanks)
TO_update2<-mutate(TO_update4, newsuit=if_else(is.na(newsuit),newsuit1, newsuit))%>%
  select(-newsuit1, -mod1)

#pull out ratings of bgcs in TO update from main  
suitx<-subset(suit, bgc %in% TO_update$bgc)
suit<-anti_join(suit, suitx)

#add back in original suitability
suitx<-select(suitx, bgc, ss_nospace, sppsplit, suitability, spp)
TO_update2<-left_join(TO_update2, suitx)

#merge back to main table 
suit<-rbind(suit, TO_update2)

#write as new version
write.csv(suit, "tables/versioned/Suitability_v13_15.csv")


