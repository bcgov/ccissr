library(tidyverse)
#Suitability----
#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_25.csv")
suit$X<-NULL
names(suit)

#read in update 
USBEC_update<-read.csv("tables/regional_updates/Feasibility_USBEC_DVM_4Jan2026.csv") 
USBEC_update<-subset(USBEC_update, mod=="dvm26")%>%rename(mod2=mod)%>%select(-Notes..Jan.Mar.2026, -X, -X.1, -X.2) 
                                                                                 
#which ones are updates to ratings currently in table? 
suit<-left_join(suit, USBEC_update)
#update ratings/info in table
suit<-mutate(suit, newsuit=if_else(!is.na(newfeas), newfeas, newsuit))
suit<-mutate(suit, suitability=if_else(!is.na(feasible), feasible, suitability))
suit<-mutate(suit, mod=if_else(!is.na(mod2), mod2, mod))

suit2<-subset(suit, !is.na(newfeas))%>%select(-suitability, -newsuit, -mod2) #pull them out 

suit$newfeas<-NULL
suit$feasible<-NULL
suit$mod2<-NULL

#which are newly added? 
USBEC_update2<-anti_join(USBEC_update, suit2)
unique(USBEC_update$bgc)
USBEC_update2<-rename(USBEC_update2, mod=mod2, newsuit=newfeas, suitability=feasible)

#rbind to main 
names(USBEC_update2)
names(suit)
suit<-rbind(suit, USBEC_update2)

#change stats
#182 ratings updated 
#407 ratings added 

#write out new version
write.csv(suit, "tables/versioned/Suitability_v13_26.csv")


#Edatopic----
edat<-read.csv("tables/versioned/Edatopic_v13_11.csv") 

#read in update
Edat_update<-read.csv("tables/regional_updates/Edatopic_USBEC_DVM_4Jan2026.csv")
Edat_update<-subset(Edat_update, Source=="Add_2026"|Source=="USA_2026")

#Remove anything already in table for overwrite
check<-unique(Edat_update$SS_NoSpace)
edat<-subset(edat, !SS_NoSpace %in% check)

#bind back to table
edat<-rbind(edat, Edat_update)

#change stats
#72 edatopes updated  
#138 edatopes added

#write out updated table
write.csv(edat, "tables/versioned/Edatopic_v13_12.csv")

#add missing PPJW units 4/21/26
edat<-read.csv("tables/versioned/Edatopic_v13_12.csv") 
edat$X<-NULL
Edat_update<-read.csv("tables/regional_updates/Edatopic_USBEC_DVM_4Jan2026.csv")
Edat_update2<-subset(Edat_update, BGC=="PPJWdm_WY"|BGC=="PPJWdw_UT"|BGC=="PPJWwm_CO"|BGC=="PPJWxh_CA"|BGC=="PPJWxw_NV")
edat<-rbind(edat, Edat_update2)

write.csv(edat, "tables/versioned/Edatopic_v13_13.csv")


#remove old US units
edat<-read.csv("tables/versioned/Edatopic_v13_13.csv") 
edat$X<-NULL
old<-c("BGmk_ID","GOun_CO" ,"ICHvk_ID" ,"JPWdm_WY" ,"JPWdw_UT" ,"JPWmk_WY" ,"JPWwm_CO" ,"JPWxh_CA" ,"JPWxw_NV", "OWdm_OR" ,"OWun_CA" )
edat<-subset(edat, !BGC %in% old)
write.csv(edat, "tables/versioned/Edatopic_v13_14.csv")
