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

#fix missing "PPJWxh_CA" ratings 
suit<-read.csv("tables/versioned/Suitability_v13_28.csv")
suit$X<-NULL
names(suit)

USBEC_update<-read.csv("tables/regional_updates/Feasibility_USBEC_DVM_4Jan2026.csv") 
USBEC_update<-subset(USBEC_update, bgc=="PPJWxh_CA")%>%select(-Notes..Jan.Mar.2026, -X, -X.1, -X.2) 
names(USBEC_update)
USBEC_update<-rename(USBEC_update, suitability=feasible, newsuit=newfeas)

suit<-rbind(suit, USBEC_update)

#write out new version
write.csv(suit, "tables/versioned/Suitability_v13_29.csv")

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


#BEC Info----
#libraries 
#library(tidyverse)
library(dplyr)
library(stringr)
#load WNA BGCs list
WNA_BGC_list<-read.csv("tables/WNA_BGCs_Info.csv")
WNA_BGC_list$X<-NULL

#read in US 2026 info  
US2026<-read.csv("tables/regional_updates/USA_2026_BGCs.csv")
US2026<-subset(US2026, BGCName!="") #remove zone headers

#filter out of main BGC list 
US_BGC_list<-subset(WNA_BGC_list, BGC %in% US2026$BGC) #135 
US_BGC_list<-tidyr::separate(US_BGC_list, Subzone, into = c("Subzone"), sep = "_", remove = F) #clean up 

new_US_BGC_list<-subset(US2026, !BGC %in% WNA_BGC_list$BGC) #28 added in 2026
new_US_BGC_list$Source<- "USA_2026"

#overwrite names/info to 2026 naming convention 
US2026<-rename(US2026, ZoneName2=ZoneName, BGCName2=BGCName, General.Location2=General.Location, SubzoneName2=SubzoneName)
US_BGC_list<-left_join(US_BGC_list, US2026)

US_BGC_list$BGCName<-US_BGC_list$BGCName2
US_BGC_list$ZoneName<-US_BGC_list$ZoneName2
US_BGC_list$SubzoneName<-US_BGC_list$SubzoneName2
US_BGC_list$General.Location<-US_BGC_list$General.Location2

US_BGC_list<-select(US_BGC_list, -General.Location2, -ZoneName2, -BGCName2, -SubzoneName2)
names(US_BGC_list)
names(new_US_BGC_list)

new_US_BGC_list$SubzoneGroup<- ""
new_US_BGC_list$OldSubzoneLabel<- ""
new_US_BGC_list$Notes<-""
new_US_BGC_list$DataSet<-"USA"
new_US_BGC_list$VariantName<- ""

new_US_BGC_list<-tidyr::separate(new_US_BGC_list, BGC, into =  c("X", "State"), sep = "_", remove = F)
new_US_BGC_list$X<-NULL

new_US_BGC_list$SubzoneGroup[new_US_BGC_list$Zone=='ESSF'& grepl('Dry', new_US_BGC_list$SubzoneName)] <-  'ESSF_dry'
new_US_BGC_list$SubzoneGroup[new_US_BGC_list$Zone=='ESSF'& grepl('Moist', new_US_BGC_list$SubzoneName)] <-  'ESSF_meso'
new_US_BGC_list$SubzoneGroup[new_US_BGC_list$Zone=='ESSF'& grepl('Wet', new_US_BGC_list$SubzoneName)] <-  'ESSF_wet'
new_US_BGC_list$SubzoneGroup[new_US_BGC_list$Zone=='MH'& grepl('Submaritime', new_US_BGC_list$SubzoneName)] <-  'MH_submaritme'
new_US_BGC_list$SubzoneGroup[new_US_BGC_list$Zone=='MH'& grepl('Maritime', new_US_BGC_list$SubzoneName)] <-  'MH_maritme'

new_US_BGC_list<-mutate(new_US_BGC_list, SubzoneGroup=if_else(SubzoneGroup=="", Zone, SubzoneGroup))
new_US_BGC_list<-mutate(new_US_BGC_list, VariantName=if_else(VariantName=="", State, VariantName))

new_US_BGC_list<- select(new_US_BGC_list, Source,Zone, SubzoneGroup, Subzone,BGC,    BGCName,ZoneName,        
                         SubzoneName, VariantName, General.Location, OldSubzoneLabel,  DataSet,State,Notes)   

US_BGC_list<-rbind(US_BGC_list, new_US_BGC_list)

#remove US units from WNA_BGC list and replace with new un
WNA_BGC_list<-subset(WNA_BGC_list, !grepl('USA', Source))#remove them
WNA_BGC_list<-subset(WNA_BGC_list, BGC!="FGff" & BGC!="SBAPfp") #remove two from AB updated here 

WNA_BGC_list<-rbind(WNA_BGC_list, US_BGC_list)#add back in

write.csv(WNA_BGC_list, "tables/versioned/WNA_BGCs_v13_4.csv")

#add colors info
cols<-read.csv("tables/metadata/WNA_SubzoneCols.csv")
cols$BGC<-cols$classification
cols$classification<-NULL

cols$BGC[cols$BGC=="ESSFun"] <-  "ESSFun2"
cols$BGC[cols$BGC=="ICHun"] <-  "ICHun1"
cols$BGC[cols$BGC=="MHun"] <-  'MHun1'
cols$BGC[cols$BGC=="MSun"] <-  'MSun1'
cols$BGC[cols$BGC=="SWBun"] <-  'SWBun1'
cols$BGC[cols$BGC=="OWdm_OR"] <-  'OWwm_OR'
cols$BGC[cols$BGC=="OWun_CA"] <-  'OWxm_CA'
cols$BGC[cols$BGC=="JPWdm_WY"] <-  'PPJWdm_WY'
cols$BGC[cols$BGC=="JPWxw_NV"] <-  'PPJWxw_NV'
cols$BGC[cols$BGC=="JPWxh_CA"] <-  'PPJWxh_CA'
cols$BGC[cols$BGC=="JPWwm_CO"] <-  'PPJWwm_CO'
cols$BGC[cols$BGC=="JPWdw_UT"] <-  'PPJWdw_UT'
cols$BGC[cols$BGC=="JPWmk_WY"] <-  'PPJWmw_UT'
cols$BGC[cols$BGC=="ESSFmw"] <-  'ESSFwc7'

WNA_BGC_list<-left_join(WNA_BGC_list, cols)

#which are missing colors? 
miss<-subset(WNA_BGC_list, is.na(colour))
sort(unique(miss2$BGC))


cols<-unique(cols$classification)
sort(setdiff(BGC_list, cols))
sort(setdiff(cols, BGC_list))

