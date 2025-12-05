library(tidyverse)
#read in current table version
suit<-read.csv("tables/versioned/suitability_v13_16.csv")
suit$X<-NULL
names(suit)

#clean up moderator/reviewer initials 
unique(suit$mod)
cts<-group_by(suit, mod)%>%summarise(cts=n())
suit<-mutate(suit, mod= case_when(mod=="BG"|mod=="BGx"~"BG",#Bob Green 
                                  mod=="AB"|mod=="Abx"|mod=="AB?"| mod=="ABx"~"AB", #Allen Banner 
                                  mod=="MRyan"|mod=="MR"~"MR", #Mike Ryan 
                                  mod=="KEI"|mod=="KI"~"KI",#Kristi Iversen
                                  mod=="WHM"~"WHM", #Will MacKenzie 
                                  mod=="DVM"~"DVM", #Del Meideinger 
                                  mod=="VC"|mod=="REVIEW"~"VC", #Vanessa Comeau
                                  #ratings added based on BEC guides for sites that did not have a suit rating
                                  mod=="DS"~"DS",#Dan Sklar
                                  mod=="KSD"~"KSD", #who is this? 
                                  mod=="SAS-HAK"|mod=="Coast"~"SAS-HAK", #Sari Saunders
                                  mod=="HK" ~"HK", #Heather Klassen
                                  TRUE~NA))
unique(suit$mod)
group_by(suit, mod)%>%summarise(cts=n())

#check that anything without a moderator isn't changed
check<-subset(suit, is.na(mod))%>%mutate(same=if_else(suitability==newsuit, T, F))
noinit<-subset(check, same=="FALSE")  #only 41 obs 

#input initials  
suitx<-left_join(suit, select(noinit, ss_nospace, sppsplit, spp, same))
suitx<-anti_join(suitx, noinit)

noinit<-mutate(noinit, mod=ifelse(grepl("CW", bgc), "DVM", "MR"))# Coast=Del, Interior=Mike 
suitx<-rbind(suitx, noinit)
suit<-select(suitx, -same)

check2<-subset(suit, is.na(mod))%>%mutate(same=if_else(suitability==newsuit, T, F))
noinit2<-subset(check2, same=="FALSE")  #0

write.csv(suit, "tables/versioned/suitability_v13_17.csv")


#read in current table 5/20/25
suit<-read.csv("tables/versioned/Suitability_v13_19.csv")
suit$X<-NULL
names(suit)

#coast update - update NAs
coast_update<-read.csv("tables/regional_updates/Coast_Feb2025.csv") 
coast_update<-subset(coast_update, !is.na(FEASFEB2025))%>%select(ss_nospace, spp)
coast_update$mod2<-"SAS-HAK"

suit<-left_join(suit, coast_update)

suit<-mutate(suit, mod=if_else(is.na(mod)&!is.na(mod2), mod2, mod))
suit$mod2<-NULL

#what has not been reviewed in BC? 
#notrev<-subset(suit, is.na(mod))
BCunits<-filter(suit, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|SBAP|SASbo', bgc))#US units
AB_units<-c("BWBScmC",  "BWBScmE",  "BWBScmNW" ,"BWBScmW",  "BWBSdmN" , "BWBSdmS" , "BWBSlbE" , "BWBSlbN",  "BWBSlbW", 
            "BWBSlf",   "BWBSnm",   "BWBSpp",   "BWBSub",   "BWBSuf")  
BCunits<-filter(BCunits, !bgc %in% AB_units) #AB units
BCunits<-mutate(BCunits, Zone= case_when(grepl('ICH', bgc)~"ICH",grepl('SWB', bgc)~"SWB", grepl('MS', bgc)~"MS",
                                         grepl('SBPS', bgc)~"SBPS", grepl('BAFA', bgc)~"BAFA", grepl('CWH', bgc)~"CWH", grepl('IDF', bgc)~"IDF",
                                         grepl('BG', bgc)~"BG", grepl('ESSF', bgc)~"ESSF", grepl('CDF', bgc)~"CDF", grepl('SBS', bgc)~"SBS", 
                                         grepl('MH', bgc)~"MH", grepl('CMA', bgc)~"CMA", grepl('PP', bgc)~"PP", grepl('BWBS', bgc)~"BWBS", TRUE~ NA))
BCunits<-subset(BCunits, !is.na(Zone)) #NAs
BCunits<-filter(BCunits, spp!="X") #no species  ~10K 
#BCunits<-filter(BCunits , !grepl("Fm|Gg|Gs|Ro|Rt|Ws", ss_nospace))#forested units only- inconsistently reviewed 
notrevBC<-subset(BCunits, is.na(mod)) #NAs


#Unresolved... 
sort(unique(notrevBC$bgc))
# "BWBSdk"  "BWBSmk"  "BWBSmw"  "BWBSvk"  "BWBSwk1" "BWBSwk3" "CWHvh3"  "CWHwh1"  "CWHwh2"  "ESSFdc2" "ESSFdcp"
# "ESSFdcw" "ESSFdk1" "ESSFdk2" "ESSFdkp" "ESSFdkw" "ESSFdvp" "ESSFdvw" "ESSFun"  "ESSFunp" "ESSFvc"  "ESSFvcp"
# "ESSFvcw" "ESSFwc4" "ESSFwcp" "ESSFwcw" "ESSFwh1" "ESSFwh2" "ESSFwh3" "ESSFwk1" "ESSFwm1" "ESSFwm2" "ESSFwm3"
# "ESSFwm4" "ESSFwmp" "ESSFwmw" "ESSFxc1" "ESSFxc3" "ESSFxcp" "ESSFxvp" "ESSFxvw" "ICHdk"   "ICHdw1"  "ICHdw3" 
# "ICHdw4"  "ICHmc1"  "ICHmc1a" "ICHmc2"  "ICHmk1"  "ICHmk3"  "ICHmk4"  "ICHmk5"  "ICHmm"   "ICHmw1"  "ICHmw2" 
# "ICHmw3"  "ICHmw4"  "ICHvc"   "ICHvk1"  "ICHvk2"  "ICHwc"   "ICHwk1"  "ICHwk2"  "ICHwk3"  "ICHxwa"  "IDFdh"  
# "IDFdk3"  "IDFdk4"  "IDFdk5"  "IDFdm1"  "IDFdm2"  "IDFxh1"  "IDFxh2"  "IDFxk"   "IDFxx1"  "IDFxx2"  "MSdc2"  
# "MSdc3"   "MSdk"    "MSdm1"   "MSdm2"   "MSdm3"   "MSdv"    "MSdw"    "MSxk1"   "MSxv"    "SBPSdc"  "SBPSmc" 
# "SBPSmk"  "SBPSxc"  "SBSdh1"  "SBSdh2"  "SBSdk"   "SBSdw2"  "SBSdw3"  "SBSmc2"  "SBSmh"   "SBSmk1"  "SBSmm"  
# "SBSun"   "SBSwk1"  "SBSwk3"  "SWBun"   "SWBvk"  


# "BWBSdk"  "BWBSmk"  "BWBSmw"  "BWBSvk"  "BWBSwk1" "BWBSwk3" #"ICHmm" "ICHwk3"  - Dan did not do some spp- ask Craig De Long?  
#"ICHmc1"  "ICHmc1a" "ICHmc2" "ESSFun" "ESSFunp "ICHvc" "ICHwc" - Erica did not do some species- ask Craig De Long? 

#"CWHvh3"  "CWHwh1"  "CWHwh2"- Haida Gwaii not reviewed by Coast - check with Sari & Heather
#"ESSFdvp" "ESSFdvw" "MSdc2"  "MSdc3"  "MSdv" - ask Coast?? 

#"ESSFdc2" "ESSFdcp"  "ESSFdcw" "ESSFxc1" "ESSFxc3" "ESSFxcp" "ICHdw3 "ICHmw3" "ICHvk1" "ICHvk2 "ICHwk1" "IDFdh"
#"IDFdm1" "IDFxh1"  "IDFxh2"  "IDFxx1"  "MSdm1"   "MSdm2"   "MSdm3" "MSxk1" - #Mike?
  
# - Deb still working on -"I need to finish the ICHmk variants plus the IDF and MS."
#"ESSFdk1" "ESSFdk2" "ESSFdkp" "ESSFdkw" "ESSFvc"  "ESSFvcp"
#"ESSFvcw" "ESSFwc4" "ESSFwcp" "ESSFwcw" "ESSFwh1" "ESSFwh2" "ESSFwh3" 
#"ESSFwk1" "ESSFwm1" "ESSFwm2" "ESSFwm3" "ESSFwm4" "ESSFwmp" "ESSFwmw"
# "ICHmk1"  "ICHmk3"  "ICHmk4"  "ICHmk5" "ICHdw1"  "ICHdw4"  "ICHmw1"  "ICHmw2" 
#"ICHmw4"  "ICHxwa"  "IDFdk5"  "IDFdm2" "IDFxk"    "IDFxx2" "MSdk"  "MSdw"   

#ESSFxvw, "ICHdk  "ICHwk2" "IDFdk3"  "IDFdk4" "MSxv"  - Kristi

unique(suit$mod)

#make summary plots---- 
BCunits<-mutate(BCunits, status=if_else(!is.na(mod), 'reviewed', 'not reviewed'))
BCunits<-mutate(BCunits, Status=case_when(suitability!=newsuit & status=='reviewed'~ 'reviewed-changed' ,
                                           suitability==newsuit & status=='reviewed'~'reviewed-not changed', 
                                           status=='not reviewed' ~'not reviewed', 
                                           TRUE~'reviewed-added'))

BCunits$ratings<-1


#separate conifers/broadleaves 
BCunits<-mutate(BCunits, treetype=if_else(spp=="Cw"|spp=="Fd"|spp=="Hw"|spp=="Yc"|spp=="Bl"|spp=="Tw"|spp=="Hm"|
                                          spp=="Ba"|spp=="Bg"| spp=="Lw"|spp=="Lt"|spp=="Jr"| spp=="Js"|spp=="Bp"|
                                            spp=="Pa"|spp=="Pf"|spp=="Pl"|spp=="Pw"|spp=="Py"|spp=="La"|spp=="Oa"|
                                            spp=="Sb"|spp=="Ss"|spp=="Sw"|spp=="Sx", "Conifer", "Broadleaf"))

cts<-group_by(BCunits, spp) %>%summarise(nrat=n())
BCunits<-mutate(BCunits, sppsplit=if_else(is.na(sppsplit), spp, sppsplit))#fill in NAs for subspp for broadleaf

#remove rare conifers from plot to make smaller 
BCunits<-subset(BCunits, spp!='Oa' &  spp!="Sw" & spp!="Jr" & spp!="Js" )
#combine BG & PP plots
BCunits<-mutate(BCunits, Zone2= ifelse(Zone=="BG"|Zone=="PP", "BG & PP", Zone))
# Stacked bar plot

library(paletteer)
ggplot(subset(BCunits,treetype=="Conifer"), aes(fill=Status, x=spp, y=ratings)) + 
  geom_bar(position="stack", stat="identity") + facet_wrap(~Zone2, scales="free_y", ncol=2)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top") +ggtitle("Conifers") +
  scale_fill_paletteer_d("wesanderson::AsteroidCity3") +ylab("N ratings")

ggplot(subset(BCunits,treetype=="Broadleaf"), aes(fill=Status, x=sppsplit, y=ratings)) + 
  geom_bar(position="stack", stat="identity") + facet_wrap(~Zone, scales="free_y", ncol=3)+ theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "top") +ggtitle("Broadleaves")+ ylab("N ratings")+
  scale_fill_paletteer_d("wesanderson::AsteroidCity3")


#August 2025
#Haida Gwaii units under review by SS/HK (CWHvh3, wh1, wh2)
notrevBC<-subset(notrevBC, Zone!="CWH"& bgc!="ESSFun" & bgc!="ESSFunp")

#add CGC imputed values to review list 
suit_add<-read.csv("tables/regional_updates/imputed_suit_ratings_CGC.csv")
suit_add$X<-NULL
suit_add<-mutate(suit_add, Zone= case_when(grepl('ICH', bgc)~"ICH",grepl('SWB', bgc)~"SWB", grepl('MS', bgc)~"MS",
                                         grepl('SBPS', bgc)~"SBPS", grepl('BAFA', bgc)~"BAFA", grepl('CWH', bgc)~"CWH", grepl('IDF', bgc)~"IDF",
                                         grepl('BG', bgc)~"BG", grepl('ESSF', bgc)~"ESSF", grepl('CDF', bgc)~"CDF", grepl('SBS', bgc)~"SBS", 
                                         grepl('MH', bgc)~"MH", grepl('CMA', bgc)~"CMA", grepl('PP', bgc)~"PP", grepl('BWBS', bgc)~"BWBS", TRUE~ NA))

notrevBC<-rbind(notrevBC, suit_add)
notrevBC<-mutate(notrevBC, Note=if_else(mod=="inputed", 'imputed', 'NA'))
notrevBC$mod<-NA
write.csv(notrevBC, "needsreview.csv")


#Sep 2025- 
# compare all outstanding reviews against needs review table 
coast_new<-read.csv("tables/regional_updates/Coast_review_missing_SCS_HAK_review.csv") 
coast_new$Note<-" "
SE_update<-read.csv("tables/regional_updates/needsreviewSE.csv") #Craig delong 9/19
SE_update$Comments<-NULL
NE_update<-read.csv("tables/regional_updates/needsreviewNE.csv") #Craig delong 9/12
NE_update$Comments<-NULL
NE_update$Footnotes<-NULL
all_new<-rbind(coast_new, NE_update, SE_update)
all_new$AugSep2025<-"Y"
all_new<-select(all_new, bgc, ss_nospace, sppsplit, spp, AugSep2025)

notrevBC<-read.csv("needsreview.csv")

notrevBC<-left_join(notrevBC, all_new)
notrevBC<-subset(notrevBC, is.na(AugSep2025))#remove all reviewed by Craig or Coast in Aug/Sep 2025
notrevBC$AugSep2025<-NULL

sort(unique(notrevBC$bgc))

#which are in SE/NE? 
notrevBCx<-subset(notrevBC, bgc=="ESSFdk2"|bgc=="ESSFdkp"|bgc=="ESSFwm1"|bgc=="ESSFmm3"|bgc=="ESSFmm1"|bgc=="ESSFwc3"|bgc=="ICHmk4"|
                    bgc=="ICHmk5"|bgc=="IDFdk5"|bgc=="IDFxk"|bgc=="IDFxx2"|bgc=="IMAunp"|bgc=="MSdw"|bgc=="SBSdh2")

notrevBCNE<-subset(notrevBCx, bgc=="ESSFwc3"|bgc=="SBSdh2"| bgc=="IMAunp"| bgc=="ESSFmm1")
notrevBCSE<-anti_join(notrevBCx, notrevBCNE)

#add back to review table 
SE_update<-rbind(SE_update, notrevBCSE)
write.csv(SE_update, "tables/regional_updates/needsreviewSE2.csv") 

write.csv(notrevBCNE, "tables/regional_updates/needsreview_misc.csv") 

#what is STILL not reviewed? 
notrevBCy<-anti_join(notrevBC, notrevBCx)
sort(unique(notrevBCy$bgc))

#add imputed ratings for coast to csv
notrevBCz<-subset(notrevBCy, Zone== "CDF"|Zone=="CMA"|Zone=="CWH"|Zone=="MH"|bgc=="ESSFmw"|bgc=="ESSFmw1"|bgc=="ESSFmwp"|
                    bgc=="ESSFun"|bgc=="IDFww"|bgc=="ESSFmkw")
coast_new<-rbind(coast_new, notrevBCz)

write.csv(coast_new, "tables/regional_updates/Coast_review_missing_SCS_HAK_review2.csv") 

#what is STILL not reviewed? Not sure who can look at these remaining units...  
#many are parkland, woodland and floodplain/wetland
#some inputed 
notrevBCy<-anti_join(notrevBCy, notrevBCz)
sort(unique(notrevBCy$bgc))
check<-group_by(notrevBCy, bgc)%>%summarise((ct=n()))

x<-read.csv("regions_tab.csv")
x<-distinct(x)
notrevBCy<-left_join(notrevBCy, x)

write.csv(notrevBCy, "needsreview_sub.csv")

#combine 
notrevBC<-read.csv("needsreview_sub.csv")
notrevBC2<-read.csv("tables/regional_updates/needsreview_misc.csv")
x<-read.csv("regions_tab.csv")
notrevBC2<-left_join(notrevBC2, x)
notrevBC<-rbind(notrevBC, notrevBC2)
notrevBC$X<-NULL

write.csv(notrevBC, "needsreview_sub.csv")

#Nov 2025- 
suit<-read.csv("tables/versioned/suitability_v13_23.csv") 
suit<-subset(suit, spp!="X")

#take out the US and alberta stuff because it won't match plot data
suit<-filter(suit, !grepl('_OC|_WC|_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|	MGPmg|
 MGPdm|SBAP|SASbo|BWBScmC|BWBScmE|BWBScmNW|BWBScmW|BWBSdmN|BWBSdmS|BWBSlbE|BWBSlbN|BWBSlbW|BWBSlf|BWBSnm|BWBSpp|BWBSub|BWBSuf', ss_nospace))

needs_rev<-subset(suit, is.na(mod)|mod=="inputed")

x<-read.csv("regions_tab.csv")
needs_rev<-left_join(needs_rev, x)
needs_rev<-mutate(needs_rev, Region= if_else(is.na(Region), "Coast", Region))

write.csv(needs_rev, "needs_review_Nov2025.csv")

