library(tidyverse)

#read in current table
suit<-read.csv("tables/versioned/Suitability_v13_19.csv")
suit$X<-NULL
names(suit)

#read in inputed ratings from BEC plot data (ss_quant repo)
suit_add<-read.csv("tables/regional_updates/inputed_suit_ratings_CGC.csv")
#select only suit data 
suit_add<-mutate(suit_add, sppsplit=case_when(spp=='Pl'& grepl("CWH|ESSFun|MH", bgc)~"Plc",
                                              spp=='Pl'& !grepl("CWH|ESSFun|MH", bgc)~"Pli", 
                                              spp=='Fd'& grepl("CWH|ESSFun|MH", bgc)~"Fdc",
                                              spp=='Fd'& !grepl("CWH|ESSFun|MH", bgc)~"Fdi", 
                                              TRUE~spp))
suit_add<-select(suit_add, bgc, ss_nospace, sppsplit, suitability, spp, newsuit,mod,outrange) %>% distinct(.)
write.csv(suit_add, "tables/regional_updates/inputed_suit_ratings_CGC2.csv")
#add to main table 
suit<-rbind(suit, suit_add)

#write as new version
write.csv(suit, "tables/versioned/suitability_v13_20.csv")
