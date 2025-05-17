library(tidyverse)
#read in current table
suit<-read.csv("tables/suitability.csv")
suit$X<-NULL
sort(unique(suit$spp))

broads<-subset(suit, spp=="At"|spp=="Act"|spp=="Ac"|spp=="Ep"|spp=="Dr"|spp=="Mb"|
                 spp=="Ra"|spp=="Wt"|spp=="Mv"|spp=="Ws"|spp=="Wp"|spp=="Qg"|spp=="Gp"|spp=="Vb"|
                 spp=="Up"|spp=="Kc"| spp=="Vp")
#remove US units 
broads<-filter(broads, !grepl('_CA|_OR|_WA|_ID|_MT|_CA|_WY|_CO|_NV|UT|BSJP|abE|abN|abS|abC|SBAP|SASbo', ss_nospace))

broads<-mutate(broads, zone=(case_when(grepl("SWB", bgc)~"SWB",
                                 grepl("BWBS", bgc)~"BWBS",
                                 grepl("SBPS", bgc)~"SBPS", 
                                 grepl("SBS", bgc)~"SBS",
                                 grepl("MH", bgc)~"MH",
                                 grepl("ESSF", bgc)~"ESSF",
                                 grepl("MS", bgc)~"MS",
                                 grepl("BG", bgc)~"BG",
                                 grepl("PP", bgc)~"PP",
                                 grepl("IDF", bgc)~"IDF",
                                 grepl("CDF", bgc)~"CDF",
                                 grepl("ICH", bgc)~"ICH",
                                 grepl("CWH", bgc)~"CWH",TRUE~NA)))
               
               
broads_tab<-group_by(broads, sppsplit, bgc)%>%summarise(nratings_ss=n())
broads_tab2<-group_by(broads, sppsplit)%>%mutate(nratings_spp=n())%>%select(sppsplit, nratings_spp)%>%distinct(.)
broads_tab<-left_join(broads_tab, broads_tab2)

#write.csv(broads_tab, "broadleaf_table.csv")
