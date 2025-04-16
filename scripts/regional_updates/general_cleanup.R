#read in current table version
suit<-read.csv("tables/versioned/suitability_v13_16.csv")
suit$X<-NULL
names(suit)

#clean up moderator/reviewer initials 
unique(suit$mod)
cts<-group_by(suit, mod)%>%summarise(cts=n())
suit<-mutate(suit, mod= case_when(mod=="BG"|mod=="BGx"~"BG",
                                  mod=="AB"|mod=="Abx"|mod=="AB?"| mod=="ABx"~"AB",
                                  mod=="MRyan"|mod=="MR"~"MR",
                                  mod=="KEI"|mod=="KI"~"KI",
                                  mod=="WHM"~"WHM",
                                  mod=="DVM"~"DVM",
                                  mod=="VC"|mod=="REVIEW"~"VC", #ratings added based on BEC guides for sites that did not have a suit rating
                                  mod=="DS"~"DS",
                                  mod=="KSD"~"KSD",
                                  mod=="SAS-HAK"|mod=="Coast"~"SAS-HAK",
                                  mod=="HK" ~"HK", 
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

