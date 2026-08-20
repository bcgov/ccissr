library(data.table)
library(usethis)


##suitability
S1 <- fread("https://raw.githubusercontent.com/bcgov/ccissr/refs/heads/feas_tables/tables/suitability.csv")
S1[,V1 := NULL]
setnames(S1, old = c("suitability","newsuit"), new = c("feasible","newfeas"))

E1 <- fread("https://raw.githubusercontent.com/bcgov/ccissr/refs/heads/feas_tables/tables/edatopic.csv", header = TRUE)
SS <- fread("https://raw.githubusercontent.com/bcgov/ccissr/refs/heads/feas_tables/tables/site_series.csv", header = TRUE, encoding = "Latin-1")
SS[SpecialCode == "", SpecialCode := NA]
SS <- SS[,.(SS_NoSpace,SpecialCode,SiteSeriesLongName)]
SS <- unique(SS)
E1[SS, SpecialCode := i.SpecialCode, on = "SS_NoSpace"]

WNA_BGCs <- fread("https://raw.githubusercontent.com/bcgov/ccissr/refs/heads/feas_tables/tables/WNA_BGCs_Info.csv", header = TRUE, encoding = "Latin-1")
WNA_BGCs <- WNA_BGCs[,.(BGC, Zone, SubzoneGroup, SubzoneColour = colour, ZoneColour, notinTHLB,ZoneName, SubzoneName, VariantName,DataSet, State)]
use_data(WNA_BGCs, overwrite = TRUE)

################
bc_units <- WNA_BGCs[DataSet == "BC", BGC]
ss_needed <- unique(E1[BGC %in% bc_units, SS_NoSpace])
ss_missing <- ss_needed[!ss_needed %in% SS$SS_NoSpace]

setcolorder(E1,c("Source","BGC","SS_NoSpace","Edatopic","SpecialCode"))
phases <- E1[grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace),]
E1 <- E1[!(grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace)),]
vars <- E1[grep("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1 <- E1[!grepl("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1_Phase <- rbind(phases,vars)
E1_Phase[,MainUnit := gsub("[a-z]$","",SS_NoSpace)]
E1_Phase[,MainUnit := gsub("\\.[1-9]$","",MainUnit)]
use_data(S1, E1, E1_Phase, SS, overwrite = TRUE)
