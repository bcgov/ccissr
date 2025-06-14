# Save internal csv to package data
# You run this manually to update package data from the csvs

library(data.table)
library(usethis)
library(readxl)
library(ccissr)

N1 <- fread("site_series.csv", encoding = "Latin-1")
N1[RealmClass == "", RealmClass := NA]
N1 <- N1[,.(SS_NoSpace,SiteSeriesLongName,RealmClass)]
use_data(N1, overwrite = T)

bgcs <- fread("tables/versioned/WNA_BGCs_Info_v13_2.csv")
bgcs_wna <- unique(bgcs[DataSet == "BC", BGC])

bgc_notin <- bgcs_wna[!bgcs_wna %in% stocking_standards$ZoneSubzone]
eda <- E1[BGC %in% bgcs_wna, ]
bc_ss <- unique(eda$SS_NoSpace)

ss_notin <- bc_ss[!bc_ss %in% stocking_standards$SS_NoSpace]
fwrite(data.table(bgc_notin),"BGC_NoStocking.csv")
fwrite(data.table(ss_notin), "SiteSeries_NoStocking.csv")

edatopic <- fread("../../../Downloads/edatopic.csv")
suit <- fread("../../../Downloads/suitability.csv")
ss <- fread("tables/versioned/SpecialSites_v13_2.csv")

#subzones_colours_ref <- parse_qml("../../../Downloads/WNAv13_v5_Subzones.qml")
sz_wna <- fread("../Common_Files/WNAv13_SubzoneCols.csv")
sz_bc <- fread("../CCISS_ShinyApp/app/BC_SubzoneColours_v13_6.csv")
sz_wna <- sz_wna[!classification %in% sz_bc$classification,]
sz_all <- rbind(sz_bc, sz_wna)
fwrite(sz_all, "../CCISS_ShinyApp/app/WNA_SZ_Cols_v13_6.csv")

SS <- ss[,.(SS_NoSpace,SpecialCode)]
SS <- SS[SpecialCode != "",]
E1 <- SS[edatopic, on = "SS_NoSpace"]
setcolorder(E1,c("Source","BGC","SS_NoSpace","Edatopic","SpecialCode"))
phases <- E1[grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace),]
E1 <- E1[!(grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace)),]
vars <- E1[grep("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1 <- E1[!grepl("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1_Phase <- rbind(phases,vars)
E1_Phase[,MainUnit := gsub("[a-z]$","",SS_NoSpace)]
E1_Phase[,MainUnit := gsub("\\.[1-9]$","",MainUnit)]

suit[,V1 := NULL]
setnames(suit, c("bgc", "ss_nospace", "sppsplit", "feasible", "spp", "newfeas", 
                 "mod", "outrange"))
S1 <- copy(suit)

zone_cols <- fread("../../../Downloads/WNAv13_Zone_colours_2.csv")
dat2 <- zone_cols[,.(ZONE,RGB)] |> unique()
setnames(dat2,c("classification","colour"))
fwrite(dat2, "WNAv13_ZoneCols.csv")

E1 <- fread("tables/versioned/Edatopic_v13_5.csv")
S1 <- fread("tables/suitability.csv")
N1 <- fread("./data-raw/data_tables/SiteSeries_names_v12_15.csv", encoding = "Latin-1")
#N1[,SiteSeriesLongName := gsub(pattern = "[\x80-\xff]", "",SiteSeriesLongName, perl = T)]


SS <- fread("tables/versioned/Special_SS_v13_1.csv")

covMat <- read.csv("data-raw/Feas_CovMat.csv", header = TRUE, row.names = 1)

S1[,Confirmed := NULL]
S1 <- S1[!is.na(Feasible),]
setnames(S1, old = "SppVar",new = "Spp")
S1[Spp %in% c("Fdi","Fdc"),Spp := "Fd"]
S1[Spp %in% c("Pli","Plc"),Spp := "Pl"]
S1[Spp %in% c("Sw","Se","Sxw"),Spp := "Sx"]
S1[Spp %in% c("Ss", "Sxl","Sxs"),Spp := "Ss"]
S1[Spp %in% c("Pyi","Pyc"),Spp := "Py"]
S1[Spp %in% c("Acb","Act"),Spp := "Ac"]

SIBEC <- fread("~/SIBEC_Modelled/PredSI_Sept2021_2.csv") 
setnames(SIBEC,old = "SppVar", new = "Spp")
SIBEC[Spp %in% c("Fdi","Fdc"),Spp := "Fd"]
SIBEC[Spp %in% c("Pli","Plc"),Spp := "Pl"]
SIBEC[Spp %in% c("Sw","Se","Sxw"),Spp := "Sx"]
SIBEC[Spp %in% c("Ss", "Sxl","Sxs"),Spp := "Ss"]
SIBEC[Spp %in% c("Pyi","Pyc"),Spp := "Py"]
SIBEC[Spp %in% c("Acb","Act"),Spp := "Ac"]
SIBECnew <- fread("~/PortfolioKiri/InputsGit/SI_to_add.csv")
SIBEC <- rbind(SIBEC, SIBECnew)
###import SI data (currently from BART)
setcolorder(SIBEC,c("SS_NoSpace","Spp","SIPred"))
setnames(SIBEC,c("SS_NoSpace","TreeSpp","MeanPlotSiteIndex"))
SIBEC <- unique(SIBEC)

TreeCols <- fread("data-raw/PortfolioSppColours.csv", header = TRUE) ##in package data

SS <- SS[,.(SS_NoSpace,SpecialCode)]
SS <- SS[SpecialCode != "",]
E1 <- SS[E1, on = "SS_NoSpace"]
setcolorder(E1,c("Source","BGC","SS_NoSpace","Edatopic","SpecialCode"))
phases <- E1[grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace),]
E1 <- E1[!(grepl("BEC",Source) & grepl("[0-9]a$|[0-9]b$|[0-9]c$",SS_NoSpace)),]
vars <- E1[grep("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1 <- E1[!grepl("\\.1$|\\.2$|\\.3$",SS_NoSpace),]
E1_Phase <- rbind(phases,vars)
E1_Phase[,MainUnit := gsub("[a-z]$","",SS_NoSpace)]
E1_Phase[,MainUnit := gsub("\\.[1-9]$","",MainUnit)]

R1 <- fread("./data-raw/data_tables/RuleTable.csv")
F1 <- fread("./data-raw/data_tables/FeasibilityLabels.csv", key = "SuitDiff")
T1 <- fread("./data-raw/data_tables/Tree speciesand codes_2.0_25Aug2021.csv", key = "TreeCode")

V1 <- fread("./data-raw/data_tables/Variables_ClimateBC.csv", key = "Code")
zones_colours_ref <- fread("./data-raw/data_tables/WNAv11_Zone_Colours.csv", key = "classification")
subzones_colours_ref <- fread("./data-raw/data_tables/WNAv12_3_SubzoneCols.csv", key = "classification")

# StockingStds
stocking_standards_v12 <- fread("./data-raw/data_tables/StockingStds/StockStands_v13_1.csv", key = c("Region", "ZoneSubzone","SiteSeries", "Species"), colClasses = c("Standard" = "numeric"))
stocking_info_v12 <- fread("./data-raw/data_tables/StockingStds/StockingInfo_v13.csv", encoding = "Latin-1", key = "Standard", colClasses = c("Standard" = "numeric"))
stocking_height_v12 <- fread("./data-raw/data_tables/StockingStds/StockingHeight_v13.csv", key = c("Standard", "Species"), colClasses = c("Standard" = "numeric"))
crosswalk <- fread("./data-raw/data_tables/StockingStds/Crosswalk.csv", key = "Modeled")
stocking_standards_v12[,V1 := NULL]

# Massaging data
# Some Standards end with CC, discarding them
stocking_info_v12[ , Standard := as.numeric(Standard)]
setkey(stocking_info_v12, "Standard")
# Remove standards with no matching Standard in main table
stocking_info_v12 <- stocking_info_v12[!is.na(Standard) & Standard %in% stocking_standards_v12$Standard]
stocking_height_v12 <- stocking_height_v12[!is.na(Standard) & Standard %in% stocking_standards_v12$Standard]

# # Duplicated pairs
# dupPairs <- function(data) {
#   data[duplicated(data[, j = key(data), with=FALSE]) | duplicated(data[, j = key(data), with=FALSE], fromLast = TRUE), j = .SD, by=key(data)]
# }
# 
# # Checks standards for duplicates
# dupPairs(stocking_standards_v12)
# dupPairs(stocking_info_v12)
# dupPairs(stocking_height_v12)

# Remove duplicates for now, keeping the first of each combination
remDups <- function(d) {
  d[!duplicated(d[, data.table::key(d), with = FALSE])]
}
stocking_standards_v12 <- remDups(stocking_standards_v12)
stocking_info_v12 <- remDups(stocking_info_v12)
stocking_height_v12 <- remDups(stocking_height_v12)


# Stocking standards formatting
stocking_standards <- data.table::copy(stocking_standards_v12)
stocking_standards[, Footnotes := list(list({x <- unname(do.call(c, .SD)); x[!x %in% c(NA, "")]})), by=1:NROW(stocking_standards), .SDcols = FN1:FN5]
stocking_standards[, c("FN1","FN2","FN3","FN4","FN5","FN6") := NULL]
# add-in crosswalk rows to complete standards dataset
# Gettings standards that would be substitute according to crosswalk
a <- stocking_standards[ZoneSubzone %chin% crosswalk$Tables]
# Generate all possible BGC that would use a substitute
a <- a[crosswalk, on = c(ZoneSubzone = "Tables"), allow.cartesian = TRUE, nomatch = NULL]
# Checking if any of those already have a match in the standards table
nrow(stocking_standards[a, on = c(Region = "Region", ZoneSubzone = "Modeled", SS_NoSpace = "SS_NoSpace", Species = "Species"), nomatch = NULL])
# Does not seems like it, so it is safe to add all of them
a[, `:=`(ZoneSubzone = Modeled, Modeled = NULL)]
k <- data.table::key(stocking_standards)
stocking_standards <- rbindlist(list(stocking_standards, a))
setkeyv(stocking_standards, k)
# Recheck for dups
#dupPairs(stocking_standards)
#stocking_standards[,SiteSeries := gsub("[^0-9.-]", "", SiteSeries)]
stocking_standards[,SS_NoSpace := paste0(ZoneSubzone,"/",SiteSeries)]
stocking_standards[,SiteSeries := NULL]
region_cw <- data.table(RegionNew = unique(stocking_standards$Region), 
                        RegionOld = c("Kamloops","Vancouver","Nelson","Nelson","Pr Rupert","Pr Rupert","Vancouver","Cariboo","Pr George","Nelson","Nelson","Nelson","Kamloops","Pr George"))
setnames(stocking_standards,old = "Region",new = "RegionNew")
stocking_standards[region_cw, Region := i.RegionOld, on = "RegionNew"]
stocking_standards[,RegionNew := NULL]

#stocking_standards[,FN6 := NULL]
temp <- stocking_standards[grep("BWBS",ZoneSubzone),]
temp[,Region := "Pr Rupert"]
stocking_standards <- rbind(stocking_standards,temp)
stocking_standards[Species %in% c("Se","Sx","Sxw","Sw","Sxs"),Species := "Sx"]
stocking_standards[Species %in% c("Act","Acb","Aca"),Species := "Ac"]

cfrg_rules <- fread("data-raw/PreferredAcceptibleRules.csv")
cfrg_rules <- melt(cfrg_rules,id.vars = "Spp",variable.name = "Feasible",value.name = "PrefAcc")
cfrg_rules[,Feasible := as.character(gsub("E","",Feasible))]

# Stocking height formatting
stocking_height <- copy(stocking_height_v12[,.(Standard, Species, Height)])
# add-in species instead of others by merging with stock standards and removing dups as
# they will occurs in a , which comes later in the order
a <- stocking_height["Others", on = .(Species)][stocking_standards, on = c(Standard = "Standard"), nomatch = NULL]
a <- a[, Species := i.Species][, .(Standard, Species, Height)]
stocking_height <- rbindlist(list(stocking_height[Species != "Others"], a))
stocking_height <- stocking_height[!duplicated(stocking_height)]
# flag heights where there is a Species with a suit 1,2,3, or PrefAccept in P,A
stocking_height[stocking_standards, Flag := !is.na(i.Species), on = c(Standard = "Standard", Species = "Species")]


# Stocking info formatting
# Replace non ascii characters in SiteSeriesName and trim
# Remove extras spaces
stocking_info <- copy(stocking_info_v12)
stocking_info[, SiteSeriesName := {
  x <- chartr("\U2013\U2019\U0024\U00A0","-'  ", SiteSeriesName)
  x <- trimws(x)
  x <- gsub("\\s+;", ";", x)
  x <- gsub("\\s{2,}", " ", x)
  x
}]

# Notes
footnotes <- read_xlsx("./data-raw/data_tables/StockingStds/Revised Reference Guide Footnotes.xlsx", "Sheet1")
setDT(footnotes)
footnotes <- footnotes[Remove %chin% c(NA, "")]
footnotes[, `Restrictive Footnote` := !`Restrictive Footnote` %chin% c(NA, "")]
footnotes[, `Advisory Footnote` := !`Advisory Footnote` %chin% c(NA, "")]
footnotes[, `Geographic Restriction` := !`Geographic Restriction` %chin% c(NA, "")]
footnotes[, `Site Condition` := !`Site Condition` %chin% c(NA, "")]
footnotes[, `Pest/Disease` := !`Pest/Disease` %chin% c(NA, "")]
footnotes[, Administrative := !Administrative %chin% c(NA, "")]

# Silvical characteristics
silvics_tol <- setDT(read_xlsx("./data-raw/data_tables/StockingStds/EcoSilvic Tables_Klinka.xlsx", range = "Tolerance!A2:I35"))
names(silvics_tol)[4:9] <- paste(names(silvics_tol)[4:9], "tolerance")
silvics_regen <- setDT(read_xlsx("./data-raw/data_tables/StockingStds/EcoSilvic Tables_Klinka.xlsx", range = "Silvics!A2:H35"))
silvics_mature <- setDT(read_xlsx("./data-raw/data_tables/StockingStds/EcoSilvic Tables_Klinka.xlsx", range = "Silvics2!A2:I35"))
silvics_resist <- setDT(read_xlsx("./data-raw/data_tables/StockingStds/EcoSilvic Tables_Klinka.xlsx", range = "Resistance!A3:I36"))
names(silvics_resist) <- c("Life form", "Tree code", "Tree species", "Snow Resistance Class", "Wind Resistance Class",
                           "Fire Risk Class", "Insect Risk Class", "Fungi Risk Class", "Other Risk Class")

names(silvics_tol) <- tools::toTitleCase(names(silvics_tol))
names(silvics_regen) <- tools::toTitleCase(names(silvics_regen))
names(silvics_mature) <- tools::toTitleCase(names(silvics_mature))
names(silvics_resist) <- tools::toTitleCase(names(silvics_resist))
names(footnotes) <- tools::toTitleCase(names(footnotes))
names(stocking_standards) <- tools::toTitleCase(names(stocking_standards))
names(stocking_info) <- tools::toTitleCase(names(stocking_info))
names(stocking_height) <- tools::toTitleCase(names(stocking_height))

stocking_info <- stocking_info[!is.na(StockingTarget),]

# models informations
models_info <- fread("./data-raw/data_tables/CCISS_DataTable_Versions.csv")
models_info[, Date := as.character(Date, format = "%Y/%m/%d")]
subzones_colours_ref <- fread("../Common_Files/WNAv13_v6_SubzoneCols.csv")

use_data(stocking_standards,stocking_info,stocking_height, overwrite = TRUE)

use_data(E1, E1_Phase, S1, SS, N1, R1, F1, T1, V1,
         cfrg_rules, SIBEC, covMat,
         zones_colours_ref, subzones_colours_ref,
         stocking_standards, stocking_info, stocking_height, footnotes,
         silvics_tol, silvics_regen, silvics_mature, silvics_resist,
         models_info, TreeCols,
         overwrite = TRUE)


use_data(N1, overwrite = T)
use_data(E1, E1_Phase, S1, overwrite = TRUE)
# see version in ?usethis::use_data, if you all use R 3.5 and up. You should bump to version 3
# use_data(E1, S1, R1, F1, zones_colours_ref, subzones_colours_ref, overwrite = TRUE, version = 3)

# This will document your dataset in R/_data.R. See https://roxygen2.r-lib.org/articles/rd.html#datasets
# if you want to document them individually
writeLines(c(
"#' Data to be included in ccissr package",
"#'",
"#' @name bccciss-data",
"#' @docType data",
"#' @keywords data",
paste("#' @aliases", paste(gsub("\\.rda", "", dir("./data")), collapse = " ")),
"NULL"
), "./R/z_data.R")


##check edatopic data
#eda <- fread("~/CommonTables/Edatopic_v12_8.csv")
eda <- E1
posPosition <- as.data.table(expand.grid(SNR = c("A","B","C","D","E"),SMR = 0:7,BGC = unique(eda$BGC)))
posPosition[,Edatopic := paste0(SNR,SMR)]
t1 <- eda[posPosition,on = c("BGC","Edatopic")]
missed <- t1[is.na(SS_NoSpace),.(BGC,Edatopic)]