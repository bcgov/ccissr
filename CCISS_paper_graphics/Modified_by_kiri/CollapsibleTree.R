###Run first part of predict script to get Y2.sub
library(collapsibleTree)
library(reshape2)

Y2.year <- Y2.sub[-c(4,8:10)]

###cast to hierarchical structure
Y2.year <- dcast(Y2.year, GCM + Scenario + BGC + SiteNo ~ FuturePeriod)
Y2.year <- Y2.year[order(Y2.year$SiteNo, Y2.year$`2025`),]

fcData <- Y2.year[Y2.year$SiteNo == 8,]
fcData$`2025` <- as.factor(fcData$`2025`)
fcData$`2055` <- as.factor(fcData$`2055`)
fcData$`2085` <- as.factor(fcData$`2085`)

###calculate percent going to each model and paste with subzone name
##2025###
percent <- data.frame(table(fcData$`2025`))
fcData <- merge(fcData, percent, by.x = "2025", by.y = "Var1", all = TRUE)
fcData$Perc25 <- fcData$Freq/length(fcData$`2025`)*100
fcData$Perc25 <- round(fcData$Perc25, digits = 0)
fcData$`2025` <- paste(fcData$Perc25, fcData$`2025`,  sep = "% ")

##2055####
fcData$Num <- rep(1, length(fcData$BGC))
fcData$Perc55 <- ave(fcData$Num,fcData$`2025`,fcData$`2055`, FUN = sum)
fcData$Count25 <- ave(fcData$Num,fcData$`2025`, FUN = sum)
fcData$Percent55 <- fcData$Perc55/fcData$Count25 * 100
fcData$Percent55 <- round(fcData$Percent55, digits = 0)
fcData$`2055` <- paste(fcData$Percent55, fcData$`2055`,  sep = "% ")

#2085##
fcData$Sum <- ave(fcData$Num,fcData$`2025`,fcData$`2055`, fcData$`2085`, FUN = sum)
fcData$Count <- ave(fcData$Num,fcData$`2025`, fcData$`2055`, FUN = sum)
fcData$Percent85 <- fcData$Sum/fcData$Count * 100
fcData$Percent85 <- round(fcData$Percent85, digits = 0)
fcData$`2085` <- paste(fcData$Percent85, fcData$`2085`,  sep = "% ")

##graph###
collapsibleTreeSummary(
  fcData,
  c("2025","2055","2085"),
  root = unique(fcData$BGC),
  collapsed = FALSE
)
