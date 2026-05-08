
library(caret)


addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  return(dat)
}



vars <- read.csv("inputs\\53Var_WNABGCv11.csv", stringsAsFactors = F)[,1] # read in Will's set of biologically relevant variables
vars <- vars[-grep("BGC", vars)] #remove BGC
vars <- vars[-grep("EMT|EXT|DD18|DD_18", vars)] #remove BGC


# ===================================
#   ecoprovince climate data
# ===================================

ecoprov.hist <- read.csv("inputs/ecoprov.ref.100pts.mean.csv")
ecoprov.proj <- read.csv("inputs/ecoprov.proj.100pts.mean.csv")
ecoprov.6190 <- read.csv("inputs/ecoprov.6190.100pts.mean.csv")
ecoprov.0118 <- read.csv("inputs/ecoprov.0118.100pts.mean.csv")

ecoprov.hist <- addVars(ecoprov.hist)
ecoprov.proj <- addVars(ecoprov.proj)
ecoprov.6190 <- addVars(ecoprov.6190)
ecoprov.0118 <- addVars(ecoprov.0118)

ecoprov.proj <- ecoprov.proj[-grep("HadGEM", ecoprov.proj$Year)]


ecoprovs <- levels(ecoprov.proj$id1)
ecoprovs.order <- c(4,3,8,9,2,7,1,10,5)
ecoprov.names <- c("Boreal Plains", "Central Interior", "Coast and Mountains", "Georgia Depression", "Northern Boreal Mountains", "Southern Alaska", "Sub-Boreal Interior", "Southern Interior Mountains", "Southern Interior", "Taiga Plains")
variables <- names(ecoprov.proj)[-c(1:6)]
variable.names <- read.csv("inputs/Variables_ClimateBC.csv")

variable.types <- rep(NA, length(variables))
variable.types[grep("PPT|DD|PAS|NFFD|Eref|FFP|CMD|MAP|MSP|AHM|SHM|Rad|MAR", variables)] <- "ratio"
variable.types[grep("Tmax|Tmin|Tave|MAT|MWMT|MCMT|TD|EMT|EXT|bFFP|eFFP", variables)] <- "interval"
variable.types[grep("RH", variables)] <- "pct"

Ystr <- strsplit(as.character(ecoprov.proj[,1]), "_")
scenario <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)
scenario[,3] <- substr(scenario[,3],1,4)
scenario[grep("Ensemble", scenario[,1]),1] <- "Ensemble-mean"

rcps <- sort(unique(scenario[,2]))
proj.years <- unique(scenario[,3])

## calculate change in each climate variable

ecoprov.change <- ecoprov.proj[-c(1:6)]
ecoprov.change[,] <- NA

for(ecoprov in ecoprovs){
  s <- which(ecoprov.proj$id1==ecoprov)
  j <- which(ecoprovs==ecoprov)
  for(variable in variables){
    i <- which(variables==variable)
    ecoprov.change[s,i] <- if(variable.types[i]%in%c("interval", "pct")) ecoprov.proj[s,i+6] - ecoprov.hist[j,i+5] else ecoprov.proj[s,i+6]/ecoprov.hist[j,i+5]
  }
  print(ecoprov)
}

ecoprov.change.0118 <- ecoprov.0118[-c(1:6)]
ecoprov.change.0118[,] <- NA

for(ecoprov in ecoprovs){
  s <- which(ecoprov.6190$id1==ecoprov)
  for(variable in variables){
    i <- which(variables==variable)
    ecoprov.change.0118[s,i] <- if(variable.types[i]%in%c("interval", "pct")) ecoprov.0118[s,i+6] - ecoprov.6190[s,i+6] else ecoprov.0118[s,i+6]/ecoprov.6190[s,i+6]
  }
  print(ecoprov)
}

# ===================================
#   RECENT: determine which variables to toss and keep, based on both spatial and temporal correlation
# ===================================

cor.spatial <- cor(ecoprov.hist[,which(names(ecoprov.hist)%in%vars)], use="complete.obs")
cor.temporal <- cor(ecoprov.change.0118[,which(names(ecoprov.change.0118)%in%vars)], use="complete.obs")

hist(abs(cor.spatial))
hist(abs(cor.temporal))

# product of spatial and temporal
cor.prod <- cor.spatial*cor.temporal
hist(abs(cor.prod))

#ordered list 
cor.list <- cor.prod
cor.list[lower.tri(cor.list,diag=TRUE)]=NA  #Prepare to drop duplicates and diagonal
cor.list=as.data.frame(as.table(cor.list))  #Turn into a 3-column table
cor.list=na.omit(cor.list)  #Get rid of the junk we flagged above
cor.list=cor.list[order(-abs(cor.list$Freq)),]    #Sort by highest correlation (whether +ve or -ve)


t=0.80 #threshold for correlation product
candidates.0118 <- unique(as.vector(unlist(cor.list[which(cor.list$Freq>t),1:2])))
toss.0118 <- vars[findCorrelation(cor.prod, cutoff = t)]
keep.0118 <- candidates[-which(candidates.0118%in%toss.0118)]
vars.final.0118 <- vars[-which(vars%in%toss.0118)]

toss.0118
keep.0118
vars.final.0118


# ===================================
#   GCMs: determine which variables to toss and keep, based on both spatial and temporal correlation
# ===================================

cor.spatial <- cor(ecoprov.hist[,which(names(ecoprov.hist)%in%vars)], use="complete.obs")
cor.temporal <- cor(ecoprov.change[,which(names(ecoprov.change)%in%vars)], use="complete.obs")

hist(abs(cor.spatial))
hist(abs(cor.temporal))

# product of spatial and temporal
cor.prod <- cor.spatial*cor.temporal
hist(abs(cor.prod))

#ordered list 
cor.list <- cor.prod
cor.list[lower.tri(cor.list,diag=TRUE)]=NA  #Prepare to drop duplicates and diagonal
cor.list=as.data.frame(as.table(cor.list))  #Turn into a 3-column table
cor.list=na.omit(cor.list)  #Get rid of the junk we flagged above
cor.list=cor.list[order(-abs(cor.list$Freq)),]    #Sort by highest correlation (whether +ve or -ve)


t=0.80 #threshold for correlation product
candidates.proj <- unique(as.vector(unlist(cor.list[which(cor.list$Freq>t),1:2])))
toss.proj <- vars[findCorrelation(cor.prod, cutoff = t)]
keep.proj <- candidates[-which(candidates.proj%in%toss.proj)]
vars.final.proj <- vars[-which(vars%in%toss.proj)]

toss.proj
keep.proj
vars.final.proj

write.csv(vars.final.proj, "inputs\\VarSet_Spatiotemporal_GCMsOnly.csv", row.names = F)


# compare 
sort(table(c(candidates.0118, candidates.proj))) # lots of overlap between the highly correlated variables
sort(table(c(toss.0118, toss.proj))) # less overlap in variables excluded. 
sort(table(c(keep.0118, keep.proj))) # little/no overlap in retained variables


# ===================================
#   COMBINED: same analysis but for product of all three correlations
# ===================================

cor.spatial <- cor(ecoprov.hist[,which(names(ecoprov.hist)%in%vars)], use="complete.obs")
cor.temporal.0118 <- cor(ecoprov.change.0118[,which(names(ecoprov.change.0118)%in%vars)], use="complete.obs")
cor.temporal.proj <- cor(ecoprov.change[,which(names(ecoprov.change)%in%vars)], use="complete.obs")

hist(abs(cor.spatial))
hist(abs(cor.temporal.0118))
hist(abs(cor.temporal.proj))

# product of spatial and temporal
cor.prod <- cor.spatial*cor.temporal.0118*cor.temporal.proj
hist(abs(cor.prod))

#ordered list 
cor.list <- cor.prod
cor.list[lower.tri(cor.list,diag=TRUE)]=NA  #Prepare to drop duplicates and diagonal
cor.list=as.data.frame(as.table(cor.list))  #Turn into a 3-column table
cor.list=na.omit(cor.list)  #Get rid of the junk we flagged above
cor.list=cor.list[order(-abs(cor.list$Freq)),]    #Sort by highest correlation (whether +ve or -ve)


t=0.80 #threshold for correlation product
candidates <- unique(as.vector(unlist(cor.list[which(cor.list$Freq>t),1:2])))
toss <- vars[findCorrelation(cor.prod, cutoff = t)]
keep <- candidates[-which(candidates%in%toss)]
vars.final <- vars[-which(vars%in%toss)]

toss
keep
vars.final

write.csv(vars.final, "inputs\\VarSet_Spatiotemporal_GCMsAndObs.csv", row.names = F)
