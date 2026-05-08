
library(car)
library(RColorBrewer)


source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

## Load the input data

ecoprov.hist <- read.csv("inputs/ecoprov.ref.100pts.mean.csv")
ecoprov.proj <- read.csv("inputs/ecoprov.proj.100pts.mean.csv")
ecoprov.6190 <- read.csv("inputs/ecoprov.6190.100pts.mean.csv")
ecoprov.0118 <- read.csv("inputs/ecoprov.0118.100pts.mean.csv")
ecoprov.1118 <- read.csv("inputs/ecoprov.1118.100pts.mean.csv")
ecoprov.1718 <- read.csv("inputs/ecoprov.1718.100pts.mean.csv")
ecoprov.ts <- read.csv("inputs/ecoprov.ts.100pts.mean.csv")
# str(ecoprov.climate)

ecoprovs <- levels(ecoprov.proj$id1)
ecoprov.names <- c("Boreal Plains", "Central Interior", "Coast and Mountains", "Georgia Depression", "Northern Boreal Mountains", "Southern Alaska", "Sub-Boreal Interior", "Southern Interior Mountains", "Southern Interior", "Taiga Plains")
variables <- names(ecoprov.proj)[-c(1:6)]
variables.select <- variables[c(grep("_wt|_sp|_sm|_at", variables), 225:247)]
variables.select <- variables.select[-grep("RH|Rad|MAR", variables.select)]

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
ecoprov.change.1118 <- ecoprov.change.0118
ecoprov.change.1718 <- ecoprov.change.0118

for(ecoprov in ecoprovs){
  s <- which(ecoprov.6190$id1==ecoprov)
  for(variable in variables){
    i <- which(variables==variable)
    ecoprov.change.0118[s,i] <- if(variable.types[i]%in%c("interval", "pct")) ecoprov.0118[s,i+6] - ecoprov.6190[s,i+6] else ecoprov.0118[s,i+6]/ecoprov.6190[s,i+6]
    ecoprov.change.1118[s,i] <- if(variable.types[i]%in%c("interval", "pct")) ecoprov.1118[s,i+6] - ecoprov.6190[s,i+6] else ecoprov.1118[s,i+6]/ecoprov.6190[s,i+6]
    ecoprov.change.1718[s,i] <- if(variable.types[i]%in%c("interval", "pct")) ecoprov.1718[s,i+6] - ecoprov.6190[s,i+6] else ecoprov.1718[s,i+6]/ecoprov.6190[s,i+6]
  }
  print(ecoprov)
}




############
## scatter plots of summer vs winter precipitation change
############ 

rcp="rcp45"
proj.year=proj.years[1]
models <- unique(scenario[which(scenario[,2]==rcp),1])
mods <- c("ENS", "ACC", "CAN", "CCSM", "CESM", "CNRM", "CSIR", "GFDL", "GISS", "HAD", "INM", "IPSL", "MIRE", "MIR5", "MPI", "MRI")[which(unique(scenario[,1])%in%models)]

ecoprovs.order <- c(3,9,2,7,1,5)

var1="Tave_wt"
var2="PPT_wt"

colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)][-1]
set.seed(2)
ColScheme <- c(brewer.pal(n=12, "Paired"),sample(colors,length(models)-12))
ColScheme[11] <- "blue"

proj <- ecoprov.change[which(scenario[,2]== rcp & scenario[,3]== proj.year ), ]
proj.model <- scenario[which(scenario[,2]== rcp & scenario[,3]== proj.year ), 1]
# str(proj)

png(filename=paste("results\\CCISS2019.ensembleScatter",var1, var2, proj.year, rcp,"png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=11, res=300)

mat <- matrix(c(1,9,3,4,5,1,9,6,7,8,9,9,2,2,2),3, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,0.15,1,1,1), heights=c(1,1,0.2))   #set up the multipanel plot

par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, paste("Change in", variable.names$Variable[which(variable.names$Code==var2)]), srt=90,cex=1.2)

plot(1, type="n", axes=F, xlab="", ylab="")
text(1,1, paste("Change in", variable.names$Variable[which(variable.names$Code==var1)]), cex=1.2)

par(mgp=c(2,0.25,0))

par(mar=c(0.1,0.1,0.1,0.1)) 

for(ecoprov in ecoprovs[ecoprovs.order]){
  
  xlim <- c(0,c(3.5, 5, 7)[which(proj.years==proj.year)])
  ylim <- if(var2=="PPT_wt") c(0.8, 1.2) else c(0.65,1.35)
  plot(0, 0,  xlim=xlim, ylim=ylim, col="white", xlab="", ylab="", tck=0, yaxt="n", xaxt="n")
  lines(c(-99,99), c(1,1), lty=2, col="gray")
  lines(c(0,0), c(-99,99), lty=2, col="gray")
  
  par(xpd=T)
  if(ecoprov%in%ecoprovs[ecoprovs.order][c(1,4)]) axis(2, at=seq(0,2,0.1), labels=paste(seq(0,2,0.1)*100-100, "%", sep=""), las=2, tck=0)
  if(ecoprov%in%ecoprovs[ecoprovs.order][c(4,5,6)]) axis(1, at=seq(0,8,1), labels=seq(0,8,1), las=1, tck=0)
  par(xpd=F)
  
  s <- which(ecoprov.proj$id1==ecoprov & scenario[,2]== rcp & scenario[,3]== proj.year )
  x <- ecoprov.change[s, which(names(ecoprov.change)==var1)]
  y <- ecoprov.change[s, which(names(ecoprov.change)==var2)]
  x.0118 <- ecoprov.change.0118[, which(names(ecoprov.change.0118)==var1)][which(ecoprovs==ecoprov)]
  y.0118 <- ecoprov.change.0118[, which(names(ecoprov.change.0118)==var2)][which(ecoprovs==ecoprov)]
  x.1118 <- ecoprov.change.1118[, which(names(ecoprov.change.1118)==var1)][which(ecoprovs==ecoprov)]
  y.1118 <- ecoprov.change.1118[, which(names(ecoprov.change.1118)==var2)][which(ecoprovs==ecoprov)]
  
  points(x.0118,y.0118, pch=16, col="gray", cex=5)
  text(x.0118,y.0118, "2001-2018", cex=1.5, font=2, pos=4, col="gray", offset=1.1)
  
  for(model in models){
    i=which(models==model)
    points(x[i],y[i], pch=21, bg=ColScheme[i], cex=if(model==models[1]) 5.5 else 4)
    text(x[i],y[i], mods[i], cex=if(model==models[1]) 1 else 0.7, font=2)
    
  }
  
  mtext(paste("(", letters[which(ecoprovs[ecoprovs.order]==ecoprov)], ") ", ecoprov.names[which(ecoprovs==ecoprov)], sep=""), side=3, line=-1.5, adj=0.025, cex=0.8, font=2)
  box()
  
  print(ecoprov)
}

dev.off()
