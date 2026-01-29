########### 
library(scales)
library(MASS)   #contains the lda; also eqscplot(x,y) function, which plots axes of equal scale
library(stats)
library(raster)
library(FNN)
library(RColorBrewer)
library(colorRamps)
library(fields) #for rdist distance calculation
library(rgdal)
# library(igraph)

type <- "SimpleMahal" #indicates that this method calculates Mahalanobis distance directly from the correlation matrix of the reference variability, without doing the PCA on the spatiotemporal variation.it follows from the PCA2only sensitivity analysis on the StationFocals type.  
VarSet <- "Seasonal"
VarCode <- "S"
Grid.fine <- "NAnaec4"
park.analogs <- "SingleParkReverse"
Grid.medium <- "parks.F32"
Grid.analog2 <- ""
Grid.analog3 <- ""
model.medium <- if(Grid.medium%in%c("parks.F16", "parks.F32", "parks.F64")) "15GCM-Ensemble" else "GlobalMean"
model.fine <- if(Grid.fine%in%c("parks.F16", "parks.F32", "parks.F64")) "15GCM-Ensemble" else "GlobalMean"
model.analog2 <- if(Grid.analog2%in%c("parks.F16", "parks.F32", "parks.F64")) "15GCM-Ensemble" else "GlobalMean"
model.analog3 <- if(Grid.analog3%in%c("parks.F16", "parks.F32", "parks.F64")) "15GCM-Ensemble" else "GlobalMean"
scenario <- "6variable"

setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")
complex.lookup <- read.csv(paste("ComplexLookup.csv", sep=""), stringsAsFactors = F, strip.white = T)
nicknames <- unique(complex.lookup$nickname)

nickname=nicknames[13]
for(nickname in nicknames){
  complex.parks <- complex.lookup$PROT_NAME[which(complex.lookup$nickname==nickname)]
  
proj.year=2085
RCP="RCP45"
analogPool="SinglePark"

      ##################
      #### Read in the input data 
      ##################
      ##these files were created in a separate script "BCparks_InputData.R"
      setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")
      
land.medium <- read.csv(paste("land.",Grid.medium,".csv", sep=""))[,1]  
nonCNA.medium <- read.csv(paste("nonCNA.",Grid.medium,".csv", sep=""))[,1]  
subsample.medium <- read.csv(paste("subsample.",Grid.medium,".csv", sep=""))[,1]  
X.grid.ref.medium <- read.csv(paste("X.",Grid.medium,".ref.csv",sep=""))
X.grid.proj.medium <- read.csv(paste("X.",Grid.medium,".proj_",model.medium,"_",RCP,"_",proj.year,".csv",sep=""))

land.fine <- read.csv(paste("land.",Grid.fine,".csv", sep=""))[,1]  
nonCNA.fine <- read.csv(paste("nonCNA.",Grid.fine,".csv", sep=""))[,1] 
subsample.fine <- read.csv(paste("subsample.",Grid.fine,".csv", sep=""))[,1] 
X.grid.ref <- read.csv(paste("X.",Grid.fine,".ref.csv",sep=""))
X.grid.proj <- read.csv(paste("X.",Grid.fine,".proj_",model.fine,"_",RCP,"_",proj.year,".csv",sep=""))

X.cru.ref <- read.csv("X.stn.ref.csv")  #note that these are the cru stations, not the cru surrogates used in previous versions of the process. these are the focal locations
X.cru.proj <- read.csv(paste("X.stn.proj_GlobalMean_",RCP,"_",proj.year,".csv",sep="")) ##note that these are the cru stations. these are the focal locations

A.grid.4nn <- read.csv(paste("grid4nn_",Grid.fine,".csv",sep=""))
A.grid.dist <- read.csv(paste("griddist_",Grid.fine,".csv",sep=""))

X.cru_detrended <- read.csv("X.stn_detrended.csv")
A.cru_detrended <- read.csv("A.stn_detrended.csv")[,1]
Year.cru_detrended <- read.csv("Year.stn_detrended.csv")[,1]


########################
### Modifications for Sensitivity Analysis
########################
# reduce analysis variables 
X.grid.ref.medium <- X.grid.ref.medium[,c(1,3,5,7,9,11)]
X.grid.proj.medium <- X.grid.proj.medium[,c(1,3,5,7,9,11)]
X.grid.ref <- X.grid.ref[,c(1,3,5,7,9,11)]
X.grid.proj <- X.grid.proj[,c(1,3,5,7,9,11)]
X.cru.ref <- X.cru.ref[,c(1,3,5,7,9,11)]
X.cru.proj <-  X.cru.proj[,c(1,3,5,7,9,11)]
X.cru_detrended <- X.cru_detrended[,c(1,3,5,7,9,11)]

#restrict analog pools to land only. 
X.grid.ref.medium <- X.grid.ref.medium[subsample.medium,]
X.grid.proj.medium <- X.grid.proj.medium[subsample.medium,]
dim(X.grid.ref.medium)

# If specified, subset the analog pool to the specified park
  parks <- read.csv(paste("BCparksSpatial.",Grid.medium,".csv", sep=""))[-nonCNA.medium,][subsample.medium,]
  dim(parks)
  X.grid.ref.medium <- X.grid.ref.medium[which(parks$PROT_NAME%in%complex.parks),]
  X.grid.proj.medium <- X.grid.proj.medium[which(parks$PROT_NAME%in%complex.parks),]

  
########################
### Pre-Processing
########################

## set the list of reference locations to be analysed
cru.id <- seq(1,length(X.cru.ref[,1]))
focals <- sort(unique(as.vector(as.matrix(A.grid.4nn))))  #gets all the cru stations that are used by grid cells in the study area

#create dummy records to ensure that the loop can operate on the full set of cru.ids in all iterations. 
A.grid.4nn <- rbind(A.grid.4nn,matrix(rep(cru.id,times=4),length(cru.id),4))    #add a complete set of CRU.ids so that there are none missing for the local stPCA loop. 
dummies <- seq(length(land.fine)+1,length(A.grid.4nn[,1]),1)
X.grid.ref[dummies,] <- 1     #add dummy records to the reference grid
X.grid.proj[dummies,] <- 1    #add dummy records to the projected grid

#initiate the data frame to store the projected novelty and disappearance distances (note that these include the dummy records created above)
NN.dist.proj <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
NN.id.proj <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
# Analog.dist.proj <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
NN.dist.ref <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
NN.id.ref <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
# Analog.dist.ref <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
NN.dist.his <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
NN.id.his <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))
# Analog.dist.his <- data.frame(matrix(rep(NA,length(A.grid.4nn[,1]),times=4), nrow=length(A.grid.4nn[,1]), ncol=4))

# initiate data frames to store the scree vectors for each station PCA
PC2.sdev <- data.frame(matrix(NA,length(cru.id),length(names(X.cru.ref))))
ccsignal <- data.frame(matrix(NA,length(cru.id),length(names(X.cru.ref))))
spatialsignal <- data.frame(matrix(NA,length(cru.id),length(names(X.cru.ref))))

#create the data frame of reference variability at each cru surrogate
X.MeanGroupStDev <- aggregate(X.cru_detrended,list(A.cru_detrended),FUN=sd, na.rm=T)[,2:(length(X.cru_detrended)+1)]


##################
##### LOCAL STPCA: Generate the grids of Novelty distance by looping local-specific spatiotemporal standardization. 
##################

trunc.SD2s <- 0.1
trunc.rule2 <- paste(trunc.SD2s,"SD2",sep="")

# Start the clock
  ptm <- proc.time()

for(i in focals){       #this loop runs at 1.5 it/sec using the 4km grid
  
  members <- unique(c(which(A.grid.4nn[1:length(land.fine),1]==i),which(A.grid.4nn[1:length(land.fine),2]==i),which(A.grid.4nn[1:length(land.fine),3]==i),which(A.grid.4nn[1:length(land.fine),4]==i))) #define the fine-grid cells that are associated with the focal cru surrogate, based on that cru surrogate being one of the four nearest neighbours of the grid cell. note that this vector doesn't and shouldn't include the dummy records, because it is used to the define the study area subset. 
  
  #### Step 2: standardize the data using detrended reference period variability for cru surrogate i (i.e. spherize the temporal climate envelope of the focal CRU surrogate)
  X.cru_detrended.std <- sweep(X.cru_detrended[which(A.cru_detrended==i),],MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)
  X.grid.ref.medium.std <- sweep(X.grid.ref.medium,MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #all medium-res grid cells in the local subset, since these constitute the available analogs and the observations for stPCA
  X.grid.proj.medium.std <- sweep(X.grid.proj.medium,MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #all medium-res grid cells in the local subset, since these constitute the available analogs and the observations for stPCA
  X.grid.ref.std <- sweep(X.grid.ref[c(members,dummies[i]),],MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #high-res grid cells within the focal CRU surrogate
  X.grid.proj.std <- sweep(X.grid.proj[c(members,dummies[i]),],MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #high-res grid cells within the focal CRU surrogate
  X.cru.ref.std <- sweep(X.cru.ref[i,],MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #high-res grid cells within the focal CRU surrogate
  X.cru.proj.std <- sweep(X.cru.proj[i,],MARGIN=2,unlist(X.MeanGroupStDev[i,]),`/`)   #high-res grid cells within the focal CRU surrogate
  
  #### Step 4b: second rotation; a pca on the truncated reference variability. 
  PCA2 <- prcomp(X.cru_detrended.std[!is.na(apply(X.cru_detrended.std,1,mean)),])   #the apply() term is there simply to select all years with complete observations. 
  PC2.sdev[i,1:length(unlist(summary(PCA2)[1]))] <- unlist(summary(PCA2)[1])
  PC2s <- max(which(PC2.sdev[i,]>trunc.SD2s))    #PC truncation rule of eigenvector stdev > a specified threshold
  
  #### Step 4c: project the data onto the PCs
  Z.cru_detrended<- as.data.frame(predict(PCA2,X.cru_detrended.std))
  Z.grid.ref.medium<- as.data.frame(predict(PCA2,X.grid.ref.medium.std))
  Z.grid.proj.medium<- as.data.frame(predict(PCA2,X.grid.proj.medium.std))
  Z.grid.ref<- as.data.frame(predict(PCA2,X.grid.ref.std))
  Z.grid.proj<- as.data.frame(predict(PCA2,X.grid.proj.std))
  Z.cru.ref<- as.data.frame(predict(PCA2,X.cru.ref.std))
  Z.cru.proj<- as.data.frame(predict(PCA2,X.cru.proj.std))
  
  ##Step 5: Standardize the Z data based on the detrended variability of CRU surrogate i (i.e. re-spherize the temporal climate envelope)
  Z.MeanGroupStDev <- apply(Z.cru_detrended,2,sd, na.rm=T)     #Create the vector of standard deviations of reference variability in each variable
  Z.grid.ref.medium.sphere <- sweep(Z.grid.ref.medium,MARGIN=2,Z.MeanGroupStDev,`/`)     
  Z.grid.proj.medium.sphere <- sweep(Z.grid.proj.medium,MARGIN=2,Z.MeanGroupStDev,`/`)     
  Z.grid.ref.sphere <- sweep(Z.grid.ref,MARGIN=2,Z.MeanGroupStDev,`/`)     
  Z.grid.proj.sphere <- sweep(Z.grid.proj,MARGIN=2,Z.MeanGroupStDev,`/`)     
  Z.cru.ref.sphere <- sweep(Z.cru.ref,MARGIN=2,Z.MeanGroupStDev,`/`)     
  Z.cru.proj.sphere <- sweep(Z.cru.proj,MARGIN=2,Z.MeanGroupStDev,`/`)     
  
  #calculate the climate change signal for the CRU station
  ccsignal[i,1:length(unlist(summary(PCA2)[1]))] <- Z.cru.proj.sphere-Z.cru.ref.sphere
  spatialsignal[i,1:length(unlist(summary(PCA2)[1]))] <- apply(Z.grid.ref.medium.sphere,2,sd, na.rm=T)
  
  #loop the Novelty measurements for each of the four sets of fine grid cells associated with the focal cru surrogate
  for(j in 1:4){ 
    NN <- get.knnx(data=Z.grid.ref.medium.sphere[,1:PC2s],query=Z.grid.proj.sphere[which(A.grid.4nn[c(members,dummies[i]),j]==i),1:PC2s],k=1,algorithm="brute")
    NN.dist.proj[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[2]])     #Novelty distance. M distance from projected condition of hi-res member cells to their nearest neighbour reference medium-res grid cell 
    NN.id.proj[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[1]])     #Analog ID. SEQUENCE NUMBER of hi-res member cells to their nearest neighbour reference medium-res grid cell.
    # Analog.dist.proj[which(A.grid.4nn[,j]==i),j] <- rdist.vec(Z.grid.ref.sphere[which(A.grid.4nn[c(members,dummies[i]),j]==i),1:PC2s],Z.grid.ref.medium.sphere[as.vector(NN[[1]]),1:PC2s])
    }

  #loop the disappearance measurements for each of the four sets of fine grid cells associated with the focal cru surrogate
  for(j in 1:4){ 
    NN <- get.knnx(data=Z.grid.proj.medium.sphere[,1:PC2s],query=Z.grid.ref.sphere[which(A.grid.4nn[c(members,dummies[i]),j]==i),1:PC2s],k=1,algorithm="brute")
    NN.dist.ref[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[2]])     #Novelty distance. M distance from projected condition of hi-res member cells to their nearest neighbour reference medium-res grid cell 
    NN.id.ref[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[1]])     #Analog ID. SEQUENCE NUMBER of hi-res member cells to their nearest neighbour reference medium-res grid cell.
    # Analog.dist.ref[which(A.grid.4nn[,j]==i),j] <- rdist.vec(Z.grid.proj.sphere[which(A.grid.4nn[c(members,dummies[i]),j]==i),1:PC2s],Z.grid.proj.medium.sphere[as.vector(NN[[1]]),1:PC2s])
  }

  #loop the reference similarity measurements for each of the four sets of fine grid cells associated with the focal cru surrogate
  for(j in 1:4){ 
    NN <- get.knnx(data=Z.grid.ref.medium.sphere[,1:PC2s],query=Z.grid.ref.sphere[which(A.grid.4nn[c(members,dummies[i]),j]==i),1:PC2s],k=1,algorithm="brute")
    NN.dist.his[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[2]])     #Novelty distance. M distance from projected condition of hi-res member cells to their nearest neighbour reference medium-res grid cell 
    NN.id.his[which(A.grid.4nn[,j]==i),j] <- as.vector(NN[[1]])     #Analog ID. SEQUENCE NUMBER of hi-res member cells to their nearest neighbour reference medium-res grid cell.
  }
  print(which(focals==i))
  
}

####DELETE THE DUMMY RECORDS FROM All DATA FRAMES
A.grid.4nn <- A.grid.4nn[-dummies,]    
X.grid.ref <- X.grid.ref[-dummies,]
X.grid.proj <- X.grid.proj[-dummies,]    
NN.dist.proj <- NN.dist.proj[-dummies,]    
NN.id.proj <- NN.id.proj[-dummies,]    
# Analog.dist.proj <- Analog.dist.proj[-dummies,]   
NN.dist.ref <- NN.dist.ref[-dummies,]    
NN.id.ref <- NN.id.ref[-dummies,]    
# Analog.dist.ref <- Analog.dist.ref[-dummies,]   
NN.dist.his <- NN.dist.his[-dummies,]    
NN.id.his <- NN.id.his[-dummies,]    
# Analog.dist.his <- Analog.dist.his[-dummies,]   

setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\OutputData\\")
dir.create(nickname)
setwd(paste("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\OutputData\\", nickname, sep=""))
write.csv(NN.dist.proj,paste("NNdistproj",Grid.fine,park.analogs , Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(NN.id.proj,paste("NNidproj",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
# write.csv(Analog.dist.proj,paste("Analogdistproj",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(NN.dist.ref,paste("NNdistref",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(NN.id.ref,paste("NNidref",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
# write.csv(Analog.dist.ref,paste("Analogdistref",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(NN.dist.his,paste("NNdisthis",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(NN.id.his,paste("NNidhis",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
# write.csv(Analog.dist.his,paste("Analogdisthis",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(PC2.sdev,paste("PC2sdev",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(ccsignal,paste("ccsignal",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
write.csv(spatialsignal,paste("spatialsignal",Grid.fine,park.analogs ,Grid.medium,Grid.analog2,Grid.analog3,model.fine,RCP,proj.year,trunc.rule2,type,nickname,scenario,".csv", sep="_"), row.names=FALSE)
setwd("D:\\Backup_Jun12\\Masters\\Research\\Data\\BCparks\\InputData")

print(nickname)
    }

