##======================================================================================
## CCISS Publication Scripts
## Step 2 - Species suitability projections
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019

# What this script does: calculate suitability outputs (but not stocking standards) equivalent to the CCISS tool, but simplified for fast running. 
# reason for script: the full CCISS script takes ~10-20 seconds per POI, meaning weeks of running time for a publishable raster grid. 
# Objective: produce results for a ~250,000-cell grid at a run time of maximum 8 hours.  
# general approach: this script achieves speed by: 
# 1. running only for selected edatopes
# 2. running on the whole grid vector, not looping on the POIs. 
# 3. outputs and inputs in many small tables, rather than a few big tables. 
# 4. running only for selected species

# Analysis Steps: 
#   /Edatope/RCP/Spp/proj.year/GCM: 
#     1. Project BGC unit 
#     2. find the site series associated with the edatope
#     3. obtain suitability of spp in site series

source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts




#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- fread(paste("./inputs/", grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
sort(table(BGC))

#===============================================================================
# Import BGC projections for each period
#===============================================================================

BGC.pred.ref <- as.character(read.csv(paste("outputs/BGC.pred", grid, "ref", model,"csv", sep="."), header = F)[,1])
zone.pred.ref <- rep(NA, length(BGC))
for(i in BGCcolors$classification){ zone.pred.ref[grep(i,BGC.pred.ref)] <- i }

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("outputs/BGC.pred", grid,hist.year, model,"csv", sep="."), header = F)[,1])
  assign(paste("BGC.pred", hist.year, sep="."), BGC.pred) #bgc projection
  print(hist.year)
}

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("outputs\\BGC.pred",grid, GCM, rcp, proj.year, model,"csv", sep="."), header = F)[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred) #bgc projection
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
      # print(GCM)
    }
    print(proj.year)
  }
  print(rcp)
}

#===============================================================================
# import will's lookup table for the site series associated with the selected edatope in each BGC unit
#===============================================================================

SiteLookup <- data.frame(BGC=unique(SiteSeries_Use$BGC))
edatopes<- c("B2", "C4", "D6")
for(edatope in edatopes){
  # SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[match(SiteLookup[,1], SiteSeries_Use$MergedBGC[which(SiteSeries_Use$Use==edatope)])])
  SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[which(SiteSeries_Use$Use==edatope)])
  names(SiteLookup)[which(edatopes==edatope)+1] <- edatope
}
str(SiteLookup)


#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

# Import suitability tables
S1 <- treesuit
S1 <- unique(S1)[,1:4]
dim(S1)
S1 <- unique(S1)
dim(S1)

## EDA: are there suitabilities for all projected units? 
NoSuit <- PredSum[-which(PredSum$BGC%in%S1$BGC),]
NoSuit[rev(order(NoSuit$count)),]

## EDA: Which site series are missing suitabilities? 
for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$SS_NoSpace),which(names(SiteLookup)==edatope)])
for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

## EDA: are there any units missing from the SiteSeries_Use table? 
BGClist <- unique(S1$BGC)
BGClist[-which(BGClist%in%SiteLookup$BGC)]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

# for(spp in c("Pl", "Fd", "Cw", "Sx")){
  for(spp in spps){
  for(edatope in edatopes){
    # get the suitability for the reference period predicted BGC.
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction
    
    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), as.character(S1$SS_NoSpace[which(S1$Spp==spp)]))]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$BGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    write.csv(Suit.ref, paste("outputs\\Suit.ref", grid, spp, edatope, "csv", sep="."), row.names = F)
    
    for(hist.year in hist.years){
      BGC.pred <- as.character(get(paste("BGC.pred", hist.year, sep=".")))
      ## identify cells with no suitability interpretations
      bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
      bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
      
      suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
      temp <- suit[match(BGC.pred, SiteLookup$BGC)]
      temp[is.na(temp)] <- 5 #set the NA values to suitability 5
      temp[temp==4] <- 5 #set 4 to suitability 5
      temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
      assign(paste("Suit", hist.year, sep="."), temp)
      write.csv(temp, paste("outputs\\Suit", grid, hist.year, spp, edatope, "csv", sep="."), row.names = F)
      # print(hist.year)
    }
    
    # get the suitability for future periods, for each projection/edatope/species combination
    for(GCM in GCMs){
      for(rcp in rcps){
        for(proj.year in proj.years){
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- as.character(get(paste("BGC.pred", GCM, rcp, proj.year, sep=".")))
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units
          
          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
          
          # get the suitability for the selected species associated with each site series
          suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$SS_NoSpace[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$BGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretation to NA
          assign(paste("Suit", GCM, rcp, proj.year, sep="."), temp)
          write.csv(temp, paste("outputs\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."), row.names = F)
          # print(proj.year)
        }
        # print(rcp)
      }
      # print(GCM)
    }
    # print(edatope)
  }
  print(paste(spp, " (", round(which(spps==spp)/length(spps)*100, 0), "%)", sep=""))
}
