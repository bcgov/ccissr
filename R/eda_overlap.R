## CCISS 2020 Step 2: Edatopic Overlap
## Kiri Daust, 2020

# Copyright 2021 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# E1 <- fread("~/CommonTables/Edatopic_v12_3.csv")
# E1[Special == "", Special := NA]
# E1[Codes == "", Codes := NA]
# E1 <- E1[Edatopic != "",]
# save(E1,file = "./data/E1.rda")
# BGC = test; Edatope = E1

#' EDA Topic Overlap
#' @param BGC BGC
#' @param E1 A data.table?
#' @param E1_Phase A phase?
#' @details What the function does
#' @return What the function returns
#' @import data.table
#' @importFrom dplyr left_join distinct
#' @importFrom stats complete.cases na.omit
#' @export
edatopicOverlap <- function(BGC,E1,E1_Phase,onlyRegular = FALSE){
  
  # Declare binding for checks
  if (FALSE) {
    SS_NoSpace <- 
      Edatopic <- 
      SpecialCode <- 
      BGC.pred <- 
      SS.pred <- 
      SiteRef <- 
      FuturePeriod <- 
      BGC.prop <- 
      SS.Curr <- 
      i.NumEdas <- 
      SSProb <- 
      SS.prob <- 
      SSProbRev <- 
      allOverlap <- 
      MainUnit <- 
      Phase <- 
      i.BGC <- 
      PhasePred <- 
      MainUnit.y <- 
      PhaseSum <- 
      phaseOverlap <- 
      allOverlap.y <- 
      Flag <- 
      SSratio <- 
      SSprob <- 
      NULL
  }
  
  SS <- E1[,list(BGC,SS_NoSpace,Edatopic)]
  edaPhase <- E1_Phase
  SS <- unique(SS)
  BGC <- unique(BGC)
  SSsp <- E1[!is.na(SpecialCode),list(BGC,SS_NoSpace,SpecialCode)]
  SSsp <- unique(SSsp)
  SSsp_phase <- edaPhase[!is.na(SpecialCode),list(BGC,SS_NoSpace,SpecialCode)]
  edaPhase <- edaPhase[is.na(SpecialCode),!"SpecialCode"]
  SSsp <- rbind(SSsp,SSsp_phase)
  
  ##Special site series edatopes
  CurrBGC <- SSsp[BGC, on = "BGC", allow.cartesian = T]
  setkey(BGC, BGC.pred)
  setkey(SSsp, BGC)
  FutBGC <- SSsp[BGC, allow.cartesian = T]
  setnames(FutBGC, old = c("BGC","SS_NoSpace","i.BGC"), 
           new = c("BGC.pred","SS.pred","BGC"))
  FutBGC <- FutBGC[!is.na(SS.pred),]
  setkey(FutBGC, SiteRef, FuturePeriod, BGC,BGC.pred, SpecialCode)
  setkey(CurrBGC,SiteRef,FuturePeriod, BGC,BGC.pred, SpecialCode)
  # "CWHvh2" BGC gives out Join results in 258 rows; more than 135 = nrow(x)+nrow(i), I'm
  # setting it to allow.cartesian, might need to investigate.
  new <- CurrBGC[FutBGC, allow.cartesian=TRUE]
  new <- new[!is.na(SS_NoSpace),] ##this removes special site series that don't align
  SSsp.out <- new[,list(allOverlap = 1/.N,SS.pred,BGC.prop), keyby = list(SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace)]
  SSsp.out <- unique(SSsp.out)
  
  ##regular site series edatopes
  SS <- E1[,list(BGC,SS_NoSpace,Edatopic)]
  temp <- rbind(SS,E1_Phase[is.na(SpecialCode),list(BGC,SS_NoSpace,Edatopic)])
  CurrBGC <- temp[BGC, on = "BGC", allow.cartesian = T]
  CurrBGC <- CurrBGC[!duplicated(CurrBGC),]
  setkey(BGC, BGC.pred)
  setkey(SS, BGC)
  FutBGC <- SS[BGC, allow.cartesian = T]
  FutBGC <- FutBGC[!duplicated(FutBGC),] 
  setnames(FutBGC, old = c("BGC","SS_NoSpace","i.BGC"), 
           new = c("BGC.pred","SS.pred","BGC"))
  FutBGC <- na.omit(FutBGC)
  
  setkey(FutBGC, SiteRef, FuturePeriod, BGC,BGC.pred, Edatopic)
  FutBGC[,BGC.prop := NULL]
 
  setkey(CurrBGC,SiteRef,FuturePeriod, BGC,BGC.pred, Edatopic)
  new <- FutBGC[CurrBGC, allow.cartesian = T]
  new <- new[!is.na(SS.pred),]
  setkey(new, SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)
  ##new <- new[complete.cases(new),]
  
  numEda <- E1[,list(NumEdas = .N), by = list(BGC,SS_NoSpace)]
  
  ###forwards overlap
  SS.out <- new[,list(SS.prob = .N,BGC.prop = BGC.prop[1]), 
                keyby = list(SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)]
  SS.out[numEda,SS.Curr := i.NumEdas, on = c(SS_NoSpace = "SS_NoSpace"), allow.cartesian = T]
  SS.out[,SSProb := SS.prob/SS.Curr]
  
  ###reverse overlap
  SS.out.rev <- new[,list(SS.prob = .N,BGC.prop = BGC.prop[1]), 
                    keyby = list(SiteRef,FuturePeriod,BGC,BGC.pred,SS.pred,SS_NoSpace)]
  SS.out.rev[numEda,SS.Curr := i.NumEdas, on = c(SS.pred = "SS_NoSpace")]
  SS.out.rev[,SSProbRev := SS.prob/SS.Curr]

  ##combine them
  combAll <- merge(SS.out,SS.out.rev,by = c("SiteRef","FuturePeriod","BGC","BGC.pred","SS_NoSpace","SS.pred"))
  combAll[,allOverlap := SSProb*SSProbRev]
  setnames(combAll, old = "BGC.prop.x",new = "BGC.prop")
  combAll <- combAll[,list(SiteRef, FuturePeriod, BGC, BGC.pred, SS_NoSpace, 
                        allOverlap, SS.pred, BGC.prop)]
  combAll <- na.omit(combAll)
  #combAllSave <- combAll
  combAll <- unique(combAll)
  combAll <- merge(combAll,SSsp.out,
                   by = c("SiteRef","FuturePeriod","BGC","BGC.pred","SS_NoSpace","SS.pred"), all = T)
  ##which ones are not current special site series
  combAll[E1[,.(SS_NoSpace,SpecialCode)], SppCode := i.SpecialCode, on = "SS_NoSpace"]
  combAll[E1[,.(SS_NoSpace,SpecialCode)], SppCodePred := i.SpecialCode, on = c(SS.pred = "SS_NoSpace")]
  combAll[SppCode != SppCodePred, Flag := F]
  combAll <- combAll[is.na(Flag),]
  combAll[,Flag := if(all(is.na(allOverlap.y))) T else F, by = list(SiteRef,FuturePeriod,BGC,SS_NoSpace,BGC.pred)]
  combAll[is.na(SppCode) & !is.na(SppCodePred), rmFlag := T]
  combAll[!is.na(SppCode) & is.na(SppCodePred) & Flag == FALSE, rmFlag := T]
  
  combAll <- combAll[is.na(rmFlag),]
  combAll[!is.na(allOverlap.y), `:=`(allOverlap.x = allOverlap.y, BGC.prop.x = BGC.prop.y)]
  combAll[,c("allOverlap.y","BGC.prop.y","Flag","rmFlag","SppCode","SppCodePred") := NULL]
  setnames(combAll,old = c("allOverlap.x","BGC.prop.x"), new = c("allOverlap","BGC.prop"))
  #combAll[!is.na(allOverlap.y), `:=`(allOverlap.x = allOverlap.y, BGC.prop.x = BGC.prop.y)]
  #temp <- combAll[!is.na(allOverlap.y),]
  #temp[,c("allOverlap.x","BGC.prop.x") := NULL]
  
  # setnames(temp,old = c("allOverlap.y","BGC.prop.y"), new = c("allOverlap","BGC.prop"))
  # combAll[,Flag := if(all(is.na(allOverlap.y))) T else F, by = list(SiteRef,FuturePeriod,BGC,SS_NoSpace,BGC.pred)]
  # combAll <- combAll[(Flag),!"Flag"]
  # combAll[,c("allOverlap.y","BGC.prop.y") := NULL]
  # setnames(combAll,old = c("allOverlap.x","BGC.prop.x"), new = c("allOverlap","BGC.prop"))
  # combAll <- rbind(combAll,temp)
  combAll[,MainUnit := gsub("[a-c]$|\\.[1-9]$","",SS.pred)]
  combAll <- combAll[!(BGC == BGC.pred  &  SS_NoSpace != MainUnit),] ### removes overlap where past BGC = future BGC
  combAll[,MainUnit := NULL]
  combAll <- unique(combAll[!is.na(SS_NoSpace),])
  
  
  ##add in BGC probability
  combAll <- na.omit(combAll)
  combAll[,SSratio := allOverlap/sum(allOverlap), by = list(SiteRef, FuturePeriod, BGC, BGC.pred,SS_NoSpace)] ##should check this?
  setorder(combAll, SiteRef, FuturePeriod, BGC, BGC.pred, SS_NoSpace)
  
  setkey(combAll, SiteRef, FuturePeriod, BGC,BGC.pred)
  temp <- unique(combAll[,list(SiteRef,FuturePeriod,BGC,BGC.pred,BGC.prop)])
  temp[,BGC.prop := BGC.prop/sum(BGC.prop), by = list(SiteRef,FuturePeriod,BGC)]
  temp <- unique(temp)
  combAll[,BGC.prop := NULL]
  combAll <- temp[combAll]
  combAll[,SSprob := SSratio*BGC.prop]
  combAll <- combAll[!duplicated(combAll),]
  noPhase <- combAll
  if(onlyRegular){
    return(noPhase)
  }
  
  ###########################################################################
  ##now redo for phases
  numEdaPh <- E1_Phase[,list(NumEdas = .N), by = list(SS_NoSpace)]
  phaseSmall <- unique(edaPhase[,list(BGC,MainUnit,Phase = SS_NoSpace)])
  combPhase <- phaseSmall[noPhase, on = c(MainUnit = "SS.pred"), allow.cartesian = T]
  justPhase <- combPhase[!is.na(Phase),]
  curr <- unique(justPhase[,list(SiteRef, FuturePeriod, BGC = i.BGC, BGC.pred, SS_NoSpace)])
  fut <- unique(justPhase[,list(SiteRef, FuturePeriod, BGC = i.BGC, BGC.pred, MainUnit, Phase)])
  phaseTemp <- E1_Phase[,list(SS_NoSpace,Edatopic)]
  curr <- E1[curr, on = "SS_NoSpace", allow.cartesian = T] 
  fut <- phaseTemp[fut, on = c(SS_NoSpace = "Phase"),allow.cartesian = T]
  setnames(fut,old = "SS_NoSpace", new = "PhasePred")
  setkey(curr, SiteRef, FuturePeriod, BGC, BGC.pred,Edatopic)
  setkey(fut,SiteRef,FuturePeriod,BGC,BGC.pred,Edatopic)
  new <- fut[curr, allow.cartesian = T]
  new <- new[!is.na(PhasePred),]
  
  ###forwards overlap
  SS.out <- new[,list(SS.prob = .N,MainUnit = MainUnit[1]), 
                keyby = list(SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,PhasePred)]
  SS.out[numEda,SS.Curr := i.NumEdas, on = c(SS_NoSpace = "SS_NoSpace"), allow.cartesian = T]
  SS.out[,SSProb := SS.prob/SS.Curr]
  
  ###reverse overlap
  SS.out.rev <- new[,list(SS.prob = .N,MainUnit = MainUnit[1]), 
                    keyby = list(SiteRef,FuturePeriod,BGC,BGC.pred,PhasePred,SS_NoSpace)]
  SS.out.rev[numEdaPh,SS.Curr := i.NumEdas, on = c(PhasePred = "SS_NoSpace"), allow.cartesian = T]
  SS.out.rev[,SSProbRev := SS.prob/SS.Curr]

  ##combine them
  combPhaseAll <- merge(SS.out,SS.out.rev,by = c("SiteRef","FuturePeriod","BGC","BGC.pred","SS_NoSpace","PhasePred"))
  combPhaseAll[,allOverlap := SSProb*SSProbRev]
  combPhaseAll <- combPhaseAll[,list(SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,
                                  SS.pred = MainUnit.y, PhasePred,phaseOverlap = allOverlap)]
  combPhaseAll <- na.omit(combPhaseAll)
  setkey(combPhaseAll,SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)
  setkey(noPhase,SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)
  combAll <- combPhaseAll[noPhase]
  ####add phase to non-phase###
  combAll[,PhaseSum := sum(phaseOverlap), by = list(SiteRef,FuturePeriod,BGC,BGC.pred,SS_NoSpace,SS.pred)]
  combAll[,phaseOverlap := phaseOverlap/PhaseSum]
  combAll[!is.na(phaseOverlap),allOverlap := allOverlap * phaseOverlap]
  combAll[!is.na(phaseOverlap),SS.pred := PhasePred]
  combAll[,c("phaseOverlap","PhaseSum","PhasePred") := NULL]
  ###done phases
  
  # combAll <- merge(combAll,SSsp.out,
  #                  by = c("SiteRef","FuturePeriod","BGC","BGC.pred","SS_NoSpace","SS.pred"), all = T)
  # temp <- combAll[!is.na(allOverlap.y),]
  # temp[,c("allOverlap.x","BGC.prop.x") := NULL]
  # setnames(temp,old = c("allOverlap.y","BGC.prop.y"), new = c("allOverlap","BGC.prop"))
  # combAll[,Flag := if(all(is.na(allOverlap.y))) T else F, by = list(SiteRef,FuturePeriod,BGC,SS_NoSpace,BGC.pred)]
  # combAll <- combAll[(Flag),!"Flag"]
  # combAll[,c("allOverlap.y","BGC.prop.y") := NULL]
  # setnames(combAll,old = c("allOverlap.x","BGC.prop.x"), new = c("allOverlap","BGC.prop"))
  # combAll <- rbind(combAll,temp)
  # combAll[,MainUnit := gsub("[a-c]$|\\.[1-9]$","",SS.pred)]
  # combAll <- combAll[!(BGC == BGC.pred  &  SS_NoSpace != MainUnit),] ### removes overlap where past BGC = future BGC
  # combAll[,MainUnit := NULL]
  # combAll <- unique(combAll[!is.na(SS_NoSpace),])

  
  ##add in BGC probability
  combAll <- na.omit(combAll)
  combAll[,SSratio := allOverlap/sum(allOverlap), by = list(SiteRef, FuturePeriod, BGC, BGC.pred,SS_NoSpace)] ##should check this?
  setorder(combAll, SiteRef, FuturePeriod, BGC, BGC.pred, SS_NoSpace)
  
  setkey(combAll, SiteRef, FuturePeriod, BGC,BGC.pred)
  temp <- unique(combAll[,list(SiteRef,FuturePeriod,BGC,BGC.pred,BGC.prop)])
  temp[,BGC.prop := BGC.prop/sum(BGC.prop), by = list(SiteRef,FuturePeriod,BGC)]
  temp <- unique(temp)
  combAll[,BGC.prop := NULL]
  combAll <- temp[combAll]
  combAll[,SSprob := SSratio*BGC.prop] ##this part may not be correct
  combAll <- combAll[!duplicated(combAll),]
  print("Done EDA")
  return(list(phase = combAll, NoPhase = noPhase))
}
