#License info ----
#Copyright 2019 Province of British Columbia
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#libraries 
#library(tidyverse)
library(dplyr)
library(stringr)

#compare the list of BGCs, edatopic space table and suitability table
#Each BGC should have data in both the edatopic and suitability table.
#Each site series listed in the edatopic table should also have species information in the suitability table (and vis versa).

#April 2026 list in CCISS 
BGC_list<- c("MSdm1", "IDFxm", "MSdm2", "MSdm3", "CWHmm1", "IMAun", "CWHmm2", "SBPSxc", "ESSFunp", "IDFdk3", "IDFdk4", "IDFdk1", "ICHmk5", 
             "IDFdk2", "ICHmk3", "ICHmk4", "IDFdk5", "ICHxwa", "MSxv", "ICHmc1","CWHvh3", "ICHmc2", "ICHdw3", "CWHvh1", "SWBvks", "CWHvh2", "ICHdw1", 
             "CWHds1", "ICHdw4", "CWHds2", "CMAwh", "ICHdk", "IDFxx1", "ESSFwmp","IDFmw2", "ESSFvc", "CWHws1", "BAFAun", "CWHws2", "PPxh1", "PPxh2", 
             "ICHdm", "SBSmw", "CWHws3", "SBSmz", "SBSvc", "ESSFdvp", "ESSFwmw", "ESSFmwp", "ESSFvcp", "SBSvk", "ESSFwcw", "ESSFdvw", "ESSFxvp", 
             "SBSwk3a", "SBSmk1", "SBSmk2", "BGxh2", "BGxh3", "ESSFxvw", "IDFdm1", "BGxh1", "IDFdm2", "SBSmc3", "BWBSdk", "SBSmc1", "SBSdw1", "ICHmk1", 
             "SBSmc2", "ICHmk2", "ICHmc1a", "IDFxw", "SBSdw2", "SBSdw3", "MHun1", "ESSFmc", "ESSFun2", "MSxk2", "MSxk3", "MSxk1", "ESSFmh", "MHms", 
             "ESSFmk", "MSdk", "CWHdm2", "CWHdm3", "BWBSmk", "CWHdm1", "ICHwc", "MSdw", "IDFxx2", "ICHxm1", "MSdv", "MHwhp", "MSun1", "ESSFxc1", 
             "ICHwk3", "ESSFxc2", "ICHwk4", "ICHwk1", "SWBmks", "ICHwk2","IDFxh1", "ICHun1", "ICHmm", "ESSFvcw", "SBSdh1", "IDFxh2", "SBSdh2", 
             "ESSFwvp", "ESSFdh2", "ESSFdh1", "ESSFmkp", "ESSFmm2", "ESSFwvw", "ESSFmm3", "MHwh", "ESSFmm1", "ICHvc", "ESSFmcp", "ESSFwh3", 
             "MSdc2", "MSdc3", "ESSFwh1", "ESSFmkw", "ESSFwh2", "MSdc1", "CWHwm","ESSFmcw", "ICHxw", "MHmsp", "ESSFwv", "MHmm1", "ESSFdk1", "ESSFdk2", 
             "MHmm2", "SBSwk1", "SWBun1", "MHvh", "SBSwk2", "ESSFdc1", "ESSFdc2","ESSFdc3", "SWBmk", "ESSFwk2", "SBSwk3", "MHunp", "CDFmm", "ESSFwk1", 
             "MHvhp", "CWHvm1", "BWBSmw", "ICHvk1", "ESSFxcp", "ESSFwc2", "ESSFwc3", "ESSFxc3", "ICHvk2", "BWBSvk", "ESSFmw2", "SBPSmc", 
             "ICHmw1", "ICHmw2", "ICHmw5", "CWHms3", "ICHmw3", "CWHms4", "ICHmw4", "SBPSmk", "ESSFmmp", "CWHms5", "ESSFmmw", "IDFdc", "SBSdk", "ESSFun1", 
             "CWHwh2", "BGxw1", "BGxw2", "CWHwh1", "ESSFmvp", "IDFdh", "ESSFmv3",  "ESSFmv4", "SBSun", "ESSFmv1", "ESSFmv2", "MHmmp", "SBSmh", "SBPSdc", 
             "ESSFdkp", "IDFdw", "SBSmm", "ESSFwm1", "ESSFdcp", "ESSFdkw","ESSFwm4", "ESSFwm2", "CMAun", "ESSFwm3", "ESSFdcw", "ESSFdv1", 
             "IDFww", "CWHvm4", "CWHvm2", "ESSFdv2", "CWHvm3", "ESSFwcp", "ESSFmw1", "IDFxc", "SWBuns", "ESSFwc7", "ESSFwc4", "ESSFxcw", 
             "CWHxs", "ESSFxv1", "BWBSwk3", "BWBSwk1", "BWBSwk2", "IDFxk","ESSFmww", "SWBvk", "ESSFxv2", "BWBScmW", "BWBSnm", "ESSFabN", 
             "ESSFabS", "BWBSlbE", "BWBScmE", "BWBScmC", "BWBSlbN", "BWBSlbW",  "ESSFabC", "BWBSlf", "ESSFabE", "MSabN", "BWBSdmN", "MSabS", 
             "BWBSdmS", "MGPdm", "FGff", "SBAPfp", "IMAab", "BWBScmNW", "BSJPpa", "BWBSub", "BSJPap", "SBAPcp", "BWBSuf", "ESSFdk_MT", #"BWBScm", 
             "MGPmg", "MSdw_MT", "IMAun_MT", "FGnf", "SASbo", "MGPmw_MT", "IDFdk_MT", "BWBSpp", "BSJPku", "OWdm_CA", "CVGdm_CA", "IGFdk_WA", 
             "BGxh_WA", "MSxh_OR", "CWHdm_OR", "ESSFxc_CO", "ICHdw_ID", "PPxw_WY",  "IDFxx_MT", "ESSFxh_WA", "CWHdm_WC", "ICHdh_ID", "ICHxw_WA", 
             "CWHdm_WA", "BGwm_WA", "PPxh_WA", "BGxxx_WY", "CWHdm_OC", "MSdh_WY",  "PPJWmw_UT", "MHms_WA", "PPxx_OR", "ESSFwm_MT", "IDFdh_WA", "CCHun_CA", 
             "BGxx_WY", "CMXmm_OR", "MSmm_ID", "PPxw_MT", "OWmm_OR", "ESSFxk_WY",  "MHds_OR", "ESSFxw_OR", "CWHms_OR", "PPJWdm_WY", "IDFdh_MT", 
             "CMXdm_OR", "CDFmm_OR", "PPxk_OR", "IDFdx_MT", "CWHvh_WA", "IDFxh_ID", "IWFdm_CO", "ESSFxx_WY", "IWFxm_CA", "MHRFmm_OR", 
             "PPJWdw_UT", "ICHmh_MT", "ICHmw_MT", "MSdh_MT", "CRFmh_CA", "IDFdh_UT", "CWHxm_OR", "MHds_WA", "IDFxm_CO", "CMXwm_OR", "IDFdk_WA", "OWxm_CA", 
             "BGdh_OR", "CWHms_WA", "CMXxm_CA", "MSxk_WA", "CWHmh_OR", "PPxm_OR", "ESSFxk_UT", "MSxx_NV", "MSdh_UT", "MHRFdm_OR", "CDFmm_WA", "PPxh_CO", 
             "ESSFxk_MT", "IDFxh_OR", "CWFmm_OR", "WJPxm_OR", "IGFmh_OR", "MDCHun_NV", "CDFxm_CA", "MHRFds_CA", "BGmw_OR", "ICHxh_ID", 
             "BGmw_WA", "CWHxm_WA", "IGFxm_OR", "PPmx_WY", "CWHvm_OR", "PPJWxw_NV", "IDFdxx_WY", "GBDdw_UT", "MSdm_WA", "OWwm_OR", "GBDxx_UT", "ESSFdw_ID", 
             "IDFxh_WA", "IDFdm_WA", "CWFdm_OR", "CRFdh_CA", "CWHds_WA", "BGmk_NV","BGdw_OR", "ESSFxc_WA", "MHmm_WA", "SGPdm_CO", "CWFds_CA", "IWFdw_OR", 
             "CWHvm_WA", "ESSFmw_WA", "CWHwh_OR", "ESSFdh_WA", "ESSFwh_MT", "IDFdm_MT", "MHdm_OR", "ICHxw_ID", "IWFxk_NV", "BGmk_MT", "CWFwm_OR", 
             "ESSFdm_ID", "MMMdk_NV", "BGdw_WA", "MSSDun_NV", "PPJWxh_CA", "PPJWwm_CO", "MSdk_CO", "IGFmm_OR", "GBDxh_NV", "IGFdw_OR", "CWHwh_WA", 
             "BGmm_MT", "IDFxx_WA", "IDFdkx_MT", "MSxm_CO")

#read in tables 
BGC_list_new<-read.csv("tables/WNA_BGCs_Info.csv")
BGC_list2<-unique(BGC_list_new$BGC)

edatopic_tab<-read.csv("tables/edatopic.csv")
suit_tab<-read.csv("tables/suitability.csv")
ss_tab<-read.csv("tables/site_series.csv")

#cross-check BGC list & WNA_BGCs 
check1<-setdiff(BGC_list, BGC_list2) 
check3<-intersect(BGC_list2, BGC_list) #ALL in both 
check2<-setdiff(BGC_list2, BGC_list) #in WNA_BGC but not CCISS list-> can leave as is (alpine and parkland units lumped in mapping)

#cross-check BGC list & edat
check0<-subset(edatopic_tab, BGC %in%BGC_list)
check0<-unique(check0$BGC)
setdiff(BGC_list, check0)

check00<-unique(edatopic_tab$BGC)
sort(setdiff(check00, BGC_list))#in edat but not CCISS list-> can leave as is (alpine and parkland units lumped in mapping)

#cross-check BGC list & suit
check<-subset(suit_tab, bgc %in%BGC_list)
check<-unique(check$bgc)
sort(setdiff(BGC_list, check)) # all of these alpine, grassland, scrub so ok to ignore 


#cross-check BGC list & ss - PAUSED- waiting for new ss CSV from Will 
BGC_list<-as.data.frame(BGC_list)
#remove non BC units- don't have ss names
BGC_list_sub<-filter(BGC_list, !grepl('_|abE|abN|abS|abC|SBAP|SASbo', BGC_list))#US units
AB_units<-c("BWBScmC",  "BWBScmE",  "BWBScmNW" ,"BWBScmW",  "BWBSdmN" , "BWBSdmS" , "BWBSlbE" , "BWBSlbN",  "BWBSlbW", "BSJPap", "BSJPku","BSJPpa",
                      "IMAab", "BWBSlf",   "BWBSnm",   "BWBSpp",   "BWBSub",   "BWBSuf",  "FGff", "FGnf",  "MGPdm","MGPmg")  
BGC_list_sub<-filter(BGC_list_sub, !BGC_list%in% AB_units) #AB units
#turn back into lists
BGC_list<-BGC_list$BGC_list
BGC_list_sub<-BGC_list_sub$BGC_list

check4<-unique(ss_tab$BGC_NoSpace)
xx<-sort(setdiff(BGC_list_sub, check4)) 
xxedat<-subset(edatopic_tab, BGC%in%xx)
sort(setdiff(check4, BGC_list_sub)) #remove/rename

#tree spp codes
treespp<-read.csv("tables/Tree_species_codes.csv")
treespp<-unique(treespp$TreeCode)
treespp2<-unique(suit_tab$sppsplit)
sort(setdiff(treespp2, treespp))

intersect(treespp, treespp2)
