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

#tree spp codes
treespp<-read.csv("tables/metadata/Tree_species_codes.csv")
treespp<-unique(treespp$TreeCode)
treespp2<-unique(suit_tab$sppsplit)
sort(setdiff(treespp2, treespp))

intersect(treespp, treespp2)


#update newest ss list from Will - Aug 2026
site_series_new<-read.csv("tables/versioned/SiteSeries_v13_4.5.csv") 
site_series<-read.csv("tables/site_series.csv")
site_series$X<-NULL

site_series_new<-left_join(select(site_series, -SiteSeriesLongName), site_series_new, by="SS_NoSpace")%>%distinct(.)

site_series_new<-mutate(site_series_new, SpecialCode2=case_when(grepl('/Wb', SS_NoSpace) ~"Wb", 
                                                               grepl('/Wm', SS_NoSpace) ~"Wm", 
                                                               grepl('/Ws', SS_NoSpace) ~"Ws", 
                                                               grepl('/Wf', SS_NoSpace) ~"Wf", 
                                                               grepl('/Wa', SS_NoSpace) ~"Wa", 
                                                               grepl('/Wb', SS_NoSpace) ~"Wb", 
                                                               grepl('/Fl', SS_NoSpace) ~"Fl", 
                                                               grepl('/Fm', SS_NoSpace) ~"Fm", 
                                                               grepl('/Ff', SS_NoSpace) ~"Ff", 
                                                               grepl('/Gg', SS_NoSpace) ~"Gg", 
                                                               grepl('/Gb', SS_NoSpace) ~"Gb", 
                                                               grepl('/Ga', SS_NoSpace) ~"Ga", 
                                                               grepl('/Gs', SS_NoSpace) ~"Gs", 
                                                               grepl('/Ro', SS_NoSpace) ~"Ro",
                                                               grepl('/Rt', SS_NoSpace) ~"Rt",
                                                               grepl('/Sc', SS_NoSpace) ~"Sc",
                                                               grepl('/Sk', SS_NoSpace) ~"Sk",
                                                               grepl('/Vs', SS_NoSpace) ~"Vs",
                                                               grepl('/Vh', SS_NoSpace) ~"Vh",
                                                               grepl('/Ag', SS_NoSpace) ~"Ag",
                                                               grepl('/Ah', SS_NoSpace) ~"Ah",
                                                               grepl('/Bb', SS_NoSpace) ~"Bb",TRUE~NA))
site_series_new<-mutate(site_series_new, SpecialCode= if_else(SpecialCode=="", NA, SpecialCode))
site_series_new<-mutate(site_series_new, Special= if_else(Special=="", NA, Special))
site_series_new<-mutate(site_series_new, SpecialCode= if_else(is.na(SpecialCode), SpecialCode2, SpecialCode))
site_series_new$SpecialCode2<-NULL

special<-read.csv("tables/metadata/Tr68_special_codes.csv")
site_series_new<-left_join(site_series_new, special)
site_series_new<-mutate(site_series_new, Special= if_else(is.na(Special), Special2, Special))
site_series_new$Special2<-NULL

site_series_new$SpecialCode[site_series_new$SS_NoSpace=="CWHvh2/120"] <- "SS"
site_series_new$SpecialCode[site_series_new$SS_NoSpace=="CWHvh2/120.1"] <- "SS"
site_series_new$SpecialCode[site_series_new$SS_NoSpace=="CWHvh2/120.2"] <- "SS"
site_series_new$Special[site_series_new$SpecialCode=="SS"] <- "salt spray"
site_series_new<-distinct(site_series_new)

write.csv(site_series_new, "tables/versioned/SiteSeries_v13_5.csv") 

#cross-check BGC list & ss 
ss_tab<-read.csv("tables/versioned/SiteSeries_v13_5.csv") 
ss_tab$X<-NULL
BGC_list<-as.data.frame(BGC_list_new)

ss_tab<-tidyr::separate(ss_tab, SS_NoSpace, into=c("BGC"), sep="/", remove=F)

#remove non BC units- don't have ss names
BGC_list_sub<-subset(BGC_list, !grepl('_|abE|abN|abS|abC|SBAP|SASbo', BGC))#US units
AB_units<-c("BWBScmC",  "BWBScmE",  "BWBScmNW" ,"BWBScmW",  "BWBSdmN" , "BWBSdmS" , "BWBSlbE" , "BWBSlbN",  "BWBSlbW", "BSJPap", "BSJPku","BSJPpa",
            "IMAab", "BWBSlf",   "BWBSnm",   "BWBSpp",   "BWBSub",   "BWBSuf",  "FGff", "FGnf",  "MGPdm","MGPmg")  
BGC_list_sub<-subset(BGC_list_sub, !BGC%in% AB_units) #remove AB units

#turn back into lists
BGC_list<-BGC_list$BGC
BGC_list_sub<-BGC_list_sub$BGC
# check what's missing 
check4<-unique(ss_tab$BGC)
xx<-sort(setdiff(BGC_list_sub, check4)) 

#Add missing site series for these units (from Kiri)
missing<-as.data.frame(c("BAFAun/00", "BAFAunp/00", "BGxh1/Gg", "BGxh2/Fm05", "BGxh2/Gg", "BGxh2/Gs", "BGxh3/Fm01", "BGxh3/Gg", 
  "BGxh3/Ro", "BGxw1/Gg", "BGxw1/Ro", "BGxw2/Gg", "BGxw2/Ro", "BWBSdk/Gb", "BWBSdk/Gg", "BWBSdk/Ro", "BWBSmk/Gg", 
  "BWBSmk/Ro", "BWBSmw/Gg", "BWBSmw/Ro", "BWBSvk/h", "BWBSvk/m", "BWBSvk/x", "BWBSwk1/Ro", "BWBSwk2/Gg", "BWBSwk2/Ro", 
  "BWBSwk3/Ro", "CDFmm/Fm51", "CDFmm/Fm52", "CDFmm/Gg", "CDFmm/Ro", "CDFmm/Wf", "CDFmm/Wm", "CMAun/00", "CMAunp/00", 
  "CMAwh/00", "CWHdm1/Fm51", "CWHdm1/Fm52", "CWHdm1/Ro", "CWHdm1/Wsxx", "CWHdm2/Fm51", "CWHdm2/Fm52", "CWHdm2/Ro", 
  "CWHdm2/Ws", "CWHdm3/Fl50", "CWHdm3/Fm51", "CWHdm3/Ro", "CWHds1/Fm51.1", "CWHds1/Ro", "CWHds1/Wb50.2", "CWHds1/Wm", 
  "CWHds1/Ws56.2", "CWHds2/Fl", "CWHds2/Fm", "CWHds2/Ro", "CWHmm1/Fl50", "CWHmm1/Fm51", "CWHmm1/Ro", "CWHmm1/Ws##",
  "CWHmm2/Ro", "CWHms3/Fm51.1", "CWHms3/Wb50.2", "CWHms3/Wf", "CWHms3/Wm", "CWHms4/Wb50.2", "CWHms4/Ws56.1", 
  "CWHms5/Fm51.1", "CWHms5/Wb50.2", "CWHms5/Wm", "CWHvh1/Bb", "CWHvh1/Br", "CWHvh1/Em", "CWHvh1/Fm52", "CWHvh2/Bb",
  "CWHvh2/Br", "CWHvh2/Em", "CWHvh2/Fm52", "CWHvm1/Fl50", "CWHvm1/Fm52", "CWHvm1/Wf", "CWHvm1/Wm", "CWHvm2/Wf", 
  "CWHvm2/Wm", "CWHvm3/Fl50", "CWHvm3/Wf", "CWHvm3/Wm", "CWHvm4/Wf", "CWHvm4/Wm", "CWHwh2/Ro", "CWHwm/Fl50", 
  "CWHwm/Wf", "CWHwm/Wm", "CWHws1/Fl50", "CWHws1/Ro", "CWHws2/Ro", "CWHws3/Fl50", "CWHxs/Fl50", "CWHxs/Fm51.2",
  "CWHxs/Ro", "CWHxs/Ws57", "ESSFdc1/Ro", "ESSFdc2/Ro", "ESSFdcp/00", "ESSFdcw/Ro", "ESSFdh2/01", "ESSFdh2/02",
  "ESSFdh2/03", "ESSFdh2/04", "ESSFdh2/05", "ESSFdh2/06", "ESSFdh2/07", "ESSFdh2/08", "ESSFdk1/Ro", "ESSFdk2/Ro",
  "ESSFdkp/00", "ESSFdkw/Ro", "ESSFdv1/01", "ESSFdv1/02", "ESSFdv1/03", "ESSFdv1/04", "ESSFdv1/05", "ESSFdv1/06", 
  "ESSFdv2/01", "ESSFdv2/02", "ESSFdv2/03", "ESSFdv2/04", "ESSFdv2/05", "ESSFdv2/06", "ESSFdvp/00", "ESSFdvw/01", 
  "ESSFdvw/02", "ESSFdvw/03", "ESSFdvw/04", "ESSFdvw/05", "ESSFdvw/06", "ESSFmcp/00", "ESSFmcw/Ro", "ESSFmh/112.1", 
  "ESSFmh/112.2", "ESSFmh/Ro", "ESSFmkp/00", "ESSFmkw/Ro", "ESSFmm1/Ro", "ESSFmm2/01", "ESSFmm2/02", "ESSFmm2/03", 
  "ESSFmm2/04", "ESSFmm2/05", "ESSFmm2/06", "ESSFmm2/07", "ESSFmm2/Ro", "ESSFmm3/Ro", "ESSFmmp/00", "ESSFmmw/01", 
  "ESSFmmw/02", "ESSFmmw/03", "ESSFmmw/04", "ESSFmmw/05", "ESSFmmw/06", "ESSFmmw/07", "ESSFmmw/Ro", "ESSFmv1/Ro", 
  "ESSFmv2/Ro", "ESSFmv3/Ro", "ESSFmv4/Ro", "ESSFmvp/00", "ESSFmw2/01", "ESSFmw2/02", "ESSFmw2/03", "ESSFmw2/04", 
  "ESSFmw2/05", "ESSFmw2/06", "ESSFmw2/07", "ESSFmw2/08", "ESSFmwp/00", "ESSFmww/101", "ESSFmww/102", "ESSFmww/103",
  "ESSFmww/110", "ESSFunp/00", "ESSFvcp/00", "ESSFvcw/01", "ESSFvcw/02", "ESSFvcw/03", "ESSFvcw/04", "ESSFvcw/05", 
  "ESSFvcw/06", "ESSFwc2/Ro", "ESSFwc3/101", "ESSFwc3/102", "ESSFwc3/103", "ESSFwc3/110", "ESSFwc3/111", "ESSFwc3/01", 
  "ESSFwc3/02", "ESSFwc3/03", "ESSFwc3/Ro", "ESSFwc4/Ro", "ESSFwcp/00", "ESSFwcw/Ro", "ESSFwh1/Ro", "ESSFwh2/Ro", 
  "ESSFwh3/Ro", "ESSFwk1/Ro", "ESSFwk2/Ro", "ESSFwm1/Ro", "ESSFwm2/Ro", "ESSFwm3/Ro", "ESSFwm4/Ro", "ESSFwmp/00", 
  "ESSFwmw/Ro", "ESSFwvp/00", "ESSFwvw/Ro", "ESSFxc1/Gg", "ESSFxc2/Gg", "ESSFxcp/00", "ESSFxcw/101", "ESSFxcw/102",
  "ESSFxcw/103", "ESSFxcw/110", "ESSFxcw/111", "ESSFxcw/112", "ESSFxcw/Gg", "ESSFxv2/Ro", "ESSFxvp/00", "ESSFxvw/01", 
  "ESSFxvw/02", "ESSFxvw/03", "ESSFxvw/04", "ESSFxvw/05", "ESSFxvw/06", "ESSFxvw/07", "ESSFxvw/08", "ESSFxvw/09", 
  "ESSFxvw/10", "ICHdk/Ro", "ICHdm/Ro", "ICHdw1/101", "ICHdw1/102", "ICHdw1/103", "ICHdw1/104", "ICHdw1/110", 
  "ICHdw1/111", "ICHdw1/112", "ICHdw1/Ro", "ICHdw3/01", "ICHdw3/02", "ICHdw3/03", "ICHdw3/04", "ICHdw3/05",
  "ICHdw3/06", "ICHdw3/07", "ICHdw3/08", "ICHdw3/09", "ICHdw3/Ro", "ICHdw4/Ro", "ICHmc1/Ro", "ICHmc1a/Ro",
  "ICHmc2/Ro", "ICHmk1/Ro", "ICHmk3/Ro", "ICHmk4/Ro", "ICHmk5/Ro", "ICHmm/Ro", "ICHmw1/Ro", "ICHmw2/Ro", 
  "ICHmw4/Ro", "ICHmw5/Ro", "ICHvc/Ro", "ICHvk2/Ro", "ICHwc/Ro", "ICHwk1/Ro", "ICHwk2/Ro", "ICHwk3/Ro", 
  "ICHwk4/Ro", "ICHxm1/Ro", "ICHxw/Ro", "ICHxwa/Ro", "IDFdc/01", "IDFdc/02", "IDFdc/03", "IDFdc/04", 
  "IDFdc/05", "IDFdc/06", "IDFdc/07", "IDFdc/08", "IDFdc/Gg", "IDFdh/Ro", "IDFdk1/Gg1", "IDFdk1/Gg2",
  "IDFdk1/Ro", "IDFdk2/Gg", "IDFdk3/Ro", "IDFdk4/Ro", "IDFdk5/Ro", "IDFdm1/Ro", "IDFdm2/Ro", "IDFww/Fl",
  "IDFww/Fm01", "IDFww/Fm02", "IDFww/Ro", "IDFww/Wb", "IDFww/Ws", "IDFxc/01", "IDFxc/02", "IDFxc/03", "IDFxc/04", 
  "IDFxc/05", "IDFxc/06", "IDFxc/07", "IDFxc/08", "IDFxh1/Gg", "IDFxh2/Fm02", "IDFxk/Ro", "IDFxm/Gg", "IDFxm/Ro", 
  "IDFxw/Ro", "IDFxx1/Ro", "IDFxx2/Ro", "IMAun/00", "IMAunp/00", "MHmm1/Wb5X", "MHmm1/Ws55.1", "MHmm2/Wb5X", 
  "MHmm2/Ws55.1", "MHmmp/h", "MHmmp/m", "MHmmp/x", "MHmsp/h", "MHmsp/m", "MHmsp/x", "MHunp/00", "MHvh/Wb", "MHvh/Wf", 
  "MHvhp/h", "MHvhp/m", "MHvhp/x", "MHwhp/00", "MSdc3/01", "MSdc3/02", "MSdc3/03", "MSdc3/04", "MSdc3/05", "MSdk/Ro", 
  "MSdm1/Ro", "MSdm2/Rt08", "MSdw/Ro", "MSxk1/Gg", "MSxk1/Gs", "MSxk2/Gg", "MSxk2/Rt08", "MSxv/Ro", "PPxh1/Gg", 
  "PPxh1/Gs", "PPxh1/Ro", "PPxh2/Gg", "SBPSdc/Ro", "SBPSmc/Ro", "SBPSmk/Ro", "SBPSxc/Ro", "SBSdh1/Ro", "SBSdh2/01",
  "SBSdh2/02", "SBSdh2/03", "SBSdh2/04", "SBSdh2/05", "SBSdh2/06", "SBSdh2/07", "SBSdh2/08", "SBSdh2/Ro", "SBSdk/Gg", 
  "SBSdk/Ro", "SBSdw1/Ro", "SBSdw2/Ro", "SBSdw3/Ro", "SBSmc1/Ro", "SBSmc2/Ro", "SBSmc3/Ro", "SBSmh/Ro", "SBSmk1/Ro",
  "SBSmk2/Ro", "SBSmm/Ro", "SBSmw/Ro", "SBSun/h", "SBSun/m", "SBSun/x", "SBSvk/Ro", "SBSwk1/Ro", "SBSwk2/Ro", 
  "SBSwk3/Ro", "SBSwk3a/Ro", "SWBmks/00", "SWBuns/00", "SWBvk/h", "SWBvk/m", "SWBvk/x", "SWBvks/00", "MHun1/h", 
  "MHun1/m", "MHun1/x", "SWBun1/h", "SWBun1/m", "SWBun1/x", "ESSFun1/m", "ESSFun1/h", "ESSFun1/x", "ESSFun2/m", 
  "ESSFun2/x", "ESSFun2/h", "ICHun1/m", "ICHun1/x", "ICHun1/h", "MSun1/m", "MSun1/x", "MSun1/h"))
names(missing)<-"SS_NoSpace"

#add special codes for these
missing<-tidyr::separate(missing, SS_NoSpace, into=c("BGC"), sep="/", remove=F)
missing$SiteSeriesLongName<-""
missing<-mutate(missing, SpecialCode=case_when(grepl('/Wb', SS_NoSpace) ~"Wb", 
                                                                grepl('/Wm', SS_NoSpace) ~"Wm", 
                                                                grepl('/Ws', SS_NoSpace) ~"Ws", 
                                                                grepl('/Wf', SS_NoSpace) ~"Wf", 
                                                                grepl('/Wa', SS_NoSpace) ~"Wa", 
                                                                grepl('/Wb', SS_NoSpace) ~"Wb", 
                                                                grepl('/Fl', SS_NoSpace) ~"Fl", 
                                                                grepl('/Fm', SS_NoSpace) ~"Fm", 
                                                                grepl('/Ff', SS_NoSpace) ~"Ff", 
                                                                grepl('/Gg', SS_NoSpace) ~"Gg", 
                                                                grepl('/Gb', SS_NoSpace) ~"Gb", 
                                                                grepl('/Ga', SS_NoSpace) ~"Ga", 
                                                                grepl('/Gs', SS_NoSpace) ~"Gs", 
                                                                grepl('/Ro', SS_NoSpace) ~"Ro",
                                                                grepl('/Rt', SS_NoSpace) ~"Rt",
                                                                grepl('/Sc', SS_NoSpace) ~"Sc",
                                                                grepl('/Sk', SS_NoSpace) ~"Sk",
                                                                grepl('/Vs', SS_NoSpace) ~"Vs",
                                                                grepl('/Vh', SS_NoSpace) ~"Vh",
                                                                grepl('/Ag', SS_NoSpace) ~"Ag",
                                                                grepl('/Ah', SS_NoSpace) ~"Ah",
                                                                grepl('/Bb', SS_NoSpace) ~"Bb",TRUE~NA))
special<-read.csv("tables/metadata/Tr68_special_codes.csv")
missing<-left_join(missing, special)%>%rename(Special=Special2)

#bind back to main table
ss_tab<-rbind(ss_tab, missing)

write.csv(ss_tab, "tables/versioned/SiteSeries_v13_6.csv") 
