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
#compare the list of BGCs, edatopic space table and feasibility table
#Each BGC should have data in both the edatopic and feasibility table.
#Each site series listed in the edatopic table should also have species information in the feasibility table (and vis versa).

#read in tables 
edatop_tab<-read.csv("tables/versioned/Edatopic_v13_10.csv")#update to local if needed
feas_tab<-read.csv("tables/versioned/suitability_v13_15.csv")
BGC_list<-read.csv("tables/versioned/BGCs_V13_22Feb2025.csv")

subzones<-as.data.frame(unique(BGC_list$BGC))#427
subzones$BGC<-subzones$`unique(BGC_list$BGC)`
subzones$`unique(BGC_list$BGC)`<-NULL


ss_e<-select(edatop_tab, BGC, SS_NoSpace)%>%rename(ss_nospace=SS_NoSpace)%>%distinct(.)
ss_f<-select(feas_tab, bgc, ss_nospace)%>%rename(BGC=bgc)%>%distinct(.)

miss_BGCs_feas<-anti_join(subzones, ss_f)
miss_SS_feas<-anti_join(ss_e, ss_f, by = 'ss_nospace')

miss_BGCs_edatope<-anti_join(subzones, ss_e)
miss_SS_edatope<-anti_join(ss_f, ss_e, by = 'ss_nospace')

write.csv(miss_BGCs_edatope, "tables/versioned/check_tables/missing_BGCs_edatope4.csv")
write.csv(miss_BGCs_feas, "tables/versioned/check_tables/missing_BGCs_feasibility4.csv")
write.csv(miss_SS_feas, "tables/versioned/check_tables/missing_SS_feasibility4.csv")
write.csv(miss_SS_edatope, "tables/versioned/check_tables/missing_SS_edatope4.csv")
