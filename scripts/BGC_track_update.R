library(dplyr)
library(compareDF)

#read in current table
BGC_list<-read.csv("tables/WNA_BGCs_Info.csv")
BGC_list$X<-NULL
#remove columns no longer using
#BGC_list$Retire<-NULL
#BGC_list$SS_class<-NULL
#BGC_list$BGC_Label<-NULL

#read in updated table
BGC_update<-read.csv("tables/versioned/WNA_BGCs_v13_8.csv")
BGC_update$X<-NULL

#run diff and write out report 
diff <- compareDF::compare_df(BGC_update, BGC_list, group_col = c("BGC", "BGCName"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

write.csv(diffrept, "tables/versioned/tracked_changes/BGCs/diff_report_v13_8.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(BGC_update, "tables/WNA_BGCs_Info.csv")


#update readme 
path <- "tables/README.md"
updated_readme<- c("# Current table versions ",  
                   "Suitability v13_30",        ""               ,
                   "Edatopic v13_17 ",    ""               ,
                   "WNA BGCs v13_8 ",     ""               ,
                   "Site series v13_2")     #update version #s 
# Write the updated content back to the README
writeLines(updated_readme, path)


