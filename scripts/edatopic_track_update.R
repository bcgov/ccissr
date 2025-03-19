library(dplyr)
library(compareDF)

#read in current table
edat<-read.csv("tables/edatopic.csv")
edat$X<-NULL
names(edat)

#read in updated table
edat_update<-read.csv("tables/versioned/Edatopic_v13_9.csv") 
edat_update$X<-NULL

#run diff and write out report 
diff <- compareDF::compare_df(edat_update, edat, group_col = c("SS_NoSpace", "Edatopic"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

write.csv(diffrept, "tables/versioned/tracked_changes/edatopic/diff_report_v13_9.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(edat_update, "tables/edatopic.csv")

#update readme 
path <- "tables/README.md"
updated_readme<- c("# Current table versions ",  "Suitability v13_13 ",        ""               ,
                   "Edatopic v13_9 ",    ""               ,
                   "WNA BGCs v13_2 ",     ""               ,
                   "Special sites v13_2")     #update version #s 
# Write the updated content back to the README
writeLines(updated_readme, path)

