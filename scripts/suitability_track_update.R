library(dplyr)
library(compareDF)

#read in current table
suit<-read.csv("tables/suitability.csv")
suit$X<-NULL
names(suit)

#read in updated table
suit_update<-read.csv("tables/versioned/suitability_v13_17.csv") 
suit_update$X<-NULL

#run diff and write out report 
diff <- compareDF::compare_df(suit_update, suit, group_col = c("ss_nospace", "spp"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

write.csv(diffrept, "tables/versioned/tracked_changes/suitability/diff_report_v13_17.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(suit_update, "tables/suitability.csv")


#update readme 
path <- "tables/README.md"
updated_readme<- c("# Current table versions ",  "Suitability v13_17 ",        ""               ,
                   "Edatopic v13_11 ",    ""               ,
                   "WNA BGCs v13_2 ",     ""               ,
                   "Site series v13_2")     #update version #s 
# Write the updated content back to the README
writeLines(updated_readme, path)



              