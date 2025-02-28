library(dplyr)
library(compareDF)

#read in current table
suit<-read.csv("tables/suitability.csv")
suit$X<-NULL
names(suit)

#read in updated table
suit_update<-read.csv("tables/versioned/suitability_v13_8.csv") 
suit_update$X<-NULL

#run diff and write out report 
diff <- compareDF::compare_df(suit_update, suit, group_col = c("ss_nospace", "spp"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

write.csv(diffrept, "tables/versioned/tracked_changes/diff_report_v13_8.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(suit_update, "tables/suitability.csv")




              