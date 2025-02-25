library(dplyr)
library(compareDF)

#read in current table
feas<-read.csv("feas_tables/feasibility.csv")
feas$X<-NULL
names(feas)

#read in updated table
feas_update<-read.csv("feas_tables/versioned/Feasibility_v13_4.csv") 

#run diff and write out report 
diff <- compareDF::compare_df(feas_update, feas, group_col = c("ss_nospace", "spp"), stop_on_error = FALSE)
#diff$change_summary
#diff$comparison_df

write.csv(diff$comparison_df, "feas_tables/versioned/tracked_changes/diff_report_v13_4.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(feas_update, "feas_tables/feasibility.csv")




              