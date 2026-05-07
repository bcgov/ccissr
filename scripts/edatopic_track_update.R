library(dplyr)
library(compareDF)

#read in current table
edat<-read.csv("tables/edatopic.csv")
edat$X<-NULL
names(edat)

#read in updated table
<<<<<<< HEAD
edat_update<-read.csv("tables/versioned/Edatopic_v13_12.csv") 
=======
edat_update<-read.csv("tables/versioned/Edatopic_v13_18.csv") 
>>>>>>> db3301780b95ec63ae65c8fcc47883a8a727c0eb
edat_update$X<-NULL

#run diff and write out report 
diff <- compareDF::compare_df(edat_update, edat, group_col = c("SS_NoSpace", "Edatopic"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

<<<<<<< HEAD
write.csv(diffrept, "tables/versioned/tracked_changes/edatopic/diff_report_v13_12.csv") #rename with version to match update 
=======
write.csv(diffrept, "tables/versioned/tracked_changes/edatopic/diff_report_v13_18.csv") #rename with version to match update 
>>>>>>> db3301780b95ec63ae65c8fcc47883a8a727c0eb

#overwrite current table with updated table 
#indicate in git commit 
write.csv(edat_update, "tables/edatopic.csv")

#update readme 
path <- "tables/README.md"
<<<<<<< HEAD
updated_readme<- c("# Current table versions ",  "Suitability v13_17 ",        ""               ,
                   "Edatopic v13_12 ",    ""               ,
                   "WNA BGCs v13_2 ",     ""               ,
=======

updated_readme<- c("# Current table versions ",  
                   "Suitability v13_30 ",        ""               ,
                   "Edatopic v13_18 ",    ""               ,
                   "WNA BGCs v13_5 ",     ""               ,
>>>>>>> db3301780b95ec63ae65c8fcc47883a8a727c0eb
                   "Site series v13_2")     #update version #s 
# Write the updated content back to the README
writeLines(updated_readme, path)

