library(dplyr)
library(compareDF)

#read in current table
site_series<-read.csv("tables/site_series.csv")
names(site_series)
site_series$X<-NULL

#read in updated table
site_series_update<-read.csv("tables/versioned/SiteSeries_v13_4.csv") 
names(site_series_update)
cols<-names(site_series_update)

#only keep columns in updated df
site_series<-site_series %>%select(all_of(cols))

#run diff and write out report 
diff <- compareDF::compare_df(site_series_update, site_series, group_col = c("SS_NoSpace"), stop_on_error = FALSE)
diff$change_summary
diffrept<-diff$comparison_df

write.csv(diffrept, "tables/versioned/tracked_changes/SiteSeries/diff_report_v13_4.csv") #rename with version to match update 

#overwrite current table with updated table 
#indicate in git commit 
write.csv(site_series_update, "tables/site_series.csv")


#update readme - make sure all are correct!!
path <- "tables/README.md"
updated_readme<- c("# Current table versions ",  
                   "Suitability v13_30",        ""               ,
                   "Edatopic v13_18 ",    ""               ,
                   "WNA BGCs v13_9 ",     ""               ,
                   "Site series v13_4")     #update version #s 
# Write the updated content back to the README
writeLines(updated_readme, path)



