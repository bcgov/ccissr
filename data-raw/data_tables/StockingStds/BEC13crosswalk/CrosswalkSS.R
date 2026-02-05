#Crosswalk stocking standards to BEC v13
#Vanessa Comeau
#June 2025

#set working directory
setwd("C:/Users/vcomeau/Documents/Github/StockingStandards")

# Load necessary package
library(dplyr)

#####CROSSWALK the StockStands table ####

# Read in the CSV files
stockstands12 <- read.csv("StockStands_v12_2.csv", stringsAsFactors = FALSE)
crosswalk <- read.csv("LMH77_BGC_su_Crosswalk_v3.csv", stringsAsFactors = FALSE)

names(read.csv("LMH77_BGC_su_Crosswalk_v3.csv"))
names(read.csv("StockStands_v12_2.csv"))

# Step 1: Create a lookup table
lookup <- crosswalk %>%
  distinct(OldBGC, OldSiteUnitLMH28, NewBGC, NewSiteUnit)

# Find NewBGC + NewSiteUnit that map from multiple OldBGC + OldSiteUnitLMH28
many_to_one_mappings <- lookup %>%
  count(NewBGC, NewSiteUnit) %>%
  filter(n > 1) %>%
  arrange(desc(n))

# View results
print(many_to_one_mappings)

# Optional: Pull all the Old → New mappings involved
many_to_one_details <- lookup %>%
  semi_join(many_to_one_mappings, by = c("NewBGC", "NewSiteUnit"))

# Optional: Write to CSV
write.csv(many_to_one_details, "many_to_one_lookup_mappings.csv", row.names = FALSE)


# Step 2: Join lookup to tableA (brings in NewBGC and NewSiteUnit where match found)
stockstands12_join <- left_join(
  stockstands12,
  lookup,
  by = c("ZoneSubzone" = "OldBGC", "SiteSeries" = "OldSiteUnitLMH28"),
  relationship = "many-to-many"  # explicitly allow many-to-many
)


# Step 3: Replace values if a match was found
stockstands13 <- stockstands12_join %>%
  mutate(
    ZoneSubzone = ifelse(!is.na(NewBGC), NewBGC, ZoneSubzone),
    SiteSeries = ifelse(!is.na(NewSiteUnit), NewSiteUnit, SiteSeries)
  )

# Step 4: Delete the temporary update columns
StockStands13 <- stockstands13 %>%
  select(-NewBGC, -NewSiteUnit)


# Step 5: Identify conflicting Suitability values for same ZoneSubzone + SiteSeries + Species
suitability_conflicts <- StockStands13 %>%
  group_by(ZoneSubzone, SiteSeries, Species) %>%
  filter(n_distinct(Suitability) > 1) %>%
  arrange(ZoneSubzone, SiteSeries, Species)

# View conflicts
print(suitability_conflicts)

# Save to CSV
write.csv(suitability_conflicts, "StockStands13_suitability_conflicts.csv", row.names = FALSE)

# Step 6: Identify conflicting Suitability values for same ZoneSubzone + SiteSeries + Standard + Species
suitability_conflicts_standard <- StockStands13 %>%
  group_by(ZoneSubzone, SiteSeries, Standard, Species) %>%
  filter(n_distinct(Suitability) > 1) %>%
  arrange(ZoneSubzone, SiteSeries, Standard, Species)

# View conflicts
print(suitability_conflicts_standard)

# Optional: Save to CSV for review
#write.csv(suitability_conflicts_standard, "StockStands13_suitability_conflicts_by_standard.csv", row.names = FALSE) #there are none!


# Step 7: Save the updated table
write.csv(StockStands13, "StockStands_v13_1.csv")


####CROSSWALK the Stocking Height table ####

# Read in the main table to update
stock_height <- read.csv("StockingHeight_v12.csv", stringsAsFactors = FALSE)

unique(stock_height$SiteSeries)
unique(stock_height$ZoneSubzone)

# Step 1: edits to the "StockingHeight" table

#Remove symbols after site series
stock_height1 <- stock_height %>%
  mutate(
    SiteSeries = as.character(SiteSeries),
    SiteSeries = iconv(SiteSeries, from = "", to = "UTF-8", sub = ""),  # remove invalid characters like \x87
    SiteSeries = gsub("\\*", "", SiteSeries),                           # remove asterisks
    SiteSeries = gsub("\u0087", "", SiteSeries),                        # explicitly remove \x87
    SiteSeries = trimws(SiteSeries, which = "right")                    # remove trailing spaces
  )

unique(stock_height1$SiteSeries)

#Remove the number 47 that shows up after some subzone variants
stock_height2 <- stock_height1 %>%
  mutate(ZoneSubzone = gsub("47", "", ZoneSubzone))

unique(stock_height2$ZoneSubzone)

#Remove the "(" that appears after IDFdc
stock_height2 <- stock_height2 %>%
  mutate(ZoneSubzone = gsub("IDFdc \\(", "IDFdc", ZoneSubzone))

# Step 2: use the lookup table created above 
print(lookup)

# Step 3: Join lookup to stock_height using old BGC and site unit
stock_height_join <- left_join(
  stock_height2,
  lookup,
  by = c("ZoneSubzone" = "OldBGC", "SiteSeries" = "OldSiteUnitLMH28"),
  relationship = "many-to-many"  # allow many-to-many
)

# Step 4: Replace old values with new values where matches were found
stock_height_updated <- stock_height_join %>%
  mutate(
    ZoneSubzone = ifelse(!is.na(NewBGC), NewBGC, ZoneSubzone),
    SiteSeries = ifelse(!is.na(NewSiteUnit), NewSiteUnit, SiteSeries)
  )

# Step 5: Remove temporary columns
stock_height_updated <- stock_height_updated %>%
  select(-NewBGC, -NewSiteUnit)

# Step 6: Identify conflicting Suitability values for same ZoneSubzone + SiteSeries + Standard + Species
height_conflicts_standard <- stock_height_updated %>%
  group_by(ZoneSubzone, SiteSeries, Standard, Species) %>%
  filter(n_distinct(Height) > 1) %>%
  arrange(ZoneSubzone, SiteSeries, Standard, Species)

# View conflicts
print(height_conflicts_standard)

# Step 7: Save the updated table
write.csv(stock_height_updated, "StockingHeight_v13.csv", row.names = FALSE)


####CROSSWALK the Stocking Info table ####

# Read in the main table to update
stock_info <- read.csv("StockingInfo_v12.csv", stringsAsFactors = FALSE)

unique(stock_info$SiteSeries)
unique(stock_info$ZoneSubzone)

# Step 1: edits to the "StockingHeight" table

#Remove symbols after site series
stock_info1 <- stock_info %>%
  mutate(
    SiteSeries = as.character(SiteSeries),
    SiteSeries = iconv(SiteSeries, from = "", to = "UTF-8", sub = ""),  # remove invalid characters like \x87
    SiteSeries = gsub("\\*", "", SiteSeries),                           # remove asterisks
    SiteSeries = gsub("\u0087", "", SiteSeries),                        # explicitly remove \x87
    SiteSeries = trimws(SiteSeries, which = "right")                    # remove trailing spaces
  )

unique(stock_info1$SiteSeries)

#Remove the number 47 that shows up after some subzone variants
stock_info2 <- stock_info1 %>%
  mutate(ZoneSubzone = gsub("47", "", ZoneSubzone))

unique(stock_info2$ZoneSubzone)

#Remove the "(" that appears after IDFdc
stock_info2 <- stock_info2 %>%
  mutate(ZoneSubzone = gsub("IDFdc \\(", "IDFdc", ZoneSubzone))

# Step 2: use the lookup table created above 
print(lookup)

# Step 3: Join lookup to stock_height using old BGC and site unit
stock_info_join <- left_join(
  stock_info2,
  lookup,
  by = c("ZoneSubzone" = "OldBGC", "SiteSeries" = "OldSiteUnitLMH28"),
  relationship = "many-to-many"  # allow many-to-many
)

# Step 4: Replace old values with new values where matches were found
stock_info_updated <- stock_info_join %>%
  mutate(
    ZoneSubzone = ifelse(!is.na(NewBGC), NewBGC, ZoneSubzone),
    SiteSeries = ifelse(!is.na(NewSiteUnit), NewSiteUnit, SiteSeries)
  )

# Step 5: Remove temporary columns
stock_info_updated <- stock_info_updated %>%
  select(-NewBGC, -NewSiteUnit)

# Step 6: Identify conflicting Suitability values for same ZoneSubzone + SiteSeries + Standard
info_conflicts_standard <- stock_info_updated %>%
  group_by(ZoneSubzone, SiteSeries, Standard) %>%
  filter(n_distinct(StockingTarget) > 1) %>%
  arrange(ZoneSubzone, SiteSeries, Standard)

# View conflicts
print(info_conflicts_standard)

# Step 7: Save the updated table
write.csv(stock_info_updated, "StockingInfo_v13.csv", row.names = FALSE)
