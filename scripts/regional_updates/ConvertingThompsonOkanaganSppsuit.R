library(tidyr)
library(dplyr)
library(stringr)

#converting Thompson Okanagan spp suitability ratings from columns to rows

# Load the CSV file
data <- read.csv("tables/regional_updates/TOSppSuit.csv", stringsAsFactors = FALSE)

# Rename the first column from BGC to bgc
colnames(data)[colnames(data) == "BGC"] <- "bgc"

# Combine bgc and SiteSeries columns with a new title
data$ss_nospace <- paste(data$bgc, data$SiteSeries, sep = "/")

# Remove SiteSeries column
data$SiteSeries <- NULL

# Convert from wide to long format using correct column names
long_data <- pivot_longer(data, cols = c("X1", "X2", "X3", "X4", "X5"),
                          names_to = "newsuit", values_to = "Species")

# Remove rows with NA or empty Species values
long_data <- long_data[!(is.na(long_data$Species) | long_data$Species == ""), ]

# Convert newsuit column from X1, X2, X3, X4, X5 to numerical values 1, 2, 3, 4, 5
long_data$newsuit <- as.numeric(sub("X", "", long_data$newsuit))

# Remove all numeric values from Species column, keeping only letters (remove footnotes keep species)
long_data$Species <- gsub("[0-9]", "", long_data$Species)

# Separate each species (two-letter codes) into separate rows
long_data <- long_data %>% mutate(Species = str_extract_all(Species, "[A-Z][a-z]")) %>% unnest(Species)

# Rename Species column to spp
colnames(long_data)[colnames(long_data) == "Species"] <- "spp"

# Update Se to Sx
long_data$spp[long_data$spp == "Se"] <- 'Sx'

# Add a new column 'mod' with every value set to 'MR'
long_data$mod <- "MR"

# Add a column for OHR, for Mike's units Lw = TRUE, all other spp = FALSE
long_data$outrange <- ifelse(long_data$spp == "Lw", TRUE, FALSE)

# Convert species value from 'Xx' to 'X'
long_data$spp <- gsub("Xx", "X", long_data$spp)

# Set newsuit to NA where spp is 'X'
long_data$newsuit[long_data$spp == "X"] <- NA

# Reorder columns so that spp comes before newsuit
long_data <- long_data %>% select(bgc, ss_nospace, spp, newsuit, mod, outrange, everything())

# Print preview of transformed data
head(long_data)

#remove variants - Mike's decisions March 2025
#ICHdw3/06 Use the Cw-dominated unit.
#ICHdw3/08 Use the mineral soil with horsetail unit
#ICHmk2/05 Use the Cw-dominated unit
#ICHmw3/08 Use the 08 (organic soils with skunk cabbage)
#IDFdc/03 Use the03 (shallow soils)
#IDFmw2/04 Use lacks abundant devil’s club
#MSdc1/02 Use the 02 (low elevations)
#MSdm3/03 Use 03 shallow soils
#MSxk3/05 steep warm slopes

long_data<-subset(long_data, ss_nospace!= "ICHdw3/06 (Sx present)"&
                            ss_nospace!= "ICHdw3/08 (organic soils with skunk cabbage)" &
                              ss_nospace!= "ICHmk2/05 (Sx-dominant)" &
                              ss_nospace!= "ICHmw3/08 (mineral soils with horsetail)"&   
                              ss_nospace!= "IDFdc/03 (very steep slopes with pinegrass)"&
                              ss_nospace!= "IDFdc/03 (very steep slopes with bluebunch wheatgrass)"&
                              ss_nospace!= "IDFmw2/04 (abundant devil's club present)"  &
                              ss_nospace!= "MSdc1/02 (high elevations)"  &
                              ss_nospace!= "MSdm3/03 (deep soils)"  &   
                              ss_nospace!= "MSxk3/05 (moderate and gentle slopes)")  

# Remove all variant info from ss column, 
long_data<-separate(long_data, ss_nospace, into=c("ss", "remove"), sep = " " , remove=F)%>%
  select(-ss_nospace, -remove)%>%rename(ss_nospace=ss)

# Save the transformed data
write.csv(long_data, "tables/regional_updates/Thompson_Okanagan_Mar2025.csv", row.names = FALSE)




