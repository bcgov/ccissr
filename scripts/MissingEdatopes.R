###Check for missing edatopic space

library(dplyr)
library(readr)
library(tidyr)

#read in edatopic table
edatopic <- read_csv("tables/edatopic.csv")

# Create the full set of edatopic codes (A0 to E7)
letters_part <- LETTERS[1:5]   # A to E
numbers_part <- 0:7            # 0 to 7
full_edatopic <- as.vector(outer(letters_part, numbers_part, paste0))

# Find missing edatopics per BGC
missing_edatopes <- edatopic %>%
  distinct(BGC, Edatopic) %>%
  group_by(BGC) %>%
  summarise(Missing = list(setdiff(full_edatopic, unique(Edatopic)))) %>%
  unnest(Missing)

# View result
print(missing_edatopes)

write_csv(missing_edatopes, "tables/versioned/check_tables/missing_edatopes.csv")
