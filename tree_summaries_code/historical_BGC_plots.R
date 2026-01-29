# Goal: create historical BGC pie chart, suitability maps for each edatope, and (RETIRED) stacked bar charts

#NOTE from Mia: This code is surely not as efficient as it could be - I'm relatively new to writing functions and am still wrapping my head around how the duckdb database querying works (as well as the data and how its structured within duckdb), so much of this code is copied from Kiri's scripts and modified to the best of my current ability (and with some AI "help"). Given this, I KNOW it could be much tighter, but I'm not sure where...and at least I'm pretty confident that it's producing correct results. Please feel welcome to help me improve this workflow! Thanks!


library(tidyverse)
library(duckdb)
library(data.table)
library(DBI)
library(here)



#Run through the CCISS_Plot_Workflow script first!

# Set up for pies ----------------------

# Function to make cciss_sum a dataframe
make_cciss_sum <- function(con, spp, edatope, fractional = TRUE, by_zone = TRUE
) {
  
  # Pull data from DuckDB
  # collapse edatope vector into SQL-friendly string - is this necessary?
  edatope_sql <- paste(sprintf("'%s'", edatope), collapse = ", ")
  
  cciss_spp <- dbGetQuery(
    con,
    sprintf(
      "SELECT * FROM cciss_res 
     WHERE Spp = '%s' 
     AND Edatope IN (%s)",
      spp, edatope_sql
    )
  ) |> as.data.table()
  
  bgc_mapped <- dbGetQuery(con, "SELECT * FROM bgc_points") |> 
    as.data.table()
  
  # Define zone or subzone
  if (by_zone) {
    bgc_mapped[, zone := regmatches(bgc, regexpr("^[A-Z]+", bgc))]
  } else {
    bgc_mapped[, zone := bgc]
  }
  
  # Join zone onto cciss data
  cciss_spp[bgc_mapped, zone := i.zone, on = c(SiteRef = "cellnum")]
  
  # Suitability
  if (fractional) {
    cciss_spp[, SppSuit := Newsuit]
    cciss_spp[is.na(SppSuit) | SppSuit > 3.5, SppSuit := 5]
    cciss_spp[, SppSuit := 1 - (SppSuit - 1) / 4]
    
    cciss_spp[, HistSuit := Curr]
    cciss_spp[is.na(HistSuit) | HistSuit > 3.5, HistSuit := 5]
    cciss_spp[, HistSuit := 1 - (HistSuit - 1) / 4]
  } else {
    cciss_spp[, SppSuit := as.integer(Newsuit <= 3)]
    cciss_spp[, HistSuit := as.integer(Curr <= 3)]
  }
  
  # Summarise
  cciss_sum <- cciss_spp[
    , .(SppArea = sum(SppSuit)),
    by = .(zone, FuturePeriod)
  ]
  
  # Add historical baseline...?
  cciss_hist <- cciss_spp[
    FuturePeriod == "2021_2040",
    .(SppArea = sum(HistSuit)),
    by = .(zone, FuturePeriod)
  ][, FuturePeriod := "1961"]
  
  cciss_sum <- rbind(cciss_sum, cciss_hist)
  
  # Year column
  cciss_sum[, Year := factor(substr(FuturePeriod, 1, 4))]
  
  return(cciss_sum)
}

#color palette
colScheme <- c(PP = "#ea7200", MH = "#6f2997", SBS = "#2f7bd2", ESSF = "#ae38b8", 
               CWH = "#488612", BWBS = "#4f54cf", CWF = "#7577e7", IGF = "#77a2eb", 
               CMX = "#71d29e", BG = "#dd1320", IDF = "#e5d521", MS = "#e44ebc", 
               SWB = "#a1dbde", CRF = "#af3a13", WJP = "#73330e", ICH = "#1fec26", 
               CDF = "#edf418", JPW = "#96b3a5", CMA = "#eae1ee", SBPS = "#6edde9", 
               IMA = "#e3f1fa", GBD = "#4d433f", OW = "#582511", BAFA = "#eee4f1", 
               MMM = "#FF00FF", MHRF = "#2612dc", MGP = "#f0aeab", FG = "#92696c", 
               SGP = "#cca261", GO = "#f0a325", SBAP = "#51d5a7", IWF = "#d44273", 
               BSJP = "#424160", MSSD = "#dac370", MDCH = "#2d0cd4", CVG = "#c9edd3", 
               SAS = "#92b1b6", CCH = "#7e22ca")





# Create dataframe and pies for each species -------------

# Use function to make dataframe for each species
cciss_sum <- make_cciss_sum(
  con = con,
  spp = "Fd", #CHANGE SPP HERE
  edatope = c("C4", "B2", "D6"),
  fractional = TRUE,
  by_zone = TRUE
)

# Keep only historical Year = 1961, and other tidying
cciss_sum <- cciss_sum |>
  filter(Year == "1961") |>
  mutate(SppArea = round(SppArea, 2)) |>
  filter(SppArea > 0) |>
  arrange(SppArea) |>
  mutate(zone = factor(zone, levels = zone))

# Pie chart
pie <- ggplot(cciss_sum, aes(x="", y=SppArea, fill=zone)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values = colScheme,
                    guide = guide_legend(reverse = TRUE)) +
  labs(fill = "BGC zone", 
       title = "Historical (1961-1990) range of Fd in British Columbia"); pie #CHANGE SPP HERE

ggsave("Fd_histpie.png", # CHANGE SPP HERE
plot = pie, path = here("tree_summaries_code", "Figures"),width = 6, height = 4, dpi = 300, bg = "white")













# Historical suitability maps -----------------

plot_historical_suit_map <- function(dbCon, spp, edatope, bgc_template, outline, save_png = TRUE, out_dir = file.path("tree_summaries_code", "Figures"), width = 6, height = 5, res = 300) 
{
  
  # Output
  if (save_png) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_file <- file.path(out_dir, paste("Historic_Suit_Map", spp, edatope, "png", sep = "."))
    png(filename = out_file, type = "cairo", units = "in", width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  
  # Query & prep
  dat_spp <- DBI::dbGetQuery(
    dbCon,
    sprintf(
      "SELECT SiteRef, Curr
       FROM cciss_res
       WHERE Spp = '%s'
         AND Edatope = '%s'",
      spp, edatope
    )
  ) |> data.table::as.data.table()
  
  # Classify Curr like the original workflow
  dat_spp[, Curr := as.integer(round(Curr))]
  dat_spp[is.na(Curr) | Curr > 3.5, Curr := 5]  #5 = "not mapped"/white
  

  # Build raster
  X <- terra::rast(bgc_template$bgc_rast)
  terra::values(X) <- NA
  
  X[dat_spp$SiteRef] <- dat_spp$Curr
  
  # Plot
  par(mar = c(0.5, 0.5, 0.5, 0.5), bty = "n", pty = "m")
  
  breakseq  <- c(0.5, 1.5, 2.5, 3.5, 5)
  ColScheme <- c("darkgreen", "dodgerblue1", "gold2", "white")
  
  terra::image(
    X,
    xlab = NA, ylab = NA,
    xaxt = "n", yaxt = "n",
    col = ColScheme,
    breaks = breakseq,
    asp = 1
  )
  
  terra::plot(outline, add = TRUE, border = "black", col = NA, lwd = 0.4)
  
  legend(
    "topright",
    legend = c("E1 (high)", "E2 (moderate)", "E3 (low)"),
    fill = ColScheme[1:3],
    bty = "n",
    cex = 1.0,
    title = "Historical suitability",
    inset = c(0.02, 0.02)
  )
  
  invisible(NULL)
}

## Plots: --------------
bc_ol <- vect("data-raw/data_tables/bc_outline.gpkg")
plot_historical_suit_map(dbCon = con, spp = "Fd", edatope = "D6", bgc_template = bgc_template, outline = bc_ol)

















#No longer using this plot, but I'm saving the code temporarily in case we decide something similar would useful:
# # RETIRED: Set up for suitability stacked bar chart ------------------
# 
# 
# # Pull in latest suitability.csv. This is saved on the feas_tables branch, so I downloaded it from github on 12/30/2025 and put it in my tree_summaries_code folder.
#     # Question: what time period are these for? Current or historical?
# suit <- read_csv(here("tree_summaries_code", "suitability.csv"))
# str(suit)
# 
# # cleaning
# suit$bgc_all <- gsub("[^A-Z]", "", suit$bgc) #create new column for only BGC zone, lumping subzones
# unique(suit$bgc_all) #need to remove the extra groups, e.g. from US states, to only include BC BGCs.
# bgc_levels <- names(colScheme)
# suit <- suit |> 
#   filter(bgc_all %in% bgc_levels)
# suit$newsuit <- as.factor(suit$newsuit) # factor for barchart
# 
# # suitability colors
# suit_cols <- c("1" = "#0B6E00", "2" = "#2196F3","3" = "#F2D200")
# 
# 
# # stacked bar for Cw
# 
# Cw_bgc_bar <- ggplot(suit[suit$spp == "Cw"& suit$newsuit %in% c(1, 2, 3), ], 
#        aes(x = reorder(bgc_all, newsuit == "1", FUN = mean), fill = factor(newsuit, levels = c("3", "2", "1")))) + #reversed level order here in order to plot E1 on the Left
#   geom_bar(position = "fill") +
#   coord_flip() +
#   scale_y_continuous( #it still thinks this is the y axis, despite coord flip.
#     breaks = c(0, 0.25, 0.5, 0.75, 1),
#     labels = c("0", "25", "50", "75", "100")
#   ) +
#   scale_fill_manual(
#     values = suit_cols,
#     labels = c("E3", "E2", "E1"),
#     name = "Suitability",
#     guide = guide_legend(reverse = TRUE)) + #double reverses the legend to be ascending correctly
#   labs(y="Baseline environmental suitability (%)", x = "BGC zone") +
#   theme_classic(); Cw_bgc_bar
# 
# # save!  
# ggsave("Cw_histbar.png", plot = Cw_bgc_bar, path = here("tree_summaries_code", "Figures"),width = 6, height = 4, dpi = 300, bg = "white")
# 
# 
# 
# 
# 
# # stacked bar for Pl
# Pl_bgc_bar <- ggplot(suit[suit$spp == "Pl"& suit$newsuit %in% c(1, 2, 3), ], 
#                      aes(x = reorder(bgc_all, newsuit == "1", FUN = mean), fill = factor(newsuit, levels = c("3", "2", "1")))) + #reversed level order here in order to plot E1 on the Left
#   geom_bar(position = "fill") +
#   coord_flip() +
#   scale_y_continuous( #it still thinks this is the y axis, despite coord flip.
#     breaks = c(0, 0.25, 0.5, 0.75, 1),
#     labels = c("0", "25", "50", "75", "100")
#   ) +
#   scale_fill_manual(
#     values = suit_cols,
#     labels = c("E3", "E2", "E1"),
#     name = "Suitability",
#     guide = guide_legend(reverse = TRUE)) + #double reverses the legend to be ascending correctly
#   labs(y="Baseline environmental suitability (%)", x = "BGC zone") +
#   theme_classic(); Pl_bgc_bar
# 
# ggsave("Pl_histbar.png", plot = Pl_bgc_bar, path = here("tree_summaries_code", "Figures"),width = 6, height = 4, dpi = 300, bg = "white")


