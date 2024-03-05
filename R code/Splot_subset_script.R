###############
# Create regional or seasonal subsets of the splot database

# set directory
setwd("~/Data Storage Folder")

# load splot

load("sPlotOpen.RData")



# Load lubridate package
library(lubridate)

# Subset the data frame based on summer months
#df_subset_summer <- subset(header.oa, month(ymd(Date_of_recording)) >= 5 & month(ymd(Date_of_recording)) <= 9)

# Keep only those entries with abundance scale "cover percentage"
df_subset_scale <- subset(DT2.oa, DT2.oa$Abundance_scale == "CoverPerc")

## Use plot IDs to subset the cover percentage dataframe
#df_subset_id <- subset(df_subset_summer, PlotObservationID %in% df_subset_scale$PlotObservationID)
df_subset_id <- header.oa


# Subset the header data frame based on grassland and forest dominant plots
subset_eu<- subset(df_subset_id, Continent == "Europe")

subset_grass<- subset(df_subset_id, Forest == "FALSE" & Grassland == "TRUE" & Shrubland == "FALSE")

subset_forest<- subset(df_subset_id, Forest == "TRUE" & Grassland == "FALSE" & Shrubland == "FALSE")
subset_shrubland<- subset(df_subset_id, Shrubland == "TRUE")

# included in grass/forest division
#subset_eu_wetland<- subset(df_subset_id, Continent == "Europe" & Wetland == "TRUE")

## Choose a number of samples per European land cover type and combine them

subset_training_grass <- subset_grass[sample(nrow(subset_grass), 7000), ]

subset_training_forest <-  subset_forest[sample(nrow(subset_forest), 9000), ]

subset_training_shrubland <-  subset_shrubland[sample(nrow(subset_shrubland), 7000), ]


## Use plot IDs to subset the cover percentage dataframe (target dataframe)
subset_training_grass_fin <- subset(df_subset_scale, PlotObservationID %in% subset_training_grass$PlotObservationID)
subset_training_forest_fin <- subset(df_subset_scale, PlotObservationID %in% subset_training_forest$PlotObservationID)
subset_training_shrubland_fin <- subset(df_subset_scale, PlotObservationID %in% subset_training_shrubland$PlotObservationID)
#subset_training_forest_shrubland_fin <- rbind(subset_training_forest_fin,subset_training_shrubland_fin)
# Save as csv file
write.csv(subset_training_grass_fin, "splot_subset_grass_balanced.csv", row.names = FALSE)
write.csv(subset_training_forest_fin, "splot_subset_forest_balanced.csv", row.names = FALSE)
write.csv(subset_training_shrubland_fin, "splot_subset_shrubland_balanced.csv", row.names = FALSE)







