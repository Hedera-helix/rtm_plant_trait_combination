
library(dplyr)
library(purrr)
library(sampling)
library(scales)



# set directory
setwd("~/Data Storage Folder")

df_eya <- read.csv("test_new_corrected.csv")

#df_eya_subset <- data.frame(df_eya[, c("dataset", "Chl.content..μg.cm..", "N.content..mg.cm..", "LMA..g.m..", "EWT..mg.cm..", "Anthocyanin.content..μg.cm..", "Carotenoid.content..μg.cm..", "LAI..m..m..")])
df_eya_subset <- data.frame(df_eya[, c("Chl.content..μg.cm..", "N.content..mg.cm..", "LMA..g.m..", "EWT..mg.cm..", "Anthocyanin.content..μg.cm..", "Carotenoid.content..μg.cm..", "LAI..m..m..")])


df_eya_subset$N.content..mg.cm.. <- df_eya_subset$N.content..mg.cm.. * 10
head(df_eya_subset$N.content..mg.cm..)

# Loading: CHOOSE ONE
df_new <- read.csv("sPlot_mean_averages_per_plot_grasslands.csv")
df_new <- read.csv("sPlot_mean_averages_per_plot_forests.csv")
df_new <- read.csv("sPlot_mean_averages_per_plot_shrublands.csv")


# Function to identify rows with outliers based on quantiles
remove_outlier_rows <- function(df) {
  
  lower_quantile <- apply(df, 2, function(x) quantile(x, 0.05, na.rm = FALSE)) # ADJUST VALUES IF THERE IS ERROR WITH ROW LENGTHS
  upper_quantile <- apply(df, 2, function(x) quantile(x, 0.95, na.rm = FALSE))
  rows_to_keep <- apply(df, 1, function(row) all(row >= lower_quantile & row <= upper_quantile))
  return(df[rows_to_keep, ])
}

# Function for NA -> 0
set_na_to_zero <- function(x) {
  x[is.na(x)] <- 0
  return(x)
}

# Custom break function to always have the specified number of bins
custom_breaks <- function(data, num_bins) {
  range_data <- range(data)
  bin_width <- (range_data[2] - range_data[1]) / num_bins
  breaks <- seq(from = range_data[1], to = range_data[2], by = bin_width)
  breaks <- c(breaks, breaks[length(breaks)] + bin_width)  # Add an extra break to cover the last bin
  return(breaks)
}


# Apply the function to the trait columns
#trait_columns <- c("chl_avg", "nit_avg", "lma_avg", "ewt_avg", "ant_avg", "car_avg", "lai_avg")
#df_inter <- remove_outlier_rows(df_new[trait_columns])



#df_new_processed <- merge(df_inter, df_new, by = trait_columns)
# Reorder columns
#df_new_active <- df_new_processed[, c(ncol(df_new_processed), 1:(ncol(df_new_processed) - 1))]

#df_new_active <- df_new

# Create an empty data frame to store all the sampled values
all_dataframe_list <- list()
x2_list <- list()

# Define the number of bins you want (n-1)
num_bins <- 99

# Define total number of samples
sample_total <- round(nrow(df_eya)/7)


## TEST ##

df_new_active <- df_new


# Loop through the selected column
for (z in 2:length(df_new_active)) {
  #z = 2
  
  #df_new_act_subset <- cbind(df_new_active[1], df_new_active[z])
  
  df_new_act_subset <- df_new_active
  
  my_vector <- df_eya_subset[z - 1]
  
  # Calculate the 0-10% and 90-100% quantiles
  quantiles_0_to_10 <- quantile(my_vector, c(0, 0.1), na.rm = TRUE)
  quantiles_90_to_100 <- quantile(my_vector, c(0.9, 1), na.rm = TRUE)
  
  # Filter the vector based on the quantiles
  x <- my_vector[my_vector > quantiles_0_to_10[2] & my_vector < quantiles_90_to_100[1]]
  
  x <- na.omit(x)
  
  
  # RESTRCT rtm values to the range of filtered original data (x)
  min_val <- min(x)
  max_val <- max(x)
  
  # Subset the dataframe to keep only rows where col2 is within the range [min_val, max_val]
  df_new_act_subset <- df_new_act_subset[df_new_act_subset[[z]] >= min_val & df_new_act_subset[[z]] <= max_val, ]
  
  # Calculate custom breaks
  breaks <- custom_breaks(x, num_bins)
  
  # Create a histogram to calculate the frequency distribution
  histogram <- hist(x, breaks = breaks, plot = FALSE)
  
  # Print the frequency distribution table
  frequency_distribution <- data.frame(interval = histogram$mids, frequency = histogram$counts)
  x2 <- 1 - frequency_distribution[, 2] / max(frequency_distribution[, 2])
  
  ###  FIX 0 and 1 values for x2!
  
  # Get the interval index for all values in df_new[[z]]
  interval_index <- findInterval(df_new_act_subset[[z]], histogram$breaks)

  # Calculate the number of values in each interval
  interval_counts <- table(interval_index)
  
  # Calculate the total number of values
  total_values <- sum(interval_counts)
  
  
  #x2 <- x2^2
  
  # Calculate the standardized probabilities
  standardized_probabilities <- x2 / sum(x2)
  
  # Calculate the number of samples to draw from each interval
  samples_per_interval <- round(standardized_probabilities * sample_total)
  
  # Reset samples list
  all_samples_list <- list()
  
  
  # Sample from each interval and combine the results
  for (i in 1:length(interval_counts)) {
    interval <- unique(interval_index)[i]
    samples_to_draw <- samples_per_interval[i]
    
    # Sample from the current interval with replacement
    if (!is.na(samples_to_draw)) { 
    sampled_rows <- slice_sample(df_new_act_subset[interval_index == interval, ], n = samples_to_draw, replace = TRUE)
    
    # Add samples to the all_samples_list
    all_samples_list[[i]] <- sampled_rows
    }
  }
  
  
  sampled_dataframe <- bind_rows(all_samples_list, .id = NULL)
  
  all_dataframe_list[[z-1]] <- sampled_dataframe
  
  x2_list[[z-1]] <- x2
}


# Test distribution of individual variables

# Save lower and upper quantile boundaries for df_eya_subset
lower_quantile <- apply(df_eya_subset, 2, function(x) quantile(x, 0.10, na.rm = TRUE)) # ADJUST VALUES IF THERE IS ERROR WITH ROW LENGTHS
upper_quantile <- apply(df_eya_subset, 2, function(x) quantile(x, 0.90, na.rm = TRUE))

par(mfrow=c(2,3)) # To plot histograms in a 2x2 grid

hist(df_eya_subset$Chl.content..μg.cm..[df_eya_subset$Chl.content..μg.cm.. > lower_quantile[[1]] & df_eya_subset$Chl.content..μg.cm.. < upper_quantile[[1]]], breaks = num_bins, main = "real_data_Chlorophyll")
hist(all_dataframe_list[[1]]$chl_avg, breaks = num_bins, main = "rtm_frequency_Chlorophyll")
plot(x2_list[[1]])

hist(df_eya_subset$N.content..mg.cm..[df_eya_subset$N.content..mg.cm.. > lower_quantile[[2]] & df_eya_subset$N.content..mg.cm.. < upper_quantile[[2]]], breaks = num_bins, main = "real_data_Nitrogen")
hist(all_dataframe_list[[2]]$nit_avg, breaks = num_bins, main = "rtm_frequency_Nitrogen")
plot(x2_list[[2]])

hist(df_eya_subset$LMA..g.m..[df_eya_subset$LMA..g.m.. > lower_quantile[[3]] & df_eya_subset$LMA..g.m.. < upper_quantile[[3]]], breaks = num_bins, main = "real_data_LMA")
hist(all_dataframe_list[[3]]$lma_avg, breaks = num_bins, main = "rtm_frequency_LMA")
plot(x2_list[[3]])

hist(df_eya_subset$EWT..mg.cm..[df_eya_subset$EWT..mg.cm.. > lower_quantile[[4]] & df_eya_subset$EWT..mg.cm.. < upper_quantile[[4]]], breaks = num_bins, main = "real_data_EWT")
hist(all_dataframe_list[[4]]$ewt_avg, breaks = num_bins, main = "rtm_frequency_EWT")
plot(x2_list[[4]])

hist(df_eya_subset$Anthocyanin.content..μg.cm..[df_eya_subset$Anthocyanin.content..μg.cm.. > lower_quantile[[5]] & df_eya_subset$Anthocyanin.content..μg.cm.. < upper_quantile[[5]]], breaks = num_bins, main = "real_data_Anthocyanin")
hist(all_dataframe_list[[5]]$ant_avg, breaks = num_bins, main = "rtm_frequency_Anthocyanin")
plot(x2_list[[5]])

hist(df_eya_subset$Carotenoid.content..μg.cm..[df_eya_subset$Carotenoid.content..μg.cm.. > lower_quantile[[6]] & df_eya_subset$Carotenoid.content..μg.cm.. < upper_quantile[[6]]], breaks = num_bins, main = "real_data_Carotenoids")
hist(all_dataframe_list[[6]]$car_avg, breaks = num_bins, main = "rtm_frequency_Carotenoids")
plot(x2_list[[6]])

hist(df_eya_subset$LAI..m..m..[df_eya_subset$LAI..m..m.. > lower_quantile[[7]] & df_eya_subset$LAI..m..m.. < upper_quantile[[7]]], breaks = num_bins, main = "real_data_LAI")
hist(all_dataframe_list[[7]]$lai_avg, breaks = num_bins, main = "rtm_frequency_LAI")
plot(x2_list[[7]])




# combine and later use in rtm, but retain trait affiliation so that values can be set to NA afterwards


# Combine dataframe list to final dataframe
final_dataframe <- bind_rows(all_dataframe_list, .id = NULL)


par(mfrow=c(2,2))
# Print the histograms of the sampled values for each column (except plot numbers)
for (col in names(final_dataframe)[-1]) {
  hist(df_new_active[[col]], breaks = num_bins, xlim = c(min(final_dataframe[[col]]), max(final_dataframe[[col]])), main = paste0(col, " rtm_original"))
  hist(final_dataframe[[col]], breaks = num_bins, xlim = c(min(final_dataframe[[col]]), max(final_dataframe[[col]])), main = paste0(col, " rtm_inverse_frequencies"))
}


#df_combined <- c(df_eya$Chl.content..μg.cm.., final_dataframe$chl_avg)
#hist(df_combined,  breaks = 20)

#dev.off()

# Saving: CHOOSE ONE 

write.csv(final_dataframe, "df_grasslands_frequencies.csv", row.names = FALSE)

# Save plot numbers of dataframe list to delete other variables later
# Save the list as a binary file using saveRDS
saveRDS(all_dataframe_list, file = "trait_frequency_rtm_dataframe_list_grasslands.rds")




write.csv(final_dataframe, "df_forests_frequencies.csv", row.names = FALSE)

# Save plot numbers of dataframe list to delete other variables later
# Save the list as a binary file using saveRDS
saveRDS(all_dataframe_list, file = "trait_frequency_rtm_dataframe_list_forests.rds")




write.csv(final_dataframe, "df_shrublands_frequencies.csv", row.names = FALSE)

# Save plot numbers of dataframe list to delete other variables later
# Save the list as a binary file using saveRDS
saveRDS(all_dataframe_list, file = "trait_frequency_rtm_dataframe_list_shrublands.rds")




# FOR RTM CREATION, TURN ALL OTHER VARIABLES TO NA

# Load the list of dataframes back into R
all_dataframe_list_grasslands <- readRDS(file = "trait_frequency_rtm_dataframe_list_grasslands.rds")

# Load the list of dataframes back into R
all_dataframe_list_forests <- readRDS(file = "trait_frequency_rtm_dataframe_list_forests.rds")

# Load the list of dataframes back into R
all_dataframe_list_shrublands <- readRDS(file = "trait_frequency_rtm_dataframe_list_shrublands.rds")



# Define the three lists of dataframes
all_dataframe_lists <- list(all_dataframe_list_grasslands, all_dataframe_list_forests, all_dataframe_list_shrublands)

# Initialize an empty dataframe to store the final result
final_dataframe <- data.frame()

# Create a function to update the dataframes
update_dataframe <- function(df, keep_columns) {
  for (i in 2:ncol(df)) {
    if (i %in% keep_columns) {
      next
    } else {
      df[, i] <- NA
    }
  }
  return(df)
}

# Define the columns to keep for each dataframe
columns_to_keep <- list(c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(1, 6), c(1, 7), c(1, 8))

# Process each dataframe list and combine them into the final dataframe
for (df_list in all_dataframe_lists) {
  # Update each dataframe in the list
  updated_dataframes <- lapply(1:length(df_list), function(i) {
    update_dataframe(df_list[[i]], columns_to_keep[[i]])
  })
  
  # Combine all dataframes into one by row binding (rbind)
  combined_dataframe <- do.call(rbind, updated_dataframes)
  
  # Append the combined dataframe to the final dataframe
  final_dataframe <- rbind(final_dataframe, combined_dataframe)
}

# Now, final_dataframe contains the combined dataframes from all three lists in one dataframe

#=> DO RTM RUN

# LOAD RTM generated data

splot_rtm_data_combined <- read.csv("splot_rtm_data_grass_forest_shrub_balanced_frequencies.csv")


# Replace columns 2-8 of large_df with corresponding columns from small_df
splot_rtm_data_combined[, 2:8] <- final_dataframe[, 2:8]


rtm_subset1 <- splot_rtm_data_combined[, 1:8]
# Change column names for splot dataframe to Eya's conventions
colnames(rtm_subset1) <- c('plot', 'Chl content (μg/cm²)', 'N content (mg/cm²)', 'LMA (g/m²)','EWT (mg/cm²)','Anthocyanin content (μg/cm²)','Carotenoid content (μg/cm²)','LAI (m²/m²)')

rtm_subset2 <- splot_rtm_data_combined[, 9:2109]
# Remove "X" before numeric columns
colnames(rtm_subset2) <- sprintf("%d", 400:2500)

splot_rtm_data_combined_end <- cbind(rtm_subset1, rtm_subset2)

# Save as csv 
write.csv(splot_rtm_data_combined_end, "splot_rtm_data_grass_forest_shrub_frequency_final.csv", row.names = FALSE)


# Load for testing
#rtm_final_test_frequency <- read.csv("splot_rtm_data_grass_forest_shrub_frequency.csv", sep = ",")

#rtm_final_test_normal <- read.csv("splot_rtm_data_grass_forest_shrub_balanced.csv", sep = ",")




