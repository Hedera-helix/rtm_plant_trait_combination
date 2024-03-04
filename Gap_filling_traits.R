## Gap filling for Plant Trait Data ##

#install.packages('devtools') #assuming it is not already installed
#library(devtools)
#install_github('andreacirilloac/updateR')
#library(updateR)
#updateR()
#library(devtools)
#install_github('fisw10/BHPMF')
#library(usethis) 
#usethis::edit_r_environ()

library(BHPMF)
library(fuzzyjoin)
library(dplyr)
library(ggplot2)
library(hexbin)
library(ggpubr)
library(gridExtra)
library(purrr)
library(RColorBrewer)
library(Hmisc)

# set directory
setwd("~/Data Storage Folder")


## Read the input data ##

# Eya's data
#df_e <- read.csv("eyadata_processed_nospecies.csv", sep = ";")  # Read the matrix X
#hierarchy_e <- read.csv("hierarchy_processed.csv") # Read the hierarchy information

# TRY
#df <- read.csv("TRY_values_processed_finished.csv", sep = ",")  # Read the matrix X
#hierarchy <- read.csv("TRY_hierarchy_processed_finished.csv") # Read the hierarchy information

## TRY and eya combined ##
df <- read.csv("data_gapfilling_LMA_EWT_imputed.csv", sep = ",")  # Read the matrix X
hierarchy <- read.csv("gapfilling_combined_hierarchy.csv") # Read the hierarchy information

# TRY TO COMBINE BOTH -> solve NA problems
#df <- cbind(df, hierarchy)



# Look at the data distribution

hist(t(df[10]), main = colnames(df[10]))



plot(t(df[2]), t(df[5]), ylim = c(0, 200))

boxplot(df[4], xlab = "Nitrogen", ylab = "g/m2", outline = TRUE) 
boxplot(df[2], xlab = "Chlorophyll", ylab = "myg/cm2", outline = TRUE) 

boxplot(df[5], xlab = "LMA", ylab = "g/m2", outline = TRUE) 
boxplot(df[5], xlab = "LMA", ylab = "g/m2", outline = FALSE) 
boxplot(df[6], xlab = "EWT", ylab = "mg/cm2", outline = TRUE) 
boxplot(df[6], xlab = "EWT", ylab = "mg/cm2", outline = FALSE) 

boxplot(df[7], xlab = "Anthocyanin", ylab = "g/cm2", outline = TRUE) 
boxplot(df[8], xlab = "Carotenoid", ylab = "g/cm2", outline = TRUE) 
boxplot(df[9], xlab = "LAI", ylab = "m2/m2", outline = TRUE) 



############ Remove outliers from dataset => improvements for gap-filling ranges

df$LMA_gm2[df$LMA_gm2 >= 200] <- NA
df$LMA_gm2[df$LMA_gm2 <= 0] <- NA

df$EWT_mgcm2[df$EWT_mgcm2 >= 60] <- NA
df$EWT_mgcm2[df$EWT_mgcm2 <= 0] <- NA

df$carbonc_mgg[df$carbonc_mgg >= 600] <- NA
df$carbonc_mgg[df$carbonc_mgg <= 300] <- NA

df$DispersalUL_mm[df$DispersalUL_mm >= 25] <- NA

df$seed_mass_mg[df$seed_mass_mg >= 5000] <- NA

df$plant_height_m[df$plant_height_m >= 75] <- NA
df$plant_height_m[df$plant_height_m <= 0] <- NA

# Remove negative values from previous pigment gap-filling

df$caroc_gcm2[df$caroc_gcm2 <= 0] <- NA
df$antc_gcm2[df$antc_gcm2 <= 0] <- NA




# also remove problematic/unnecessary columns => esp. non-normalized traits like leaf area
df <- df[, c(1, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17:21)] 
#df <- df[, c(1, 4, 5, 8, 9, 11, 15, 17:21)] 
#df <- df[, c(3, 5, 8, 17:21)] 


df_select <- df[apply(is.na(df), 1, all),]

if (nrow(df_select) > 0){
  # get the row numbers (they remain from the original data frame)
  del_rows <- as.integer(rownames(df_select))
  
  # remove the rows from df _and_ the hierarchy
  df <- df[-c(del_rows), ] 
  hierarchy <- hierarchy[-c(del_rows), ]
}



######################################################
# ONLY FOR PROCEDURE TESTING: Validation 

# Plan: Randomly set 10% of variable to NA, then fill gaps and compare with real data

# Temporarily cbind dataframes to make validation samples
df <- cbind(df, hierarchy)

sample_indices_list <- list()
sample_data_list <- list()

# Set a seed for reproducibility
set.seed(123)

# always adjust for loop to subsets! 
for (i in colnames(df[, c(3, 4, 8:12)])) {
  

non_na_indices <- which(!is.na(df[[i]]))
sample_size <- floor(0.1 * length(non_na_indices))
sample_indices <- sample(non_na_indices, size = sample_size, replace = FALSE) # indices for row removal
sample_indices_list[[i]] <- df$X[sample_indices]
sample_data <- df[sample_indices, i]
sample_data_list[[i]] <- sample_data

# Set them to NA in the original dataframe
df[sample_indices, i] <- NA
non_na_indices2 <- which(!is.na(df[[i]])) # to test if it worked

}


# Separate dataframes again to remove NAs
hierarchy <- df[(ncol(df)-3):ncol(df)]
df <- df[1:(ncol(df)-4)]
#########################################################################




# make new data frame with all rows that now have only NA
df_select <- df[apply(is.na(df), 1, all),]

if (nrow(df_select) > 0){
  # get the row numbers (they remain from the original data frame)
  del_rows <- as.integer(rownames(df_select))
  
  # remove the rows from df3 _and_ the hierarchy
  df <- df[-c(del_rows), ] 
  hierarchy <- hierarchy[-c(del_rows), ]
}



# Remove columns with all NA
not_all_na <- function(x) any(!is.na(x))

df1 <- df %>% select(where(not_all_na))

# quickly cbind to remove rows
df1 <- cbind(df1, hierarchy)

# settings should ensure correct partitioning
df2 <- df1 %>% filter(!if_all(.cols = 1:(ncol(df)-4), is.na))

# and separate them again
hierarchy <- df2[(ncol(df2)-3):ncol(df2)]
df2 <- df2[1:(ncol(df2)-4)]

#df2 <- df2[row.names(df2) %in% row.names(hierarchy2), ]


# log transform the data
df3 <- log(df2)

###################### Pre-process error fixing
# => Remove rows with inf values produced by log() 
# make new data frame with all rows that have only NA
df_find <- df3[apply(is.na(df3), 1, all),]

if (nrow(df_find) > 0){
  # get the row numbers (they remain from the original data frame)
  del_rows <- as.integer(rownames(df_find))
  
  # remove the rows from df3 _and_ the hierarchy
  df3 <- df3[-c(del_rows), ] 
  hierarchy <- hierarchy[-c(del_rows), ]
}


# now scale the data
df4 <- scale(df3)


# option to test usefulness of log and scale
#df4 <- df2


###################
# Gap Filling
###################
GapFilling(as.matrix(df4), as.matrix(hierarchy), tuning = FALSE,
           verbose = TRUE, 
           mean.gap.filled.output.path = "mean_gap_filled_combined.txt",
           std.gap.filled.output.path="std_gap_filled_combined.txt")


#############################
## Load the finished data: ##

gap_filled <- read.csv("mean_gap_filled_combined.txt", sep = "")

meanv <- attributes(df4)$`scaled:center`

stdv <- attributes(df4)$`scaled:scale`

# Back-transform all columns of gap_filled in a loop

df_bef = NULL

for (i in 1:ncol(gap_filled)) {
  
  col <- gap_filled[, i]
  result <-  col * stdv[[i]] + meanv[[i]]
  df_bef <- cbind(df_bef, result)
}

df_new <- as.data.frame(df_bef) 
df_final <- exp(df_new)

# Rename all of the columns at once using `setNames()`
new_col_names <- colnames(df)
df_final <- setNames(df_final, new_col_names)    
#df_final <- setNames(gap_filled, new_col_names)         

# View the updated dataframe
head(df_final)


### SAVE DATA AND HIERARCHY

write.csv(df_final, "gap_filled_final.csv", row.names = FALSE)

write.csv(hierarchy, "hierarchy_final.csv", row.names = FALSE)

#
#
#
#
#
#
#
#
#
#
#

# Test data distribution of gap-filled dataframe

# histogram of one variable
hist(t(df_final[4]), main = colnames(df_final[17]))

# scatterplot of two variables
# Chlorophyll - Carotenoids
plot(t(df_final[3]), t(df_final[11]), xlim = c(0, 100), ylim = c(0,20),
     xlab = colnames(df_final[3]), ylab = colnames(df_final[11]))
# LMA - EWT
plot(t(df_final[8]), t(df_final[9]), xlim = c(0, 80), ylim = c(0,80),
     xlab = colnames(df_final[8]), ylab = colnames(df_final[9]))
# Nitrogen - LMA
plot(t(df_final[4]), t(df_final[8]), xlim = c(0, 80), ylim = c(0,80),
     xlab = colnames(df_final[4]), ylab = colnames(df_final[8]))


#AFTER gap filling ####################################################


# select all matching entries of the gap-filled "chlc_mygcm2" column based on the sample_indices
# => finished data frame necessary!(see below for code)

#### Make a loop for indices and data!

# Cbind dataframes again to link original and gap-filled values over X-column
df_bind <- cbind(df_final, hierarchy)

z <- 1

combined_val_list = list()

for (i in names(sample_indices_list)) {
  
  sample_indices_list[[z]]
  sample_data_list[[z]]
  
  # dataframe of new validation values
  new_val_df <- df_bind[df_bind$X %in% sample_indices_list[[z]], c(i, "X")]
  
  # dataframe of old validation values
  old_val_df <- data.frame(cbind(sample_indices_list[[z]], sample_data_list[[z]]))
  
  # merge both dataframes by using the X column of the hierarchy
  combined_val_df <- merge(new_val_df, old_val_df, by.x = "X", by.y = "X1")
  
  # Delete all rows where the new values are greater than the highest old value
  combined_val_list[[i]] <- combined_val_df[combined_val_df[,2] <= max(combined_val_df[,3]), 2:3]  
  
  
  z <- z + 1 
}


# the combined dataframe list
#combined_val_list

df_list_new <- map(combined_val_list, function(df) {
  names(df) <- c("y", "x")
  df
})

#df_list_new <- df_list_new[c(1,2,4,5,6,7)]

plot_list <- list()
q = 1

for (i in df_list_new) {
  
  # Fit a linear model and calculate RMSE and NRMSE
  model <- lm(y ~ x, data = i)
  rmse <- sqrt(mean(model$residuals^2))
  nrmse <- rmse / diff(range(i[1]))

  plot_list[[q]] <- ggplot(i, aes(x = x, y = y)) +
    geom_point(alpha = 0.2) +
    labs(title = names(df_list_new[q]), x = "original", y = "gapfilled") +
    scale_y_continuous(limits = c(0, round(max(df_list_new[[q]][2]))) ) +
    scale_x_continuous(limits = c(0, round(max(df_list_new[[q]][2])))) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    geom_text(x = min(i[2]) + 0.1 * diff(range(i[2])), y = max(i[1]) - 0.1 * diff(range(i[1])),
              label = paste("R2 =", round(summary(model)$r.squared, 2), 
                            "\nNRMSE =", round(nrmse, 2)), size = 2) +
    coord_fixed()  # Set aspect ratio to be 1:1
  
  
  q = q + 1
  
}

# geom_text(x = min(i[2]) + 0.1 * diff(range(i[2])), y = max(i[2])

ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)


# Let's look at hexbin plots for the larger samples!

# the combined dataframe list
#df_list_new
p <- 1
hexbin_list <- list()

for (u in 1: (length(df_list_new))) {
  
 # hexbin_list[[p]] <- hexbinplot(t(df_list_new[[u]][2]) ~ t(df_list_new[[u]][1]), xlab = paste0(names(df_list_new[u]), "_original"), ylab = paste0(names(df_list_new[u]), "_gapfilled"), xbnds = c (0, max(df_list_new[[u]][2])), ybnds = c(0, max(df_list_new[[u]][2])))
  
  #   define breaks and labels for the legend
  #
  brks <- c(0, 1, 10, 20, 50, 100, 200, 300, 400, 500, 600, Inf)
  n.br <- length(brks)
  labs <- c(paste('<=', brks[2:(n.br-1)]), paste('>=', brks[n.br-1]))
  
  
  hexbin_list[[p]] <- ggplot(df_list_new[[u]], aes(x = x, y = y)) +
    geom_hex(aes(fill=cut(..count.., breaks=brks)), color='grey80') +
    scale_fill_manual(name='Count', values = rev(brewer.pal(11, 'Spectral')), labels = labs) +
    scale_y_continuous(limits = c(0, round(max(df_list_new[[u]][2])))) +
    scale_x_continuous(limits = c(0, round(max(df_list_new[[u]][2])))) +
    labs(c(paste0(names(df_list_new[u]), "_original"), paste0(names(df_list_new[u]), "_gapfilled"))) 
  
p = p + 1
}

# Arrange all hexbin plots in a 2x2 grid
ggarrange(plotlist = hexbin_list, 
          #labels = c("A", "B", "C"),
          ncol = 3, nrow = 3)



#######################################################################
# Test area



#######################################################################


gap_filled_final <- read.csv("gap_filled_final.csv", sep = ",")

gap_filled_final <- df_final

# Re-load hierarchy

hierarchy <- read.csv("hierarchy_final.csv", sep = ",")

# Merge gap-filled data with hierarchy

combined_df <- bind_cols(hierarchy, gap_filled_final)

write.csv(combined_df, "gap_filled_final_with_names.csv", row.names = FALSE)


##### Post-processing: remove unrealistic values from gap-filled data (with rounding buffer) #####
gap_filled_final <- gap_filled_final[, c(3, 12:17)]

# Merge gap-filled data with hierarchy

combined_df <- bind_cols(hierarchy, gap_filled_final)

# Rounding principles: to first digit, except when 0, then to the next digit / based on values in Eya's data

combined_df3 <- subset(combined_df, LMA_gm2 >= 5 & LMA_gm2 <= 400 & EWT_mgcm2 >= 0.2 & EWT_mgcm2 <= 81 & nitrogenc_gm2 >= 0.06 & nitrogenc_gm2 <= 10 & antc_gcm2 >= 0.5 & antc_gcm2 <= 3 & caroc_gcm2 >= 1 & caroc_gcm2 <= 22 & LAI_m2m2 >= 0.05 & LAI_m2m2 <= 8 & chlc_mygcm2 >= 4 & chlc_mygcm2 <= 106)
# changed LAI max to 8 for closer simulation of grasslands

#combined_df2 <- subset(combined_df, LMA_gm2 >= 10 & LMA_gm2 <= 200 & EWT_mgcm2 >= 5 & EWT_mgcm2 <= 75 & nitrogenc_gm2 >= 0.5 & nitrogenc_gm2 <= 20 & antc_gcm2 >= 0.1 & antc_gcm2 <= 5 & caroc_gcm2 >= 1 & caroc_gcm2 <= 25 & LAI_m2m2 >= 0.1 & LAI_m2m2 <= 8 & chlc_mygcm2 >= 5 & chlc_mygcm2 <= 120)

write.csv(combined_df3, "gap_filled_final_7targets_with_ranges.csv", row.names = FALSE)



#######################
# load finished dataframe to check values

#gap_filled_final_names <- combined_df

gap_filled_final_test <- read.csv("gap_filled_final_7targets_with_ranges.csv", sep = ",")










