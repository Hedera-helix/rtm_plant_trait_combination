library(data.table)

# set directory
setwd("~/Data Storage Folder")

# TODO: load eyadata and fill up pigment values before gap-filling with correlation coefficients
df <- read.csv("gapfilling_combined_data.csv", sep = ",")

# Rename columns to be filled up for code structure below
colnames(df)[colnames(df) == "chlc_mygcm2"] <- "r"
colnames(df)[colnames(df) == "caroc_gcm2"] <- "m"

###
#SECOND ROUND ONLY
###

df <- as.data.frame(dt)
colnames(df)[colnames(df) == "chlc_mygcm2"] <- "r"
colnames(df)[colnames(df) == "antc_gcm2"] <- "m"

###
###

dt <- data.table(df)
#eyadata <- read.csv("test_new.csv")


# Create a new dataframe with the selected columns
#df_selected <- eyadata[, c("dataset","Chl.content..μg.cm..", "Carotenoid.content..μg.cm..", "Anthocyanin.content..μg.cm..")]
#df_selected <- df[, c("chlc_mygcm2","caroc_gcm2", "antc_gcm2")]

# Remove rows with NAs in chlc
#df_selected <- df_selected[complete.cases(df_selected$chlc_mygcm2), ]


#######
# Method for Imputation
#######

set.seed(123)

rho = cor(dt$m,dt$r,'pairwise')

# calculate linear regression of original data
fit1 = lm(r ~ m, data=dt)
fit2 = lm(m ~ r, data=dt)
# extract the standard errors of regression intercept (in each m & r direction)
# and multiply s.e. by sqrt(n) to get standard deviation 
sd1 = summary(fit1)$coefficients[1,2] * sqrt(dt[!is.na(r), .N])
sd2 = summary(fit2)$coefficients[1,2] * sqrt(dt[!is.na(m), .N])

# find where data points with missing values lie on the regression line
dt[is.na(r), r.imp := coefficients(fit1)[1] + coefficients(fit1)[2] * m] 
dt[is.na(m), m.imp := coefficients(fit2)[1] + coefficients(fit2)[2] * r]

# generate randomised residuals for the missing data, using the s.d. calculated above
dt[is.na(r), r.ran := rnorm(.N, sd=sd1)] 
dt[is.na(m), m.ran := rnorm(.N, sd=sd2)] 

# function that scales the residuals by a factor x, then calculates how close correlation of imputed data is to that of original data
obj = function(x, dt, rho) {
  dt[, r.comp := r][, m.comp := m]
  dt[is.na(r), r.comp := r.imp + r.ran*x] 
  dt[is.na(m), m.comp := m.imp + m.ran*x] 
  rho2 = cor(dt$m.comp, dt$r.comp,'pairwise')
  (rho-rho2)^2
}

# find the value of x that minimizes the discrepancy of imputed versus original correlation
fit = optimize(obj, c(-5,5), dt, rho)

x=fit$minimum
dt[, r.comp := r][, m.comp := m]
dt[is.na(r), r.comp := r.imp + r.ran*x] 
dt[is.na(m), m.comp := m.imp + m.ran*x] 
rho2 = cor(dt$m.comp, dt$r.comp,'pairwise')
(rho-rho2)^2  # check that rho2 is approximately equal to rho


fit.comp = lm(r.comp ~ m.comp, data=dt)
plot(dt$m.comp, dt$r.comp)
points(dt$m, dt$r, col="red")
abline(fit1, col="green")
abline(fit.comp, col="blue")
mtext(paste(" Rho =", round(rho,5)), at=-1)
mtext(paste(" Rho2 =", round(rho2, 5)), at=6)


### Rename .comp column to become the new imputed variable columns

###
# CHOOSE ONE: Rename columns using :=
###

dt[, c("chlc_mygcm2", "caroc_gcm2") := .(r.comp, m.comp)]

dt[, c("chlc_mygcm2", "antc_gcm2") := .(r.comp, m.comp)]

# Remove unused columns of data.table by name
columns_to_remove <- c("r", "m", "r.imp", "m.imp", "r.ran", "m.ran", "r.comp", "m.comp")
dt[, (columns_to_remove) := NULL]



# AFTER SECOND ROUND

# Check correlations of all three pigments

df_fin <- as.data.frame(dt)

df_fin_select <- df_fin[, 19:21]

df_fin_select <- na.omit(df_fin_select)

cor_matrix <- cor(df_fin_select)


# Print the correlation matrix
print(cor_matrix)



# Save the finished dataframe

write.csv(df_fin, "data_gapfilling_pigments_imputed.csv", row.names = FALSE)






