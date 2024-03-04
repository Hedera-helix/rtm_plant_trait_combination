# install.packages("git2r")
# install.packages("liquidSVM", repos="http://pnp.mathematik.uni-stuttgart.de/isa/steinwart/software/R")
# devtools::install_gitlab('jbferet/prospect')
# devtools::install_gitlab('jbferet/prosail')


library(devtools)
library(liquidSVM)
library(prospect)
library(prosail)


# set directory
setwd("~/Data Storage Folder")

######## Data Loading ########
################################################
################################################

### SPLOT DATA GRASSLANDS ###



splot <- read.csv("df_grasslands_frequencies.csv")
#OR
#splot <- read.csv("sPlot_mean_averages_per_plot_grasslands.csv")



splot1 <- na.omit(splot) # check for NAs

# Get the right value ranges for LMA, EWT and nitrogen!

splot$lma_avg <- splot$lma_avg / 10000

splot$ewt_avg <- splot$ewt_avg / 1000

splot$nit_avg <- splot$nit_avg / 10000

splot$pro_avg <- splot$nit_avg * 4.43 # conversion factor from prospect-pro paper


# Create CBC column
splot$cbc_avg <- splot$lma_avg - splot$pro_avg


# subset for testing / specific datasets
#splot_reduced <- splot[1:1000, ]
splot_reduced <- splot


# define input variables for PROSPECT. 
CHL <- splot_reduced[, 2]; CAR <- splot_reduced[, 7]; ANT <- splot_reduced[, 6]; 
EWT <- splot_reduced[, 5]; 
N = rnorm(nrow(splot_reduced), mean = 2, sd = 0.25); # conservative sd = 0.25
PROT = splot_reduced[, 9];
CBC <- splot_reduced[, 10]; 
BROWN = 0.1 ; # 0.1 as fixed value, default to unrealistic
alpha = 40

# define input variables for SAIL. 
lai <- splot_reduced[, 8];   # LAI
q <- 0.01;      # Hot spot parameter

### Inverse correlation for leaf angle and lai

# Define the original variable
#old_variable <- lai

# Define the minimum and maximum values of the original variable
#min_old <- min(old_variable)
#max_old <- max(old_variable)

# Define the minimum and maximum values of the new variable
#min_new <- 10
#max_new <- 80

#### REDUCE variation 35 - 55 runiform

# Create the new variable
#leaf_angle_cor <- ((max_new - min_new)/(max_old - min_old)) * (max_old - old_variable) + min_new
leaf_angle_cor <- runif(nrow(splot_reduced), min = 35, max = 55)

#######

TypeLidf <- 2;  LIDFa <- leaf_angle_cor;    #mean was 45 before
LIDFb <- NULL;  # leaf inclination distribution function parameters
tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 8); # keep, but for now narrow angles down for error checking
#tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 5)

tto <- runif(nrow(splot_reduced), min = 0.001, max = 5);      

# psi
min <- 0
max <- 180
psi <- min + (max-min)*rbeta(nrow(splot_reduced), 0.2, 0.2);      # geometry of acquisition

### INCREASE VARIATION for extreme run // Insert correlation with LAI

factor = sample(c(0.5,0.5), 1)
rsoil = SpecSOIL$Dry_Soil * factor + SpecSOIL$Wet_Soil * abs(factor-1)

# Soil reflectance for simulations
SpecInputATM <- SpecATM

final_reflectance = list()

for (i in 1:nrow(splot_reduced)) {
  
  # run PROSAIL(Prospect-Pro) with 4SAIL
  Ref_4SAIL <- PRO4SAIL(Spec_Sensor = SpecPROSPECT,
                        CHL = CHL[1], CAR = CAR[i], ANT = ANT[i], 
                        EWT = EWT[i], 
                        N = N[1], PROT = PROT[1], CBC = CBC[1],
                        TypeLidf = TypeLidf,LIDFa = LIDFa[i],LIDFb = LIDFb,
                        lai = lai[1], BROWN = BROWN, alpha = alpha,
                        q = q,tts = tts[i],
                        tto = tto[i],psi = psi[i],rsoil = rsoil)
  
  
  
  # Compute final reflectance from output and soil reflectance
  final_reflectance[[i]] <- as.numeric(Ref_4SAIL$rdot * SpecInputATM$Diffuse_Light + Ref_4SAIL$rsot * SpecInputATM$Direct_Light) / (SpecInputATM$Direct_Light + SpecInputATM$Diffuse_Light)
  
}  


# test how many entries are without NAs
final_reflectance_clean <- Filter(function(a) !any(is.na(a)), final_reflectance)


# Save spectra as R object (insurance)
saveRDS(final_reflectance_clean, file="spectra_grasslands.RData")
# Load again
spectra <- readRDS("spectra_grasslands.RData")

### LAST STEP: combine splot trait values and spectra ###
# convert list to dataframe
spectra_df <- data.frame(matrix(unlist(spectra), nrow=length(spectra), byrow=TRUE))
# rename columns to correct wavelengths
colnames(spectra_df) <- sprintf("%d", 400:2500)
#write.csv(spectra_df, "spectra_df_lai1_10.csv", row.names = FALSE)


# remove protein and cbc columns from splot dataframe
splot_reduced <- splot_reduced[, 1:8]
# Change column names for splot dataframe to Eya's conventions
colnames(splot_reduced) <- c('plot', 'Chl content (μg/cm²)', 'N content (mg/cm²)', 'LMA (g/m²)','EWT (mg/cm²)','Anthocyanin content (μg/cm²)','Carotenoid content (μg/cm²)','LAI (m²/m²)')

# Set correct units for traits
splot_reduced$`LMA (g/m²)` <- splot_reduced$`LMA (g/m²)` * 10000 # to g/m2

splot_reduced$`EWT (mg/cm²)` <- splot_reduced$`EWT (mg/cm²)` * 1000 # to mg/cm2

splot_reduced$`N content (mg/cm²)` <- splot_reduced$`N content (mg/cm²)` * 1000 # to mg/cm2


splot_rtm_data_grasslands <- cbind(splot_reduced, spectra_df)

# small name change to enable merge with Eya's data
names(splot_rtm_data_grasslands)[names(splot_rtm_data_grasslands) == 'plot'] <- 'PlotObservationID'

















### SPLOT DATA FORESTS ###



splot <- read.csv("df_forests_frequencies.csv")
#OR
#splot <- read.csv("sPlot_mean_averages_per_plot_forests.csv")

splot1 <- na.omit(splot) # check for NAs

# Get the right value ranges for LMA, EWT and nitrogen!

splot$lma_avg <- splot$lma_avg / 10000

splot$ewt_avg <- splot$ewt_avg / 1000

splot$nit_avg <- splot$nit_avg / 10000

splot$pro_avg <- splot$nit_avg * 4.43 # conversion factor from prospect-pro paper


# Create CBC column
splot$cbc_avg <- splot$lma_avg - splot$pro_avg


# subset for testing / specific datasets
#splot_reduced <- splot[1:1000, ]
splot_reduced <- splot


# define input variables for PROSPECT. 
# refer to prospect tutorial for default values corresponding to undefined PROSPECT variables
CHL <- splot_reduced[, 2]; CAR <- splot_reduced[, 7]; ANT <- splot_reduced[, 6]; 
EWT <- splot_reduced[, 5]; 
N = rnorm(nrow(splot_reduced), mean = 2, sd = 0.25); # conservative sd = 0.25
PROT = splot_reduced[, 9];
CBC <- splot_reduced[, 10]; 
BROWN = 0.1 ; # 0.1 as fixed value, default to unrealistic
alpha = 40

# define input variables for SAIL. 
lai <- splot_reduced[, 8];   # LAI
q <- 0.01;      # Hot spot parameter

### Inverse correlation for leaf angle and lai

# Define the original variable
#old_variable <- lai

# Define the minimum and maximum values of the original variable
#min_old <- min(old_variable)
#max_old <- max(old_variable)

# Define the minimum and maximum values of the new variable
#min_new <- 10
#max_new <- 80

# Create the new variable
#leaf_angle_cor <- ((max_new - min_new)/(max_old - min_old)) * (max_old - old_variable) + min_new
leaf_angle_cor <- runif(nrow(splot_reduced), min = 35, max = 55)

#######

TypeLidf <- 2;  LIDFa <- leaf_angle_cor;    #mean was 45 before
LIDFb <- NULL;  # leaf inclination distribution function parameters
tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 8); # keep, but for now narrow angles down for error checking
#tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 5)

tto <- runif(nrow(splot_reduced), min = 0.001, max = 5);     

# psi
min <- 0
max <- 180
psi <- min + (max-min)*rbeta(nrow(splot_reduced), 0.2, 0.2);      # geometry of acquisition

#rsoil <- SpecSOIL$Dry_Soil     # use weighted average instead!     

factor = sample(c(0.5,0.5), 1)
rsoil = SpecSOIL$Dry_Soil * factor + SpecSOIL$Wet_Soil * abs(factor-1)

# Soil reflectance for simulations
SpecInputATM <- SpecATM

# ADDITIONAL PARAMETER FOR 4SAIL2

fraction_brown <- rep(0, nrow(splot_reduced))
#fraction_brown <- runif(nrow(splot_reduced), min = 0, max = 0.1) # the fraction of LAI corresponding to brown leaf area, between 0 and 1
diss <- 0.5 #rnorm(2101, mean = 0.5, sd = 0.1) # layer dissociation factor, measure for how heterogeneous the distribution of green and brown leaves in the top and bottom canopy layer is
Cv <- runif(nrow(splot_reduced), min = 0.5, max = 1) # vertical crown cover percentage (= % ground area covered with crowns as seen from nadir direction)
Zeta <- rnorm(nrow(splot_reduced), mean = 0.5, sd = 0.25) # tree shape factor (= ratio of crown diameter to crown height)

# DEFINE DECAY FACTORS FOR 4SAIL2 
chl_decay <- 0.125
car_decay <- 0.625
ant_decay <- 3
ewt_decay <- 0.5
n_decay <- 1.333
prot_decay <- 0.5
cbc_decay <- 0.888


# DEFINE DECAY FACTORS FOR 4SAIL2 
#chl_decay <- 1
#car_decay <- 1
#ant_decay <- 1
#ewt_decay <- 1
#n_decay <- 1
#prot_decay <- 1
#cbc_decay <- 1


final_reflectance = list()

for (i in 1:nrow(splot_reduced)) {
  
  # run PROSAIL(Prospect-Pro) with 4SAIL
  Ref_4SAIL <- PRO4SAIL(SAILversion = '4SAIL2', Spec_Sensor = SpecPROSPECT,
                        CHL = c(CHL[1],CHL[1]*chl_decay), CAR = c(CAR[i],CAR[i]*car_decay), 
                        ANT = c(ANT[i],ANT[i]*ant_decay), EWT = c(EWT[i],EWT[i]*ewt_decay), 
                        N = c(N[1],N[1]*n_decay), PROT = c(PROT[1],PROT[1]*prot_decay), CBC = c(CBC[1],CBC[1]*cbc_decay),
                        TypeLidf = TypeLidf,LIDFa = LIDFa[i],LIDFb = LIDFb,
                        lai = lai[1], BROWN = BROWN, alpha = alpha,
                        q = q, tts = tts[i],
                        tto = tto[i], psi = psi[i], rsoil = rsoil, 
                        fraction_brown = fraction_brown[i], diss = diss, Cv = Cv[i], Zeta = Zeta[i])
  
  
  
  
  # Compute final reflectance from output and soil reflectance
  final_reflectance[[i]] <- as.numeric(Ref_4SAIL$rdot * SpecInputATM$Diffuse_Light + Ref_4SAIL$rsot * SpecInputATM$Direct_Light) / (SpecInputATM$Direct_Light + SpecInputATM$Diffuse_Light)
  
}  


# test how many entries are without NAs
final_reflectance_clean <- Filter(function(a) !any(is.na(a)), final_reflectance)


# Save spectra as R object (insurance)
saveRDS(final_reflectance_clean, file="spectra_forests.RData")
# Load again
spectra <- readRDS("spectra_forests.RData")

### LAST STEP: combine splot trait values and spectra ###
# convert list to dataframe
spectra_df <- data.frame(matrix(unlist(spectra), nrow=length(spectra), byrow=TRUE))
# rename columns to correct wavelengths
colnames(spectra_df) <- sprintf("%d", 400:2500)
#write.csv(spectra_df, "spectra_df_lai1_10.csv", row.names = FALSE)


# remove protein and cbc columns from splot dataframe
splot_reduced <- splot_reduced[, 1:8]
# Change column names for splot dataframe to Eya's conventions
colnames(splot_reduced) <- c('plot', 'Chl content (μg/cm²)', 'N content (mg/cm²)', 'LMA (g/m²)','EWT (mg/cm²)','Anthocyanin content (μg/cm²)','Carotenoid content (μg/cm²)','LAI (m²/m²)')

# Set correct units for traits
splot_reduced$`LMA (g/m²)` <- splot_reduced$`LMA (g/m²)` * 10000 # to g/m2

splot_reduced$`EWT (mg/cm²)` <- splot_reduced$`EWT (mg/cm²)` * 1000 # to mg/cm2

splot_reduced$`N content (mg/cm²)` <- splot_reduced$`N content (mg/cm²)` * 1000 # to mg/cm2


splot_rtm_data_forests <- cbind(splot_reduced, spectra_df)

# small name change to enable merge with Eya's data
names(splot_rtm_data_forests)[names(splot_rtm_data_forests) == 'plot'] <- 'PlotObservationID'






### SPLOT DATA SHRUBLAND ###



splot <- read.csv("df_shrublands_frequencies.csv")
#OR
#splot <- read.csv("sPlot_mean_averages_per_plot_shrublands.csv")


splot1 <- na.omit(splot) # check for NAs

# Get the right value ranges for LMA, EWT and nitrogen!

splot$lma_avg <- splot$lma_avg / 10000

splot$ewt_avg <- splot$ewt_avg / 1000

splot$nit_avg <- splot$nit_avg / 10000

splot$pro_avg <- splot$nit_avg * 4.43 # conversion factor from prospect-pro paper


# Create CBC column
splot$cbc_avg <- splot$lma_avg - splot$pro_avg


# subset for testing / specific datasets
#splot_reduced <- splot[1:1000, ]
splot_reduced <- splot


# define input variables for PROSPECT. 
# refer to prospect tutorial for default values corresponding to undefined PROSPECT variables
CHL <- splot_reduced[, 2]; CAR <- splot_reduced[, 7]; ANT <- splot_reduced[, 6]; 
EWT <- splot_reduced[, 5]; 
N = rnorm(nrow(splot_reduced), mean = 2, sd = 0.25); # conservative sd = 0.25
PROT = splot_reduced[, 9];
CBC <- splot_reduced[, 10]; 
BROWN = 0.1 ; # 0.1 as fixed value, default to unrealistic
alpha = 40

# define input variables for SAIL. 
lai <- splot_reduced[, 8];   # LAI
q <- 0.01;      # Hot spot parameter

### Inverse correlation for leaf angle and lai

# Define the original variable
#old_variable <- lai

# Define the minimum and maximum values of the original variable
#min_old <- min(old_variable)
#max_old <- max(old_variable)

# Define the minimum and maximum values of the new variable
#min_new <- 10
#max_new <- 80

# Create the new variable
#leaf_angle_cor <- ((max_new - min_new)/(max_old - min_old)) * (max_old - old_variable) + min_new
leaf_angle_cor <- runif(nrow(splot_reduced), min = 35, max = 55)

#######

TypeLidf <- 2;  LIDFa <- leaf_angle_cor;    #mean was 45 before
LIDFb <- NULL;  # leaf inclination distribution function parameters
tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 8); # keep, but for now narrow angles down for error checking
#tts <- rnorm(nrow(splot_reduced), mean = 36, sd = 5)

tto <- runif(nrow(splot_reduced), min = 0.001, max = 5);      

# psi
min <- 0
max <- 180
psi <- min + (max-min)*rbeta(nrow(splot_reduced), 0.2, 0.2);      # geometry of acquisition


factor = sample(c(0.5,0.5), 1)
rsoil = SpecSOIL$Dry_Soil * factor + SpecSOIL$Wet_Soil * abs(factor-1)

# Soil reflectance for simulations
SpecInputATM <- SpecATM

# ADDITIONAL PARAMETER FOR 4SAIL2

fraction_brown <- rep(0, nrow(splot_reduced))
#fraction_brown <- runif(nrow(splot_reduced), min = 0, max = 0.1) # the fraction of LAI corresponding to brown leaf area, between 0 and 1
diss <- 0.5 #rnorm(2101, mean = 0.5, sd = 0.1) # layer dissociation factor, measure for how heterogeneous the distribution of green and brown leaves in the top and bottom canopy layer is
Cv <- runif(nrow(splot_reduced), min = 0.1, max = 0.6) # vertical crown cover percentage (= % ground area covered with crowns as seen from nadir direction)
Zeta <- rnorm(nrow(splot_reduced), mean = 2, sd = 0.5) # tree shape factor (= ratio of crown diameter to crown height)

# DEFINE DECAY FACTORS FOR 4SAIL2 
chl_decay <- 0.125
car_decay <- 0.625
ant_decay <- 3
ewt_decay <- 0.5
n_decay <- 1.333
prot_decay <- 0.5
cbc_decay <- 0.888

# DEFINE DECAY FACTORS FOR 4SAIL2 
#chl_decay <- 1
#car_decay <- 1
#ant_decay <- 1
#ewt_decay <- 1
#n_decay <- 1
#prot_decay <- 1
#cbc_decay <- 1


final_reflectance = list()

for (i in 1:nrow(splot_reduced)) {
  
  # run PROSAIL(Prospect-Pro) with 4SAIL
  Ref_4SAIL <- PRO4SAIL(SAILversion = '4SAIL2', Spec_Sensor = SpecPROSPECT,
                        CHL = c(CHL[1],CHL[1]*chl_decay), CAR = c(CAR[i],CAR[i]*car_decay), 
                        ANT = c(ANT[i],ANT[i]*ant_decay), EWT = c(EWT[i],EWT[i]*ewt_decay), 
                        N = c(N[1],N[1]*n_decay), PROT = c(PROT[1],PROT[1]*prot_decay), CBC = c(CBC[1],CBC[1]*cbc_decay),
                        TypeLidf = TypeLidf,LIDFa = LIDFa[i],LIDFb = LIDFb,
                        lai = lai[1], BROWN = BROWN, alpha = alpha,
                        q = q, tts = tts[i],
                        tto = tto[i], psi = psi[i], rsoil = rsoil, 
                        fraction_brown = fraction_brown[i], diss = diss, Cv = Cv[i], Zeta = Zeta[i])
  
  
  
  
  # Compute final reflectance from output and soil reflectance
  final_reflectance[[i]] <- as.numeric(Ref_4SAIL$rdot * SpecInputATM$Diffuse_Light + Ref_4SAIL$rsot * SpecInputATM$Direct_Light) / (SpecInputATM$Direct_Light + SpecInputATM$Diffuse_Light)
  
}  


# test how many entries are without NAs
final_reflectance_clean <- Filter(function(a) !any(is.na(a)), final_reflectance)


# Save spectra as R object (insurance)
saveRDS(final_reflectance_clean, file="spectra_shrublands.RData")
# Load again
spectra <- readRDS("spectra_shrublands.RData")

### LAST STEP: combine splot trait values and spectra ###
# convert list to dataframe
spectra_df <- data.frame(matrix(unlist(spectra), nrow=length(spectra), byrow=TRUE))
# rename columns to correct wavelengths
colnames(spectra_df) <- sprintf("%d", 400:2500)
#write.csv(spectra_df, "spectra_df_lai1_10.csv", row.names = FALSE)


# remove protein and cbc columns from splot dataframe
splot_reduced <- splot_reduced[, 1:8]
# Change column names for splot dataframe to Eya's conventions
colnames(splot_reduced) <- c('plot', 'Chl content (μg/cm²)', 'N content (mg/cm²)', 'LMA (g/m²)','EWT (mg/cm²)','Anthocyanin content (μg/cm²)','Carotenoid content (μg/cm²)','LAI (m²/m²)')

# Set correct units for traits
splot_reduced$`LMA (g/m²)` <- splot_reduced$`LMA (g/m²)` * 10000 # to g/m2

splot_reduced$`EWT (mg/cm²)` <- splot_reduced$`EWT (mg/cm²)` * 1000 # to mg/cm2

splot_reduced$`N content (mg/cm²)` <- splot_reduced$`N content (mg/cm²)` * 1000 # to mg/cm2


splot_rtm_data_shrublands <- cbind(splot_reduced, spectra_df)

# small name change to enable merge with Eya's data
names(splot_rtm_data_shrublands)[names(splot_rtm_data_shrublands) == 'plot'] <- 'PlotObservationID'






##############################
##############################


# Combine the finished rtm datasets
splot_rtm_data_combined <- rbind(splot_rtm_data_grasslands, splot_rtm_data_forests, splot_rtm_data_shrublands)


# Save as csv 
#write.csv(splot_rtm_data_combined, "splot_rtm_data_grass_forest_shrub_balanced.csv", row.names = FALSE)
#OR
write.csv(splot_rtm_data_combined, "splot_rtm_data_grass_forest_shrub_balanced_frequencies.csv", row.names = FALSE)


# Load for testing purposes
#splot_rtm_data_combined <- read.csv("splot_rtm_data_grass_forest_shrub_balanced.csv")



