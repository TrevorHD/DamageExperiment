##### Load libraries and data -----------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(grid)
library(gridBase)
library(lme4)
library(lmerTest)
library(survival)
library(emmeans)
library(vcd)
library(coxme)
library(xlsx)

# Load data from local copy of datasheet
Data_GE <- read.xlsx("Data/ThistleData.xlsx", sheetName = "General")
Data_TR <- read.xlsx("Data/ThistleData.xlsx", sheetName = "Trimming")





##### Combine general and trimming data -------------------------------------------------------------------

# Drop redundant rows  where plant is dead (and thus no data exist)
Data_TR %>% drop_na(Height) -> Data_TR

# Merge trimming data and general plot data
Data <- merge(Data_GE, Data_TR, by = c("Row", "Group", "Plant", "Species"))

# Set Boolean warming indicator based on OTC presence
Data %>% mutate(Warmed = case_when(is.na(OTC.On) == FALSE ~ 1, TRUE ~ 0)) -> Data

# Add unique plot/plant identifiers for use in transforms, random effects, etc.
Data$PlotID <- cumsum(!duplicated(Data[, 1:2]))
Data$PlantID <- cumsum(!duplicated(Data[, 1:3]))

# Re-order columns
Data <- Data[, c(1:3, 22:23, 4, 13, 21, 6, 12, 15:19)]





##### Transform data for survival analyses ----------------------------------------------------------------

# Construct placeholder data frame
# ToD is time of death; Cens is censor status (1 = dead, 0 = censored and survived until end)
Data_Alt <- data.frame(matrix(ncol = 11, nrow = 0))
names(Data_Alt) <- c("Row", "Group", "Plant", "PlotID", "PlantID", "Species", "Treatment",
                     "Warmed", "DM_t", "ToD", "Cens")

# Populate data frame
for(i in 1:nrow(Data)){
  if(i < nrow(Data)){
    if(Data$Week[i] >= Data$Week[i + 1]){
      Data_Alt[i, 1:10] <- Data[i, 1:10]}}
  if(i == nrow(Data)){
    Data_Alt[i, 1:10] <- Data[i, 1:10]}
  if(Data$Week[i] == 65){Data_Alt[i, 11] <- 0} else {Data_Alt[i, 11] <- 1}}
Data_Alt <- drop_na(Data_Alt, Row)

# Fix data issue with single observation where time of death is listed as week zero
# This is likely a data entry error and should be week one
Data_Alt$ToD[Data_Alt$ToD == 0] <- 1

# Split the survival data into 2 years (i.e. before and after winter gap at 30 weeks)
# Within each year, survivors at end of year are censored
# This data view will likely not be used in analyses
Data_Alt_Y1 <- Data_Alt
Data_Alt_Y1[Data_Alt_Y1$ToD >= 30, ]$ToD <- 30
Data_Alt_Y1[Data_Alt_Y1$ToD == 30, ]$Cens <- 0
Data_Alt_Y2 <- Data_Alt
Data_Alt_Y2 <- Data_Alt_Y2[Data_Alt_Y2$ToD > 30, ]





##### Transform data for height/stem models ---------------------------------------------------------------

# Calculate regrowth since previous trim
# For TRT 1 (control, no trim), simply use height delta since last measurement
# For TRTs 2 or 3, same as above if previous was uncut, or subtract cut height if previous was cut
# For TRT 4, cut height is zero, so height is equal to growth
Data$HGain <- rep(0, nrow(Data))
for(i in 1:length(1:nrow(Data))){
  if(Data$Week[i] != 0){
    if(Data$Treatment[i] == 1){
      Data$HGain[i] <- Data$Height[i] - Data$Height[i - 1]}
    if(Data$Treatment[i] == 2){
      if(Data$Height[i - 1] <= 10){
        Data$HGain[i] <- Data$Height[i] - Data$Height[i - 1]}
      if(Data$Height[i - 1] > 10){
        Data$HGain[i] <- Data$Height[i] - 10}}
    if(Data$Treatment[i] == 3){
      if(Data$Height[i - 1] <= 5){
        Data$HGain[i] <- Data$Height[i] - Data$Height[i - 1]}
      if(Data$Height[i - 1] > 5){
        Data$HGain[i] <- Data$Height[i] - 5}}
    if(Data$Treatment[i] == 4){
      Data$HGain[i] <- Data$Height[i]}}}

# Calculate stem gain since previous trim
# Subtract new total stems from previous baseline (previous total minus stems trimmed)
# Accounts for fact that main stem cannot count as extra stem
# For untrimmed plants, just use current minus previous; will be zero unless new stem pops up
Data$SGain <- rep(0, nrow(Data))
for(i in 1:length(1:nrow(Data))){
  if(Data$Week[i] != 0){
    if(Data$NStems[i - 1] == Data$NStemsT[i - 1]){
      ns_prev <- 1} else {ns_prev <- Data$NStems[i - 1] - Data$NStemsT[i - 1]}
    Data$SGain[i] <- Data$NStems[i] - ns_prev
    if(Data$Treatment[i] == 1){
      Data$SGain[i] <- Data$NStems[i] - Data$NStems[i - 1]}}}

# Create dataframe of summary measures for each plant
# Then add Boolean indicators for whether or not plant budded and flowered
Data %>% 
  group_by(PlantID) %>% 
  summarise(AvgHGain = mean(HGain),    # Average weekly height change (cm)
            AvgSGain = mean(SGain),    # Average weekly stem count change
            MaxHGain = max(HGain),     # Maximum weekly height change (cm)
            MaxSGain = max(SGain),     # Maximum weekly stem count change
            MaxH = max(Height),        # Maximum height observed (cm)
            MaxS = max(NStems),        # Maximum stem count observed
            MaxHeads = max(Heads),     # Maximum flower count observed
            MaxBuds = max(Buds)) -> Data_Sum
Data_Sum <- merge(unique(Data[, 1:9]), Data_Sum, by = "PlantID")
Data_Sum$Flowered <- ifelse(Data_Sum$MaxHeads > 0, 1, 0)
Data_Sum$Budded <- ifelse(Data_Sum$MaxBuds > 0, 1, 0)

# Remove unused/temp variables
remove(i, ns_prev)





##### Set up plotting functions ---------------------------------------------------------------------------

# Function to package survival data in a form that can be easily plotted
km.curve <- function(comp, subs, species, data){
  if(comp == "trim"){
    return(survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                  data = subset(data, Species == species & Warmed == subs)))}
  if(comp == "warm"){
    return(survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                  data = subset(data, Species == species & Treatment == subs)))}}

# Function to plot survival curves comparing trimming treatments
km.plot <- function(km_mod, bottom, left, atext){
  ltypes <- c("solid", "3333", "2323", "1212")
  colours <- c("black", "purple", "green", "orange")
  par(cex.axis = 0.38, cex.lab = 0.5, tcl = -0.15)
  if(bottom == FALSE & left == TRUE){
    par(mar = c(0, 1, 0.15, 0.1))
    plot(km_mod, xlim = c(0, 30), yaxt = "n", xaxt = "n",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    axis(2, at = seq(0, 1, by = 0.2), mgp = c(0.5, 0.05, 0))
    mtext(side = 2, line = 0.53, "Probability of Survival", cex = 0.5)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == TRUE & left == TRUE){
    par(mar = c(1, 1, 0, 0.1))
    plot(km_mod, xlim = c(0, 30), yaxt = "n", xlab = "Weeks", 
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    axis(2, at = seq(0, 1, by = 0.2), mgp = c(0.5, 0.05, 0))
    mtext(side = 2, line = 0.53, "Probability of Survival", cex = 0.5)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == FALSE & left == FALSE){
    par(mar = c(0, 0.1, 0.15, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n",
         lwd = 1.1, lty = ltypes, col = colours)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == TRUE & left == FALSE){
    par(mar = c(1, 0.1, 0, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n",
         lwd = 1.1, lty = ltypes, col = colours)
    axis(1, at = seq(0, 65, by = 5), mgp = c(0.2, -0.3, 0))
    mtext(side = 1, line = 0.20, "Weeks", cex = 0.5)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}}

# Function to plot survival curves comparing warming treatments
km.plot2 <- function(km_mod, row, left, atext){
  ltypes <- c("solid", "3333")
  colours <- c("blue", "red")
  par(cex.axis = 0.38, cex.lab = 0.5, tcl = -0.15)
  if(row == 1 & left == TRUE){
    par(mar = c(0, 1, 0.15, 0.1))
    plot(km_mod, xlim = c(0, 30), yaxt = "n", xaxt = "n",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    axis(2, at = seq(0, 1, by = 0.2), mgp = c(0.5, 0.05, 0))
    mtext(side = 2, line = 0.53, "Probability of Survival", cex = 0.5)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if((row == 2 | row == 3) & left == TRUE){
    par(mar = c(0, 1, 0.15, 0.1))
    plot(km_mod, xlim = c(0, 30), yaxt = "n", xaxt = "n",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    axis(2, at = seq(0, 1, by = 0.2), mgp = c(0.5, 0.05, 0))
    mtext(side = 2, line = 0.53, "Probability of Survival", cex = 0.5)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(row == 4 & left == TRUE){
    par(mar = c(1, 1, 0, 0.1))
    plot(km_mod, xlim = c(0, 30), yaxt = "n", xlab = "Weeks", 
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    axis(2, at = seq(0, 1, by = 0.2), mgp = c(0.5, 0.05, 0))
    mtext(side = 2, line = 0.53, "Probability of Survival", cex = 0.5)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(row == 1 & left == FALSE){
    par(mar = c(0, 0.1, 0.15, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n",
         lwd = 1.1, lty = ltypes, col = colours)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}
  if((row == 2 | row == 3) & left == FALSE){
    par(mar = c(0, 0.1, 0.15, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n",
         lwd = 1.1, lty = ltypes, col = colours)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(row == 4 & left == FALSE){
    par(mar = c(1, 0.1, 0, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n", xlim = c(0, 65),
         lwd = 1.1, lty = ltypes, col = colours)
    axis(1, at = seq(0, 65, by = 5), mgp = c(0.2, -0.3, 0))
    mtext(side = 1, line = 0.20, "Weeks", cex = 0.5)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}}

# Function to (crudely) plot demographic data for prelim analyses
dd.plot <- function(data, yvar, fac, type, xlim){
  if(type == "hist"){
    gplot <- ggplot(data, aes(.data[[yvar]], fill = factor(.data[[fac]]))) +
      geom_histogram(binwidth = 1) + xlim(xlim[1], xlim[2])}
  if(type == "dens"){
    gplot <- ggplot(data, aes(.data[[yvar]], colour = factor(.data[[fac]]))) +
      geom_density(size = 1.2) + xlim(xlim[1], xlim[2])}
  return(gplot)}

