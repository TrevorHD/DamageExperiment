##### Load libraries and data -----------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(grid)
library(gridBase)
library(lme4)
library(survival)
library(sjPlot)
library(sjmisc)
library(xlsx)

# Load data from local copy of datasheet
Data_GE <- read.xlsx("Data/ThistleData.xlsx", sheetName = "General")
Data_TR <- read.xlsx("Data/ThistleData.xlsx", sheetName = "Trimming")





##### Tidy data and calculate derived quantities ----------------------------------------------------------

# Drop redundant rows where plant is recorded as dead
Data_TR %>% drop_na(Height) -> Data_TR

# Merge trimming data and general plot data
Data <- merge(Data_GE, Data_TR, by = c("Row", "Group", "Plant", "Species"))

# Set Boolean warming indicator based on OTC presence
Data %>% mutate(Warmed = case_when(is.na(OTC.On) == FALSE ~ 1, TRUE ~ 0)) -> Data

# Calculate regrowth since previous trim
# For TRT 1 (control, no trim), simply use growth since last measurement
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

# Select columns used for analyses  
Data %>% select(Row, Group, Plant, Species, Warmed, Treatment, DM_t, Week, Height,
                HGain, Buds, Heads, NStems, SGain) -> Data

# Construct concise representation of original dataset
# ToD indicates time of death
# Cens indicates censor status (1 = dead, 0 = censored and survived until end)
Data_Alt <- data.frame(matrix(ncol = 9, nrow = 0))
names(Data_Alt) <- c("Row", "Group", "Plant", "Species", "Warmed", "Treatment", "DM_t", "ToD", "Cens")
for(i in 1:nrow(Data)){
  if(i < nrow(Data)){
    if(Data$Week[i] >= Data$Week[i + 1]){
      Data_Alt[i, 1:8] <- Data[i, 1:8]}}
  if(i == nrow(Data)){
    Data_Alt[i, 1:8] <- Data[i, 1:8]}
  if(Data$Week[i] == 65){Data_Alt[i, 9] <- 0} else {Data_Alt[i, 9] <- 1}}
Data_Alt <- drop_na(Data_Alt, Row)

# Same as above, but split into before and after winter census gap
Data_Alt_1 <- data.frame(matrix(ncol = 9, nrow = 0))
Data_Alt_2 <- Data_Alt_1
names(Data_Alt_1) <- names(Data_Alt)
names(Data_Alt_2) <- names(Data_Alt)
Data_sub <- subset(Data, Week <= 30)
for(i in 1:nrow(Data_sub)){
  if(i < nrow(Data_sub)){
    if(Data_sub$Week[i] >= Data_sub$Week[i + 1]){
      Data_Alt_1[i, 1:8] <- Data_sub[i, 1:8]}}
  if(i == nrow(Data_sub)){
    Data_Alt_1[i, 1:8] <- Data_sub[i, 1:8]}
  if(Data_sub$Week[i] == 30){Data_Alt_1[i, 9] <- 0} else {Data_Alt_1[i, 9] <- 1}}
Data_sub <- subset(Data, Week >= 50)
for(i in 1:nrow(Data_sub)){
  if(i < nrow(Data_sub)){
    if(Data_sub$Week[i] >= Data_sub$Week[i + 1]){
      Data_Alt_2[i, 1:8] <- Data_sub[i, 1:8]}}
  if(i == nrow(Data_sub)){
    Data_Alt_2[i, 1:8] <- Data_sub[i, 1:8]}
  if(Data_sub$Week[i] == 50){Data_Alt_2[i, 9] <- 0} else {Data_Alt_2[i, 9] <- 1}}
Data_Alt_1 <- drop_na(Data_Alt_1, Row)
Data_Alt_2 <- drop_na(Data_Alt_2, Row)

# Remove temporary variables since they will no longer be used
remove(i, ns_prev, Data_sub)

# Construct dataframes used in growth models [WIP]
Data %>% 
  group_by(Row, Group, Plant, Species, Warmed, Treatment, DM_t) %>% 
  summarise(HG_TA = mean(HGain),
            SG_TA = mean(SGain)) -> Data_TA
names(Data_TA) <- names(Data)[c(1:7, 10, 14)]
Data %>% 
  filter(Week <= 30) %>% 
  group_by(Row, Group, Plant, Species, Warmed, Treatment, DM_t) %>% 
  summarise(HG_TA = mean(HGain),
            SG_TA = mean(SGain)) -> Data_TA_1
names(Data_TA_1) <- names(Data)[c(1:7, 10, 14)]
Data %>% 
  filter(Week >= 50) %>% 
  group_by(Row, Group, Plant, Species, Warmed, Treatment, DM_t) %>% 
  summarise(HG_TA = mean(HGain),
            SG_TA = mean(SGain)) -> Data_TA_2
names(Data_TA_2) <- names(Data)[c(1:7, 10, 14)]





##### Set up plotting functions ---------------------------------------------------------------------------

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

