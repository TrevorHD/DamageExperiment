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
Data$SGain <- rep(0, nrow(Data))
for(i in 1:length(1:nrow(Data))){
  if(Data$Week[i] != 0){
    Data$SGain[i] <- Data$NStems[i] - Data$NStems[i - 1]}}

# Select columns used for analyses  
Data %>% select(Row, Group, Plant, Species, Warmed, Treatment, DM_t, Week, Height,
                HGain, Buds, Heads, NStems, SGain) -> Data

