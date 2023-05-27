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



# note: figure out what to do with the ToD>0 bit
Surv1_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CN" & ToD > 0))
summary(Surv1_CN)
AIC(Surv1_CN)


# note: figure out what to do with the ToD>0 bit
Surv1_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CA" & ToD > 0))
summary(Surv1_CA)
AIC(Surv1_CA)


# Models for regrowth
mod1 <- lm(HGain ~ factor(Warmed) + factor(Treatment) + DM_t, data = subset(Data_TA, Species == "CN"))
summary(mod1)
mod2 <- lm(HGain ~ factor(Warmed) + factor(Treatment) + DM_t, data = subset(Data_TA, Species == "CA"))
summary(mod2)

# Models for stem count regrowth
mod1 <- lm(SGain ~ factor(Warmed) + factor(Treatment) + DM_t, data = subset(Data_TA, Species == "CN"))
summary(mod1)
mod2 <- lm(SGain ~ factor(Warmed) + factor(Treatment) + DM_t, data = subset(Data_TA, Species == "CA"))
summary(mod2)







# Survival analysis
# Models comparing probability of additional stems for each treatment combination
# Models comparing number of additional stems (only for individuals where observed for each treatment combination)
# Models comparing probability of reproductive structures for each treatment combination
# Models comparing number of reproductive structures (only for individuals where observed for each treatment combination)


Data %>%
  group_by(Row, Group, Plant, Species, Warmed, Treatment, DM_t) %>%
  summarise(meanG = mean(HGain), sdG = sd(HGain)) -> test
test %>% 
  filter(Species == "CN", Warmed == 0, Treatment == 1)







