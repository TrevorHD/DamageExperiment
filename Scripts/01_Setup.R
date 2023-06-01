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



##### [F1] Plot survival curves by trimming treatment -----------------------------------------------------

# Generate Kaplan-Meier survival curves for comparison by treatment
# Note: not using the full time horizon for CN because most of them died before winter
km_CN_NW <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Warmed == 0))
km_CN_W <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                   data = subset(Data_Alt_1, Species == "CN" & Warmed == 1))
km_CA_NW <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Warmed == 0))
km_CA_W <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                   data = subset(Data_Alt, Species == "CA" & Warmed == 1))

# Function to plot survival curves
km.plot <- function(km_mod, bottom, left, atext){
  ltypes <- c("solid", "3333", "2323", "1212")
  colours <- c("black", "purple", "green", "orange")
  par(cex.axis = 0.38, cex.lab = 0.5, tcl = -0.15)
  if(bottom == FALSE & left == TRUE){
    par(mar = c(0, 0.8, 0.15, 0.1))
    plot(km_mod, xlim = c(0, 30), xaxt = "n", ylab = "Probability of Survival",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == TRUE & left == TRUE){
    par(mar = c(1, 0.8, 0, 0.1))
    plot(km_mod, xlim = c(0, 30), xlab = "Weeks", ylab = "Probability of Survival", 
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    text(x = 30, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == FALSE & left == FALSE){
    par(mar = c(0, 0.1, 0.15, 0.3))
    plot(km_mod, xaxt = "n", yaxt = "n",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}
  if(bottom == TRUE & left == FALSE){
    par(mar = c(1, 0.1, 0, 0.3))
    plot(km_mod, xlab = "Weeks", yaxt = "n",
         mgp = c(0.2, -0.3, 0), lwd = 1.1, lty = ltypes, col = colours)
    text(x = 66.2, y = 0.98, atext, adj = 1, cex = 0.39)}}

# Prepare graphics device
tiff(filename = "Figure 1.tif", width = 2700, height = 2000, units = "px", res = 800, compression = "lzw")

# Create blank page
grid.newpage()
plot.new()

# Set grid layout and activate it
gly <- grid.layout(2000, 2700)
pushViewport(viewport(layout = gly))

# Plot unwarmed CN
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CN_NW, bottom = FALSE, left = TRUE, atext = "CN Unwarmed"))
popViewport()

# Plot warmed CN
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 50:975))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CN_W, bottom = TRUE, left = TRUE, atext = "CN Warmed"))
popViewport()

# Plot unwarmed CA
pushViewport(vp = viewport(layout.pos.row = 25:925, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CA_NW, bottom = FALSE, left = FALSE, atext = "CA Unwarmed"))
popViewport()

# Plot warmed CA
pushViewport(vp = viewport(layout.pos.row = 975:1975, layout.pos.col = 1000:2675))
par(fig = gridFIG())
par(new = TRUE)
print(km.plot(km_CA_W, bottom = TRUE, left = FALSE, atext = "CA Warmed"))
popViewport()

# Create legend
grid.text(label = c("Control", "Trim to 10 cm", "Trim to 5 cm", "Trim to Ground"),
          x = rep(0.908, 4), y = seq(0.913, 0.840, length.out = 4),
          hjust = 1, gp = gpar(cex = 0.35))
grid.segments(x0 = rep(0.920, 4), y0 = seq(0.913, 0.840, length.out = 4),
              x1 = rep(0.956, 4), y1 = seq(0.913, 0.840, length.out = 4),
              gp = gpar(lty = c(1, 5, 2, 3), lty = c("solid", "3333", "2323", "1212"), 
                        col = c("black", "purple", "green", "orange")))

# Deactivate grid layout; finalise graphics save
popViewport()
dev.off()



















km5 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt_1, Species == "CN" & Treatment == 1))
km6 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt_1, Species == "CN" & Treatment == 2))
km7 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt_1, Species == "CN" & Treatment == 3))
km8 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt_1, Species == "CN" & Treatment == 4))
plot(km5, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km6, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km7, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km8, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")


km9 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt, Species == "CA" & Treatment == 1))
km10 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt, Species == "CA" & Treatment == 2))
km11 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt, Species == "CA" & Treatment == 3))
km12 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier", data = subset(Data_Alt, Species == "CA" & Treatment == 4))
plot(km9, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km10, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km11, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")
plot(km12, xlab = "Weeks", ylab = "Probability of Survival", col = c("blue", "red"))
legend ("topright", legend = c("Unwarmed", "Warmed"), fill = c("blue", "red"), bty = "n")