##### [Supporting] Set up function to approximate dispersion parameter ------------------------------------

# Function from B. Bolker (lme4 co-author) to approximate and test the dispersion factor
# Source: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
# Refer to source link for caveats and details
dispersion <- function(model){
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)}





##### [Supporting] Fit bud count models with Poisson regression -------------------------------------------

# Select data for max bud count analysis
# Note: unlike main analysis, keeping all individuals rather than just ones that budded
# Otherwise it would ruin any potential Poisson distribution and this exercise would be moot
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, sans interactions, and sans interactions/size
Mod_Form1 <- MaxBuds ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- MaxBuds ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- MaxBuds ~ factor(Warmed) + factor(Treatment) + (1|PlotID)

# For CN, model max bud count as function of trimming and warming treatments
Mod_Buds1_CN_P <- glmer(Mod_Form1, data = Data_Temp_CN, family = "poisson")
summary(Mod_Buds1_CN_P)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_Buds2_CN_P <- glmer(Mod_Form2, data = Data_Temp_CN, family = "poisson")
summary(Mod_Buds2_CN_P)
anova(Mod_Buds1_CN_P, Mod_Buds2_CN_P)

# Fit model without rosette size; rosette size is significant, so keep it
Mod_Buds3_CN_P <- glmer(Mod_Form3, data = Data_Temp_CN, family = "poisson")
summary(Mod_Buds3_CN_P)
anova(Mod_Buds2_CN_P, Mod_Buds3_CN_P)

# Chi-squared test on dispersion doesn't show a significant difference from 1
dispersion(Mod_Buds2_CN_P)

# GOF tests show that bud count does not reliably follow a Poisson distribution
# Seems to be deviation from Poisson in aggregate, and within select trimming and warming treatments
temp1 <- subset(Data_Sum, Species == "CN")
summary(goodfit(subset(temp1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 2)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 3)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 0)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 1)$MaxBuds, "poisson"))

# For CA, model max bud count as function of trimming and warming treatments
Mod_Buds1_CA_P <- glmer(Mod_Form1, data = Data_Temp_CA, family = "poisson")
summary(Mod_Buds1_CA_P)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_Buds2_CA_P <- glmer(Mod_Form2, data = Data_Temp_CA, family = "poisson")
summary(Mod_Buds2_CA_P)
anova(Mod_Buds1_CA_P, Mod_Buds2_CA_P)

# Fit model without rosette size; rosette size is not significant, so drop it
Mod_Buds3_CA_P <- glmer(Mod_Form3, data = Data_Temp_CA, family = "poisson")
summary(Mod_Buds3_CA_P)
anova(Mod_Buds2_CA_P, Mod_Buds3_CA_P)

# Chi-squared test on dispersion doesn't show a significant difference from 1
dispersion(Mod_Buds3_CA_P)

# GOF tests show that bud count does not reliably follow a Poisson distribution
# Seems to be deviation from Poisson in aggregate, and within select trimming and warming treatments
temp1 <- subset(Data_Sum, Species == "CA")
summary(goodfit(subset(temp1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 2)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 3)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 0)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 1)$MaxBuds, "poisson"))

# Remove temporary variables that are no longer needed
remove(temp1, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Buds1_CN_P, Mod_Buds3_CN_P,
       Mod_Buds1_CA_P, Mod_Buds2_CA_P)





##### [Supporting] Fit stem count models with Poisson regression ------------------------------------------

# Select data for max stem count analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, sans interactions, and sans interactions/size
Mod_Form1 <- MaxS ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- MaxS ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- MaxS ~ factor(Warmed) + factor(Treatment) + (1|PlotID)

# For CN, model max stem count as function of trimming and warming treatments
Mod_SMax1_CN_P <- glmer(Mod_Form1, data = Data_Temp_CN, family = "poisson")
summary(Mod_SMax1_CN_P)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_SMax2_CN_P <- glmer(Mod_Form2, data = Data_Temp_CN, family = "poisson")
summary(Mod_SMax2_CN_P)
anova(Mod_SMax1_CN_P, Mod_SMax2_CN_P)

# Fit model without rosette size; rosette size is significant, so keep it
Mod_SMax3_CN_P <- glmer(Mod_Form3, data = Data_Temp_CN, family = "poisson")
summary(Mod_SMax3_CN_P)
anova(Mod_SMax2_CN_P, Mod_SMax3_CN_P)

# Approximation of dispersion factor
# Chi-squared test doesn't show a significant difference from 1
dispersion(Mod_SMax2_CN_P)

# GOF tests show that stem count does not reliably follow a Poisson distribution
# Seems to be deviation from Poisson in aggregate, and within select trimming and warming treatments
temp1 <- subset(Data_Sum, Species == "CN")
summary(goodfit(subset(temp1)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 1)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 2)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 3)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 4)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Warmed == 0)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Warmed == 1)$MaxS, "poisson"))

# For CA, model max stem count as function of trimming and warming treatments
Mod_SMax1_CA_P <- glmer(Mod_Form1, data = Data_Temp_CA, family = "poisson")
summary(Mod_SMax1_CA_P)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_SMax2_CA_P <- glmer(Mod_Form2, data = Data_Temp_CA, family = "poisson")
summary(Mod_SMax2_CA_P)
anova(Mod_SMax1_CA_P, Mod_SMax2_CA_P)

# Fit model without rosette size; rosette size is not significant, so drop it
Mod_SMax3_CA_P <- glmer(Mod_Form3, data = Data_Temp_CA, family = "poisson")
summary(Mod_SMax3_CA_P)
anova(Mod_SMax2_CA_P, Mod_SMax3_CA_P)

# Chi-squared test on dispersion doesn't show a significant difference from 1
dispersion(Mod_SMax3_CA_P)

# Again, Poisson models may not be a good fit, as GOF tests show deviation from Poisson distribution
# Seems to be deviation from Poisson in aggregate, and within select trimming and warming treatments
temp1 <- subset(Data_Sum, Species == "CA")
summary(goodfit(subset(temp1)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 1)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 2)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 3)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Treatment == 4)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Warmed == 0)$MaxS, "poisson"))
summary(goodfit(subset(temp1, Warmed == 1)$MaxS, "poisson"))

# Remove temporary variables that are no longer needed
remove(temp1, Mod_Form1, Mod_Form2, Mod_Form3, Mod_SMax1_CN_P, Mod_SMax3_CN_P,
       Mod_SMax1_CA_P, Mod_SMax2_CA_P)

