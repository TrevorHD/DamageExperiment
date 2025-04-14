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

# For CN, model max bud count as function of trimming and warming treatments
Mod_Buds1_CN_P <- glmer(MaxBuds ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CN"), family = "poisson")
summary(Mod_Buds1_CN_P)

# Fit the same model as above, but without the interaction term
# ANOVA shows that the interaction term is not significant, so we drop it
Mod_Buds2_CN_P <- glmer(MaxBuds ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CN"), family = "poisson")
summary(Mod_Buds2_CN_P)
anova(Mod_Buds1_CN_P, Mod_Buds2_CN_P)

# Approximation of dispersion factor
# Chi-squared test doesn't show a significant difference from 1
dispersion(Mod_Buds2_CN_P)

# While models converge, Poisson distribution may not be a good fit
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
Mod_Buds1_CA_P <- glmer(MaxBuds ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CA"), family = "poisson")
summary(Mod_Buds1_CA_P)

# Fit the same model as above, but without the interaction term
# ANOVA shows that the interaction term is not significant, so we drop it
Mod_Buds2_CA_P <- glmer(MaxBuds ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CA"), family = "poisson")
summary(Mod_Buds2_CA_P)
anova(Mod_Buds1_CA_P, Mod_Buds2_CA_P)

# Approximation of dispersion factor
# Chi-squared test doesn't show a significant difference from 1
dispersion(Mod_Buds2_CA_P)

# Again, Poisson models may not be a good fit, as GOF tests show deviation from Poisson distribution
# Seems to be deviation from Poisson in aggregate, and within select trimming and warming treatments
temp1 <- subset(Data_Sum, Species == "CA")
summary(goodfit(subset(temp1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 1)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 2)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Treatment == 3)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 0)$MaxBuds, "poisson"))
summary(goodfit(subset(temp1, Warmed == 1)$MaxBuds, "poisson"))

# Remove intermediate models that are no longer used
remove(Mod_Buds1_CN_P, Mod_Buds1_CA_P)





##### [Supporting] Fit stem count models with Poisson regression ------------------------------------------

# For CN, model max stem count as function of trimming and warming treatments
Mod_SMax1_CN_P <- glmer(MaxS ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CN"), family = "poisson")
summary(Mod_SMax1_CN_P)

# Fit the same model as above, but without the interaction term
# ANOVA shows that the interaction term is not significant, so we drop it
Mod_SMax2_CN_P <- glmer(MaxS ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CN"), family = "poisson")
summary(Mod_SMax2_CN_P)
anova(Mod_SMax1_CN_P, Mod_SMax2_CN_P)

# Approximation of dispersion factor
# Chi-squared test doesn't show a significant difference from 1
dispersion(Mod_SMax2_CN_P)

# While models converge, Poisson distribution may not actually be a good fit
# GOF tests show that bud count does not reliably follow a Poisson distribution
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
Mod_SMax1_CA_P <- glmer(MaxS ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CA"), family = "poisson")
summary(Mod_SMax1_CA_P)

# Fit the same model as above, but without the interaction term
# ANOVA shows that the interaction term is not significant, so we drop it
Mod_SMax2_CA_P <- glmer(MaxS ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                        data = subset(Data_Sum, Species == "CA"), family = "poisson")
summary(Mod_SMax2_CA_P)
anova(Mod_SMax1_CA_P, Mod_SMax2_CA_P)

# Approximation of dispersion factor
# Chi-squared test doesn't show a significant difference from 1
dispersion(Mod_SMax2_CA_P)

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

# Remove intermediate models that are no longer used
remove(Mod_SMax1_CN_P, Mod_SMax1_CA_P)

