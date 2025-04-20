##### Fit budding probability models ----------------------------------------------------------------------

# Select data for flowering probability analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formula for full model
Mod_Form1 <- Budded ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)

# For CN, create tables of number of individuals that produced at least one bud
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
table(Data_Temp_CN$Treatment, Data_Temp_CN$Budded)
table(Data_Temp_CN$Warmed, Data_Temp_CN$Budded)
table(Data_Temp_CN$Treatment, Data_Temp_CN$Budded, Data_Temp_CN$Warmed)

# The aforementioned separation causes boundary issues and greatly inflates standard errors
Mod_BudP_CN <- glmer(Mod_Form1, data = Data_Temp_CN, family = "binomial")
summary(Mod_BudP_CN)

# For CA, create same tables/breakouts as those done for CN
table(Data_Temp_CA$Treatment, Data_Temp_CA$Budded)
table(Data_Temp_CA$Warmed, Data_Temp_CA$Budded)
table(Data_Temp_CA$Treatment, Data_Temp_CA$Budded, Data_Temp_CA$Warmed)

# Extreme separation again causes boundary issues and greatly inflates standard errors
Mod_BudP_CA <- glmer(Mod_Form1, data = Data_Temp_CA, family = "binomial")
summary(Mod_BudP_CA)

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_BudP_CN, Mod_BudP_CA)





##### Fit flowering probability models --------------------------------------------------------------------

# Select data for flowering probability analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formula for full model
Mod_Form1 <- Flowered ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)

# For CN, create tables of number of individuals that produced at least one flower
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
# Note: for the one trimmed individual that flowered, that one flower was destroyed a week later
table(Data_Temp_CN$Treatment, Data_Temp_CN$Flowered)
table(Data_Temp_CN$Warmed, Data_Temp_CN$Flowered)
table(Data_Temp_CN$Treatment, Data_Temp_CN$Flowered, Data_Temp_CN$Warmed)

# The aforementioned separation results in the model failing to converge
Mod_FloP_CN <- glmer(Mod_Form1, data = Data_Temp_CN, family = "binomial")
summary(Mod_FloP_CN)

# For CA, create same tables/breakouts as those done for CN
table(Data_Temp_CA$Treatment, Data_Temp_CA$Flowered)
table(Data_Temp_CA$Warmed, Data_Temp_CA$Flowered)
table(Data_Temp_CA$Treatment, Data_Temp_CA$Flowered, Data_Temp_CA$Warmed)

# Extreme separation again causes boundary issues and greatly inflates standard errors
Mod_FloP_CA <- glmer(Mod_Form1, data = Data_Temp_CA, family = "binomial")
summary(Mod_FloP_CA)

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_FloP_CN, Mod_FloP_CA)





##### Fit max bud count models ----------------------------------------------------------------------------

# Select data for max bud count analysis
# Note: data restricted to individuals that budded (>0 buds)
Data_Temp_CN <- subset(Data_Sum, Species == "CN" & MaxBuds > 0)
Data_Temp_CA <- subset(Data_Sum, Species == "CA" & MaxBuds > 0)

# Create formulas for full model, log full model, and log sans interaction
Mod_Form1 <- MaxBuds ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxBuds) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxBuds) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)

# For CN, examine distribution of max bud count (>0 buds) by trimming and warming treatments
# Unsurprisingly, untrimmed plants tend to have more buds than trimmed ones
# No obvious differences for warmed vs unwarmed plants, though
dd.plot(Data_Temp_CN, "MaxBuds", "Treatment", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxBuds", "Treatment", "dens", c(0, 15))
dd.plot(Data_Temp_CN, "MaxBuds", "Warmed", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxBuds", "Warmed", "dens", c(0, 15))

# For CN, model max bud count as function of trimming and warming treatments
Mod_Buds1_CN <- lmer(Mod_Form1, data = Data_Temp_CN)
summary(Mod_Buds1_CN)

# Shapiro-Wilk rejects residual normality, likely due to right tail
# Residuals are heteroskedastic, mostly due to the control group
qqnorm(resid(Mod_Buds1_CN))
qqline(resid(Mod_Buds1_CN))
shapiro.test(resid(Mod_Buds1_CN))
data.frame("Fitted" = fitted(Mod_Buds1_CN), "Resid" = resid(Mod_Buds1_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_Buds2_CN <- lmer(Mod_Form2, data = Data_Temp_CN)
summary(Mod_Buds2_CN)

# Seems to work; residuals appear normal with little to no heteroskedasticity
qqnorm(resid(Mod_Buds2_CN))
qqline(resid(Mod_Buds2_CN))
shapiro.test(resid(Mod_Buds2_CN))
data.frame("Fitted" = fitted(Mod_Buds2_CN), "Resid" = resid(Mod_Buds2_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_Buds3_CN <- lmer(Mod_Form3, data = Data_Temp_CN)
summary(Mod_Buds3_CN)
anova(Mod_Buds2_CN, Mod_Buds3_CN)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_Buds3_CN))
qqline(resid(Mod_Buds3_CN))
shapiro.test(resid(Mod_Buds3_CN))
data.frame("Fitted" = fitted(Mod_Buds3_CN), "Resid" = resid(Mod_Buds3_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# For CA, similar to CN, untrimmed plants tend to have more buds than trimmed ones
# No obvious differences for warmed vs unwarmed plants, though
dd.plot(Data_Temp_CA, "MaxBuds", "Treatment", "hist", c(0, 50))
dd.plot(Data_Temp_CA, "MaxBuds", "Treatment", "dens", c(0, 50))
dd.plot(Data_Temp_CA, "MaxBuds", "Warmed", "hist", c(0, 50))
dd.plot(Data_Temp_CA, "MaxBuds", "Warmed", "dens", c(0, 50))

# For CA, model max bud count as function of trimming and warming treatments
Mod_Buds1_CA <- lmer(Mod_Form1, data = Data_Temp_CA)
summary(Mod_Buds1_CA)

# Residuals appear normally distributed, but heteroskedastic
qqnorm(resid(Mod_Buds1_CA))
qqline(resid(Mod_Buds1_CA))
shapiro.test(resid(Mod_Buds1_CA))
data.frame("Fitted" = fitted(Mod_Buds1_CA), "Resid" = resid(Mod_Buds1_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_Buds2_CA <- lmer(Mod_Form2, data = Data_Temp_CA)
summary(Mod_Buds2_CA)

# Seems to work; residuals appear normal with little to no heteroskedasticity
qqnorm(resid(Mod_Buds2_CA))
qqline(resid(Mod_Buds2_CA))
shapiro.test(resid(Mod_Buds2_CA))
data.frame("Fitted" = fitted(Mod_Buds2_CA), "Resid" = resid(Mod_Buds2_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_Buds3_CA <- lmer(Mod_Form3, data = Data_Temp_CA)
summary(Mod_Buds3_CA)
anova(Mod_Buds2_CA, Mod_Buds3_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_Buds3_CA))
qqline(resid(Mod_Buds3_CA))
shapiro.test(resid(Mod_Buds3_CA))
data.frame("Fitted" = fitted(Mod_Buds3_CA), "Resid" = resid(Mod_Buds3_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Note: residual plots above will show banding patterns due to discrete nature of data
# Boundary effects also visible on bottom left of plots since log(y) cannot be less than zero
# Neither of these are a cause for concern, nor do they violate model assumptions

# Pairwise comparison of (marginal mean) max bud count between trimming treatments
# For each species, trimmed plants produce significantly fewer buds than the control (as seen earlier)
# Trimming treatments 2/3 within a species do not differ from each other, though
pairs(emmeans(Mod_Buds3_CN, "Treatment"))
pairs(emmeans(Mod_Buds3_CA, "Treatment"))

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3,
       Mod_Buds1_CN, Mod_Buds2_CN, Mod_Buds1_CA, Mod_Buds2_CA)





##### Fit max height models -------------------------------------------------------------------------------

# Select data for max height analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, log full model, and log sans interaction
Mod_Form1 <- MaxH ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxH) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxH) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)

# For CN, examine distribution of max height by trimming and warming treatments
# Unsurprisingly, untrimmed plants tend to be taller than trimmed ones
# Might be a difference between warmed and unwarmed, but need to investigate further
dd.plot(Data_Temp_CN, "MaxH", "Treatment", "hist", c(0, 200))
dd.plot(Data_Temp_CN, "MaxH", "Treatment", "dens", c(0, 200))
dd.plot(Data_Temp_CN, "MaxH", "Warmed", "hist", c(0, 200))
dd.plot(Data_Temp_CN, "MaxH", "Warmed", "dens", c(0, 200))

# For CN, model max height as function of trimming and warming treatments
Mod_HMax1_CN <- lmer(Mod_Form1, data = Data_Temp_CN)
summary(Mod_HMax1_CN)

# Shapiro-Wilk rejects residual normality, likely due to tails
# Residuals are heteroskedastic, mostly due to the control group
qqnorm(resid(Mod_HMax1_CN))
qqline(resid(Mod_HMax1_CN))
shapiro.test(resid(Mod_HMax1_CN))
data.frame("Fitted" = fitted(Mod_HMax1_CN), "Resid" = resid(Mod_HMax1_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_HMax2_CN <- lmer(Mod_Form2, data = Data_Temp_CN)
summary(Mod_HMax2_CN)

# Seems to work; residuals appear normal with little to no heteroskedasticity
qqnorm(resid(Mod_HMax2_CN))
qqline(resid(Mod_HMax2_CN))
shapiro.test(resid(Mod_HMax2_CN))
data.frame("Fitted" = fitted(Mod_HMax2_CN), "Resid" = resid(Mod_HMax2_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_HMax3_CN <- lmer(Mod_Form3, data = Data_Temp_CN)
summary(Mod_HMax2_CN); summary(Mod_HMax3_CN)
anova(Mod_HMax2_CN, Mod_HMax3_CN)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_HMax3_CN))
qqline(resid(Mod_HMax3_CN))
shapiro.test(resid(Mod_HMax3_CN))
data.frame("Fitted" = fitted(Mod_HMax3_CN), "Resid" = resid(Mod_HMax3_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# For CA, similar to CN, untrimmed plants tend to be taller than trimmed ones
# Might be a difference between warmed and unwarmed, but need to investigate further
dd.plot(Data_Temp_CA, "MaxH", "Treatment", "hist", c(0, 200))
dd.plot(Data_Temp_CA, "MaxH", "Treatment", "dens", c(0, 200))
dd.plot(Data_Temp_CA, "MaxH", "Warmed", "hist", c(0, 200))
dd.plot(Data_Temp_CA, "MaxH", "Warmed", "dens", c(0, 200))

# For CA, model max height as function of trimming and warming treatments
Mod_HMax1_CA <- lmer(Mod_Form1, data = Data_Temp_CA)
summary(Mod_HMax1_CA)

# Shapiro-Wilk rejects residual normality, likely due to tails
# Residuals are heteroskedastic, mostly due to the control group
qqnorm(resid(Mod_HMax1_CA))
qqline(resid(Mod_HMax1_CA))
shapiro.test(resid(Mod_HMax1_CA))
data.frame("Fitted" = fitted(Mod_HMax1_CA), "Resid" = resid(Mod_HMax1_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_HMax2_CA <- lmer(Mod_Form2, data = Data_Temp_CA)

# Seems to work; residuals appear normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on tails; not a cause for concern
summary(Mod_HMax2_CA)
qqnorm(resid(Mod_HMax2_CA))
qqline(resid(Mod_HMax2_CA))
shapiro.test(resid(Mod_HMax2_CA))
data.frame("Fitted" = fitted(Mod_HMax2_CA), "Resid" = resid(Mod_HMax2_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_HMax3_CA <- lmer(Mod_Form3, data = Data_Temp_CA)
summary(Mod_HMax2_CA); summary(Mod_HMax3_CA)
anova(Mod_HMax2_CA, Mod_HMax3_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on left; not a cause for concern
qqnorm(resid(Mod_HMax3_CA))
qqline(resid(Mod_HMax3_CA))
shapiro.test(resid(Mod_HMax3_CA))
data.frame("Fitted" = fitted(Mod_HMax3_CA), "Resid" = resid(Mod_HMax3_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Pairwise comparison of (marginal mean) max height between trimming treatments
# For both species, trimmed plants are significantly shorter than the control (as seen earlier)
# The only trimming treatments that different from each other were 2/3 for CN and 2/4 for CA
pairs(emmeans(Mod_HMax3_CN, "Treatment"))
pairs(emmeans(Mod_HMax3_CA, "Treatment"))

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3,
       Mod_HMax1_CN, Mod_HMax2_CN,  Mod_HMax1_CA, Mod_HMax2_CA)





##### Fit max stem count models ---------------------------------------------------------------------------

# Select data for max stem count analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, log full model, and log sans interaction
Mod_Form1 <- MaxS ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxS) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxS) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)

# For CN, examine distribution of max stem count by trimming and warming treatments
# Trimming appears to disrupt apical dominance and results in more stems
# No obvious differences for warmed vs unwarmed plants, though
dd.plot(Data_Temp_CN, "MaxS", "Treatment", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Treatment", "dens", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Warmed", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Warmed", "dens", c(0, 15))

# For CN, model max stem count as function of trimming and warming treatments
Mod_SMax1_CN <- lmer(Mod_Form1, data = Data_Temp_CN)
summary(Mod_SMax1_CN)

# Residuals appear normally distributed, but heteroskedastic
qqnorm(resid(Mod_SMax1_CN))
qqline(resid(Mod_SMax1_CN))
shapiro.test(resid(Mod_SMax1_CN))
data.frame("Fitted" = fitted(Mod_SMax1_CN), "Resid" = resid(Mod_SMax1_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_SMax2_CN <- lmer(Mod_Form2, data = Data_Temp_CN)
summary(Mod_SMax2_CN)

# Seems to work; residuals appear normal with little to no heteroskedasticity
# Shapiro-Wilk rejects normality on tails; not a cause for concern
qqnorm(resid(Mod_SMax2_CN))
qqline(resid(Mod_SMax2_CN))
shapiro.test(resid(Mod_SMax2_CN))
data.frame("Fitted" = fitted(Mod_SMax2_CN), "Resid" = resid(Mod_SMax2_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_SMax3_CN <- lmer(Mod_Form3, data = Data_Temp_CN)
summary(Mod_SMax3_CN)
anova(Mod_SMax2_CN, Mod_SMax3_CN)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on left tail; not a cause for concern
qqnorm(resid(Mod_SMax3_CN))
qqline(resid(Mod_SMax3_CN))
shapiro.test(resid(Mod_SMax3_CN))
data.frame("Fitted" = fitted(Mod_SMax3_CN), "Resid" = resid(Mod_SMax3_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# For CA, similar to CN, trimming results in more stems
# No obvious differences for warmed vs unwarmed plants, though
dd.plot(temp2, "MaxS", "Treatment", "hist", c(0, 15))
dd.plot(temp2, "MaxS", "Treatment", "dens", c(0, 15))
dd.plot(temp2, "MaxS", "Warmed", "hist", c(0, 15))
dd.plot(temp2, "MaxS", "Warmed", "dens", c(0, 15))

# For CA, model max stem count as function of trimming and warming treatments
Mod_SMax1_CA <- lmer(Mod_Form1, data = Data_Temp_CA)
summary(Mod_SMax1_CA)

# Residuals appear normally distributed, but heteroskedastic
qqnorm(resid(Mod_SMax1_CA))
qqline(resid(Mod_SMax1_CA))
shapiro.test(resid(Mod_SMax1_CA))
data.frame("Fitted" = fitted(Mod_SMax1_CA), "Resid" = resid(Mod_SMax1_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Log-transform model to fix heteroskedasticity and ensure preds >0
Mod_SMax2_CA <- lmer(Mod_Form2, data = Data_Temp_CA)
summary(Mod_SMax2_CA)

# Seems to work; residuals appear normal with little to no heteroskedasticity
qqnorm(resid(Mod_SMax2_CA))
qqline(resid(Mod_SMax2_CA))
shapiro.test(resid(Mod_SMax2_CA))
data.frame("Fitted" = fitted(Mod_SMax2_CA), "Resid" = resid(Mod_SMax2_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without interactions; interaction term is not significant, so drop it
Mod_SMax3_CA <- lmer(Mod_Form3, data = Data_Temp_CA)
summary(Mod_SMax3_CA)
anova(Mod_SMax2_CA, Mod_SMax3_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_SMax3_CA))
qqline(resid(Mod_SMax3_CA))
shapiro.test(resid(Mod_SMax3_CA))
data.frame("Fitted" = fitted(Mod_SMax3_CA), "Resid" = resid(Mod_SMax3_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Note: residual plots above will show banding patterns due to discrete nature of data
# Boundary effects also visible on bottom left of plots since log(y) cannot be less than zero
# Neither of these are a cause for concern, nor do they violate model assumptions

# Pairwise comparison of (marginal mean) max stem count between trimming treatments
# For both species, trimmed plants have significantly more stems than the control (as seen earlier)
# But none of the trimming treatments are significantly different from each other
pairs(emmeans(Mod_SMax3_CN, "Treatment"))
pairs(emmeans(Mod_SMax3_CA, "Treatment"))

# Remove intermediate models that are no longer used
remove(temp1, temp2, Mod_SMax1_CN, Mod_SMax2_CN, Mod_SMax1_CA, Mod_SMax2_CA)





##### Fit survival regressions ----------------------------------------------------------------------------

# Select data for flowering probability analysis
Data_Temp_CN_Y1 <- subset(Data_Alt_Y1, Species == "CN")
Data_Temp_CA_Y1 <- subset(Data_Alt_Y1, Species == "CA")
Data_Temp_CN_Y2 <- subset(Data_Alt_Y2, Species == "CN")
Data_Temp_CA_Y2 <- subset(Data_Alt_Y2, Species == "CA")

# For CN (Year 1), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CN_Y1$Treatment, Data_Temp_CN_Y1$Cens)
table(Data_Temp_CN_Y1$Warmed, Data_Temp_CN_Y1$Cens)
table(Data_Temp_CN_Y1$Treatment, Data_Temp_CN_Y1$Cens, Data_Temp_CN_Y1$Warmed)

# For CN (Year 1), model hazard as function of trimming and warming treatments
# 2 out of 60 individuals survived until end of trimming season
Mod_Surv1_CN <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CN"))

# Fit the same model above, but without interaction effect
Mod_Surv2_CN <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CN"))

# Examine model outputs
summary(Mod_Surv1_CN)
summary(Mod_Surv2_CN)

# Interaction term is not significant, so drop it
anova(Mod_Surv2_CN, Mod_Surv1_CN)

# Examine diagnostics of full model
# Note: hazard for treatment is mostly constant, but decreases near the end
# This is what is likely driving a significant value here for treatment
# Not a cause for concern
cox.zph(Mod_Surv2_CN)
plot(cox.zph(Mod_Surv2_CN))

# For CA (Year 1), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CA_Y1$Treatment, Data_Temp_CA_Y1$Cens)
table(Data_Temp_CA_Y1$Warmed, Data_Temp_CA_Y1$Cens)
table(Data_Temp_CA_Y1$Treatment, Data_Temp_CA_Y1$Cens, Data_Temp_CA_Y1$Warmed)

# For CA (Year 1), model hazard as function of trimming and warming treatments
# 23 out of 60 individuals survived until end of trimming season
Mod_Surv1_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CA"))

# Fit the same model above, but without interaction effect
Mod_Surv2_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CA"))

# Examine model outputs
summary(Mod_Surv1_CA)
summary(Mod_Surv2_CA)

# Interaction term is significant, so keep it
anova(Mod_Surv2_CA, Mod_Surv1_CA)

# Examine diagnostics of reduced model; looks good
cox.zph(Mod_Surv1_CA)
plot(cox.zph(Mod_Surv1_CA))

# For CA (Year 2), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CA_Y2$Treatment, Data_Temp_CA_Y2$Cens)
table(Data_Temp_CA_Y2$Warmed, Data_Temp_CA_Y2$Cens)
table(Data_Temp_CA_Y2$Treatment, Data_Temp_CA_Y2$Cens, Data_Temp_CA_Y2$Warmed)

# For CA (Year 2), model hazard as function of trimming and warming treatments
# 3 out of 19 individuals survived until end of trimming season
Mod_Surv3_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y2, Species == "CA"))

# Fit the same model above, but without interaction effect
Mod_Surv4_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y2, Species == "CA"))

# Examine model outputs
summary(Mod_Surv3_CA)
summary(Mod_Surv4_CA)

# Interaction term is not significant, so drop it
anova(Mod_Surv4_CA, Mod_Surv3_CA)

# Examine diagnostics of reduced model; looks good
cox.zph(Mod_Surv2_CA)
plot(cox.zph(Mod_Surv2_CA))

