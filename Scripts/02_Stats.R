##### Fit survival regressions ----------------------------------------------------------------------------

# Note: interactions excluded as statistical power is insufficient due to small sample size
# Using Cox model since we don't need to predict outside of observation period
# Plus, unlike survreg, Cox supports random effects

# Fit survival model, using initial diameter as covariate [CN]
Surv1_CN <- coxph(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) +
                    frailty(PlotID, "gaussian"),
                  data = subset(Data_Alt, Species == "CN"))

# View model results
Surv1_CN

# Rosette diameter not useful as a covariate; drop and re-fit model
Surv2_CN <- coxph(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + frailty(PlotID, "gaussian"),
                  data = subset(Data_Alt, Species == "CN"))

# View new model results
# High concordance and significant p-value for LRT indicate good model fit
Surv2_CN; summary(Surv2_CN)

# View other model diagnostics
# No clear +/- trend in Schoenfeld residuals for factors, indicating approximately proportional hazard
# While both factor p-values are significant, it's likely just from the high/low blips in last few weeks
# This is probably fine given that hazard is proportional over most of the trimming time
# Global test fails to reject null hypothesis of proportional hazards for model as a whole
plot(cox.zph(Surv2_CN))
cox.zph(Surv2_CN)

# Fit survival model, using initial diameter as covariate [CA]
Surv1_CA <- coxph(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) +
                    frailty(PlotID, "gaussian"),
                  data = subset(Data_Alt, Species == "CA"))

# View model results
Surv1_CA

# Rosette diameter not useful as a covariate; drop and re-fit model
Surv2_CA <- coxph(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + frailty(PlotID, "gaussian"),
                  data = subset(Data_Alt, Species == "CA"))

# View new model results
# High concordance and significant p-value for LRT indicate good model fit
Surv2_CA; summary(Surv2_CA)

# View other model diagnostics
# No clear +/- trend in Schoenfeld residuals for factors, indicating approximately proportional hazard
# While TRT factor p-values are significant, it's likely just from the high/low blips in last few weeks
# This is probably fine given that hazard is proportional over most of the trimming time
# Global test fails to reject null hypothesis of proportional hazards for model as a whole
plot(cox.zph(Surv2_CA))
cox.zph(Surv2_CA)





Surv1_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + factor(PlotID),
                   data = subset(Data_Alt, Species == "CN"))
summary(Surv1_CN)

Surv1_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + factor(PlotID),
                    data = subset(Data_Alt, Species == "CN"))
summary(Surv1_CN)



Surv1_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + factor(PlotID),
                    data = subset(Data_Alt, Species == "CA"))
summary(Surv1_CA)

Surv1_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + factor(PlotID),
                    data = subset(Data_Alt, Species == "CA"))
summary(Surv1_CA)





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

# For CN model max bud count as function of trimming and warming treatments
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

