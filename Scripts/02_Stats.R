##### Fit budding probability models ----------------------------------------------------------------------

# Select data for budding probability analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formula for full model
Mod_Form1 <- Budded ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)

# Create formulas for full model, sans interactions, and sans interactions/size
# These will be used for Firth regression since glmer likely won't work
# Note the lack of random effects, as most implementations of Firth don't support them
Mod_Form2 <- Budded ~ factor(Warmed)*factor(Treatment) + scale(DM_t)
Mod_Form3 <- Budded ~ factor(Warmed) + factor(Treatment) + scale(DM_t)
Mod_Form4 <- Budded ~ factor(Warmed) + factor(Treatment)

# For CN, create tables of number of individuals that produced at least one bud
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
table(Data_Temp_CN$Treatment, Data_Temp_CN$Budded)
table(Data_Temp_CN$Warmed, Data_Temp_CN$Budded)
table(Data_Temp_CN$Treatment, Data_Temp_CN$Budded, Data_Temp_CN$Warmed)

# The aforementioned separation causes boundary issues and greatly inflates standard errors
Mod_BudP1_CN <- glmer(Mod_Form1, data = Data_Temp_CN, family = "binomial")
summary(Mod_BudP1_CN)

# Try Firth logistic regression using full model
Mod_BudP2_CN <- logistf(Mod_Form2, data = Data_Temp_CN)
summary(Mod_BudP2_CN)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_BudP3_CN <- logistf(Mod_Form3, data = Data_Temp_CN)
summary(Mod_BudP3_CN)
anova(Mod_BudP2_CN, Mod_BudP3_CN)

# Fit model without size term; size term is not significant, so drop it
Mod_BudP4_CN <- logistf(Mod_Form4, data = Data_Temp_CN)
summary(Mod_BudP4_CN)
anova(Mod_BudP3_CN, Mod_BudP4_CN)

# For CA, create tables of number of individuals that produced at least one bud
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
table(Data_Temp_CA$Treatment, Data_Temp_CA$Budded)
table(Data_Temp_CA$Warmed, Data_Temp_CA$Budded)
table(Data_Temp_CA$Treatment, Data_Temp_CA$Budded, Data_Temp_CA$Warmed)

# The aforementioned separation causes boundary issues and greatly inflates standard errors
Mod_BudP1_CA <- glmer(Mod_Form1, data = Data_Temp_CA, family = "binomial")
summary(Mod_BudP1_CA)

# Try Firth logistic regression using full model
Mod_BudP2_CA <- logistf(Mod_Form2, data = Data_Temp_CA)
summary(Mod_BudP2_CA)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_BudP3_CA <- logistf(Mod_Form3, data = Data_Temp_CA)
summary(Mod_BudP3_CA)
anova(Mod_BudP2_CA, Mod_BudP3_CA)

# Fit model without size term; size term is not significant, so drop it
Mod_BudP4_CA <- logistf(Mod_Form4, data = Data_Temp_CA)
summary(Mod_BudP4_CA)
anova(Mod_BudP3_CA, Mod_BudP4_CA)

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Form4,
       Mod_BudP1_CN, Mod_BudP2_CN, Mod_BudP3_CN, Mod_BudP1_CA, Mod_BudP2_CA, Mod_BudP3_CA)





##### Fit flowering probability models --------------------------------------------------------------------

# Select data for flowering probability analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formula for full model
Mod_Form1 <- Flowered ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)

# Create formulas for full model, sans interactions, and sans interactions/size
# These will be used for Firth regression since glmer likely won't work
# Note the lack of random effects, as most implementations of Firth don't support them
Mod_Form2 <- Flowered ~ factor(Warmed)*factor(Treatment) + scale(DM_t)
Mod_Form3 <- Flowered ~ factor(Warmed) + factor(Treatment) + scale(DM_t)
Mod_Form4 <- Flowered ~ factor(Warmed) + factor(Treatment)

# For CN, create tables of number of individuals that produced at least one flower
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
table(Data_Temp_CN$Treatment, Data_Temp_CN$Flowered)
table(Data_Temp_CN$Warmed, Data_Temp_CN$Flowered)
table(Data_Temp_CN$Treatment, Data_Temp_CN$Flowered, Data_Temp_CN$Warmed)

# The aforementioned separation causes the model to fail entirely
Mod_FloP1_CN <- glmer(Mod_Form1, data = Data_Temp_CN, family = "binomial")
summary(Mod_FloP1_CN)

# Try Firth logistic regression using full model
Mod_FloP2_CN <- logistf(Mod_Form2, data = Data_Temp_CN)
summary(Mod_FloP2_CN)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_FloP3_CN <- logistf(Mod_Form3, data = Data_Temp_CN)
summary(Mod_FloP3_CN)
anova(Mod_FloP2_CN, Mod_FloP3_CN)

# Fit model without size term; size term is not significant, so drop it
Mod_FloP4_CN <- logistf(Mod_Form4, data = Data_Temp_CN)
summary(Mod_FloP4_CN)
anova(Mod_FloP3_CN, Mod_FloP4_CN)

# For CA, create tables of number of individuals that produced at least one flower
# Break data out by trimming, warming, and trimming x warming; trimming shows extreme separation
table(Data_Temp_CA$Treatment, Data_Temp_CA$Flowered)
table(Data_Temp_CA$Warmed, Data_Temp_CA$Flowered)
table(Data_Temp_CA$Treatment, Data_Temp_CA$Flowered, Data_Temp_CA$Warmed)

# The aforementioned separation causes boundary issues and greatly inflates standard errors
Mod_FloP1_CA <- glmer(Mod_Form1, data = Data_Temp_CA, family = "binomial")
summary(Mod_FloP1_CA)

# Try Firth logistic regression using full model
Mod_FloP2_CA <- logistf(Mod_Form2, data = Data_Temp_CA)
summary(Mod_FloP2_CA)

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_FloP3_CA <- logistf(Mod_Form3, data = Data_Temp_CA)
summary(Mod_FloP3_CA)
anova(Mod_FloP2_CA, Mod_FloP3_CA)

# Fit model without size term; size term is not significant, so drop it
Mod_FloP4_CA <- logistf(Mod_Form4, data = Data_Temp_CA)
summary(Mod_FloP4_CA)
anova(Mod_FloP3_CA, Mod_FloP4_CA)

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Form4,
       Mod_FloP1_CN, Mod_FloP2_CN, Mod_FloP3_CN, Mod_FloP1_CA, Mod_FloP2_CA, Mod_FloP3_CA)





##### Fit max bud count models ----------------------------------------------------------------------------

# Select data for max bud count analysis
# Note: data restricted to individuals that budded (>0 buds)
Data_Temp_CN <- subset(Data_Sum, Species == "CN" & MaxBuds > 0)
Data_Temp_CA <- subset(Data_Sum, Species == "CA" & MaxBuds > 0)

# Create formulas for full model, log full model, log sans interactions, and log sans interactions/size
Mod_Form1 <- MaxBuds ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxBuds) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxBuds) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form4 <- log(MaxBuds) ~ factor(Warmed) + factor(Treatment) + (1|PlotID)

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

# Before potentially dropping interactions, test for warming effects within treatments
# No difference between warmed and unwarmed groups within any treatment
pairs(emmeans(Mod_Buds2_CN, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
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

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
pairs(emmeans(Mod_Buds3_CN, "Treatment"))

# Fit model without rosette size; rosette size is significant, so keep it
Mod_Buds4_CN <- lmer(Mod_Form4, data = Data_Temp_CN)
summary(Mod_Buds4_CN)
anova(Mod_Buds3_CN, Mod_Buds4_CN)

# For CA, examine distribution of max bud count (>0 buds) by trimming and warming treatments
# Unsurprisingly, untrimmed plants tend to have more buds than trimmed ones
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

# Before potentially dropping interactions, test for warming effects within treatments
# No difference between warmed and unwarmed groups within any treatment
pairs(emmeans(Mod_Buds2_CA, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
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

# Fit model without rosette size; rosette size is not significant, so drop it
Mod_Buds4_CA <- lmer(Mod_Form4, data = Data_Temp_CA)
summary(Mod_Buds4_CA)
anova(Mod_Buds3_CA, Mod_Buds4_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_Buds4_CA))
qqline(resid(Mod_Buds4_CA))
shapiro.test(resid(Mod_Buds4_CA))
data.frame("Fitted" = fitted(Mod_Buds4_CA), "Resid" = resid(Mod_Buds4_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
pairs(emmeans(Mod_Buds4_CA, "Treatment"))

# Note: all residual plots above will show banding patterns due to discrete nature of data
# Boundary effects also visible on bottom left of plots since log(y) cannot be less than zero
# Neither of these are a cause for concern, nor do they violate model assumptions

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Form4,
       Mod_Buds1_CN, Mod_Buds2_CN, Mod_Buds4_CN, Mod_Buds1_CA, Mod_Buds2_CA, Mod_Buds3_CA)





##### Fit max height models -------------------------------------------------------------------------------

# Select data for max height analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, log full model, log sans interactions, and log sans interactions/size
Mod_Form1 <- MaxH ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxH) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxH) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form4 <- log(MaxH) ~ factor(Warmed) + factor(Treatment) + (1|PlotID)

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

# Before potentially dropping interactions, test for warming effects within treatments
# Significant difference between warmed and unwarmed groups within some treatments
pairs(emmeans(Mod_HMax2_CN, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_HMax3_CN <- lmer(Mod_Form3, data = Data_Temp_CN)
summary(Mod_HMax3_CN)
anova(Mod_HMax2_CN, Mod_HMax3_CN)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_HMax3_CN))
qqline(resid(Mod_HMax3_CN))
shapiro.test(resid(Mod_HMax3_CN))
data.frame("Fitted" = fitted(Mod_HMax3_CN), "Resid" = resid(Mod_HMax3_CN),
           "Treatment" = factor(Data_Temp_CN$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
# Difference between 10 cm and 5 cm trims is significant, though
pairs(emmeans(Mod_HMax3_CN, "Treatment"))

# Fit model without rosette size; rosette size is significant, so keep it
Mod_HMax4_CN <- lmer(Mod_Form4, data = Data_Temp_CN)
summary(Mod_HMax4_CN)
anova(Mod_HMax3_CN, Mod_HMax4_CN)

# For CA, examine distribution of max height by trimming and warming treatments
# Unsurprisingly, untrimmed plants tend to be taller than trimmed ones
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
summary(Mod_HMax2_CA)

# Seems to work; residuals appear normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on tails; not a cause for concern
qqnorm(resid(Mod_HMax2_CA))
qqline(resid(Mod_HMax2_CA))
shapiro.test(resid(Mod_HMax2_CA))
data.frame("Fitted" = fitted(Mod_HMax2_CA), "Resid" = resid(Mod_HMax2_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Before potentially dropping interactions, test for warming effects within treatments
# Significant difference between warmed and unwarmed groups within some treatments
pairs(emmeans(Mod_HMax2_CA, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_HMax3_CA <- lmer(Mod_Form3, data = Data_Temp_CA)
summary(Mod_HMax3_CA)
anova(Mod_HMax2_CA, Mod_HMax3_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on tails; not a cause for concern
qqnorm(resid(Mod_HMax3_CA))
qqline(resid(Mod_HMax3_CA))
shapiro.test(resid(Mod_HMax3_CA))
data.frame("Fitted" = fitted(Mod_HMax3_CA), "Resid" = resid(Mod_HMax3_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Fit model without rosette size; rosette size is not significant, so drop it
Mod_HMax4_CA <- lmer(Mod_Form4, data = Data_Temp_CA)
summary(Mod_HMax4_CA)
anova(Mod_HMax3_CA, Mod_HMax4_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
# Shapiro-Wilk still rejects normality on tails; not a cause for concern
qqnorm(resid(Mod_HMax4_CA))
qqline(resid(Mod_HMax4_CA))
shapiro.test(resid(Mod_HMax4_CA))
data.frame("Fitted" = fitted(Mod_HMax4_CA), "Resid" = resid(Mod_HMax4_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
# Difference between 10 cm and 0 cm trims is significant, though
pairs(emmeans(Mod_HMax4_CA, "Treatment"))

# Remove temporary variables that are no longer needed
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Form4,
       Mod_HMax1_CN, Mod_HMax2_CN, Mod_HMax4_CN, Mod_HMax1_CA, Mod_HMax2_CA, Mod_HMax3_CA)





##### Fit max stem count models ---------------------------------------------------------------------------

# Select data for max stem count analysis
Data_Temp_CN <- subset(Data_Sum, Species == "CN")
Data_Temp_CA <- subset(Data_Sum, Species == "CA")

# Create formulas for full model, log full model, log sans interactions, and log sans interactions/size
Mod_Form1 <- MaxS ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form2 <- log(MaxS) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form3 <- log(MaxS) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID)
Mod_Form4 <- log(MaxS) ~ factor(Warmed) + factor(Treatment) + (1|PlotID)

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

# Before potentially dropping interactions, test for warming effects within treatments
# No difference between warmed and unwarmed groups within any treatment
pairs(emmeans(Mod_SMax2_CN, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
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

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
pairs(emmeans(Mod_SMax3_CN, "Treatment"))

# Fit model without rosette size; rosette size is significant, so keep it
Mod_SMax4_CN <- lmer(Mod_Form4, data = Data_Temp_CN)
summary(Mod_SMax4_CN)
anova(Mod_SMax3_CN, Mod_SMax4_CN)

# For CA, examine distribution of max stem count by trimming and warming treatments
# Trimming appears to disrupt apical dominance and results in more stems
# No obvious differences for warmed vs unwarmed plants, though
dd.plot(Data_Temp_CN, "MaxS", "Treatment", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Treatment", "dens", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Warmed", "hist", c(0, 15))
dd.plot(Data_Temp_CN, "MaxS", "Warmed", "dens", c(0, 15))

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

# Before potentially dropping interactions, test for warming effects within treatments
# No difference between warmed and unwarmed groups within any treatment
pairs(emmeans(Mod_SMax2_CA, ~ Warmed | Treatment))

# Fit model without interaction term; interaction term is not significant, so drop it
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

# Fit model without rosette size; rosette size is not significant, so drop it
Mod_SMax4_CA <- lmer(Mod_Form4, data = Data_Temp_CA)
summary(Mod_SMax4_CA)
anova(Mod_SMax3_CA, Mod_SMax4_CA)

# Diagnostics still fine; residuals normal with little to no heteroskedasticity
qqnorm(resid(Mod_SMax4_CA))
qqline(resid(Mod_SMax4_CA))
shapiro.test(resid(Mod_SMax4_CA))
data.frame("Fitted" = fitted(Mod_SMax4_CA), "Resid" = resid(Mod_SMax4_CA),
           "Treatment" = factor(Data_Temp_CA$Treatment)) %>%
  ggplot(aes(x = Fitted, y = Resid, colour = Treatment)) +
  geom_point()

# Test for any significant differences between trimming treatments
# Trimming treatments significantly different from control, but not from each other
# Difference between 5 cm and 0 cm trims is significant, though
pairs(emmeans(Mod_SMax4_CA, "Treatment"))

# Note: all residual plots above will show banding patterns due to discrete nature of data
# Boundary effects also visible on bottom left of plots since log(y) cannot be less than zero
# Neither of these are a cause for concern, nor do they violate model assumptions

# Remove intermediate models that are no longer used
remove(Data_Temp_CN, Data_Temp_CA, Mod_Form1, Mod_Form2, Mod_Form3, Mod_Form4,
       Mod_SMax1_CN, Mod_SMax2_CN, Mod_SMax4_CN, Mod_SMax1_CA, Mod_SMax2_CA, Mod_SMax3_CA)





##### Fit survival regressions ----------------------------------------------------------------------------

# Select data for survival analysis
Data_Temp_CN_Y1 <- subset(Data_Alt_Y1, Species == "CN")
Data_Temp_CA_Y1 <- subset(Data_Alt_Y1, Species == "CA")
Data_Temp_CN_Y2 <- subset(Data_Alt_Y2, Species == "CN")
Data_Temp_CA_Y2 <- subset(Data_Alt_Y2, Species == "CA")

# Note: we can't just list out the formulas here then use them later like we did above
# For some reason, it throws an error when using ANOVA to compare nested models
# Thus, we'll specify the formula in the same code that actually fits the model

# For CN (Year 1), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CN_Y1$Treatment, Data_Temp_CN_Y1$Cens)
table(Data_Temp_CN_Y1$Warmed, Data_Temp_CN_Y1$Cens)
table(Data_Temp_CN_Y1$Treatment, Data_Temp_CN_Y1$Cens, Data_Temp_CN_Y1$Warmed)

# Model hazard as function of trimming and warming treatments
# 2 out of 60 individuals survived until end of trimming season
Mod_Surv1_CN <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CN"))
summary(Mod_Surv1_CN)

# Get hazard ratios and calculate 95% CI for full model
summary(Mod_Surv1_CN)$coef[, 2]
exp(summary(Mod_Surv1_CN)$coef[, 1] - summary(Mod_Surv1_CN)$coef[, 3])
exp(summary(Mod_Surv1_CN)$coef[, 1] + summary(Mod_Surv1_CN)$coef[, 3])

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_Surv2_CN <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CN"))
summary(Mod_Surv2_CN)
anova(Mod_Surv1_CN, Mod_Surv2_CN)

# Fit model without size term; size term is not significant, so drop it
Mod_Surv3_CN <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CN"))
summary(Mod_Surv3_CN)
anova(Mod_Surv2_CN, Mod_Surv3_CN)

# Get hazard ratios and calculate 95% CI for reduced model
summary(Mod_Surv3_CN)$coef[, 2]
exp(summary(Mod_Surv3_CN)$coef[, 1] - summary(Mod_Surv3_CN)$coef[, 3])
exp(summary(Mod_Surv3_CN)$coef[, 1] + summary(Mod_Surv3_CN)$coef[, 3])

# Examine diagnostics of reduced model
# Note: hazard for treatment is mostly constant, but decreases near the end
# This is what is likely driving a significant value here for treatment, and is not a cause for concern
cox.zph(Mod_Surv3_CN)
plot(cox.zph(Mod_Surv3_CN))

# For CN (Year 2), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CN_Y2$Treatment, Data_Temp_CN_Y2$Cens)
table(Data_Temp_CN_Y2$Warmed, Data_Temp_CN_Y2$Cens)
table(Data_Temp_CN_Y2$Treatment, Data_Temp_CN_Y2$Cens, Data_Temp_CN_Y2$Warmed)

# We will not fit a survival model for CN in year 2
# Only two individuals remain (not enough to fit a model), and both survive until the end anyways

# For CA (Year 1), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CA_Y1$Treatment, Data_Temp_CA_Y1$Cens)
table(Data_Temp_CA_Y1$Warmed, Data_Temp_CA_Y1$Cens)
table(Data_Temp_CA_Y1$Treatment, Data_Temp_CA_Y1$Cens, Data_Temp_CA_Y1$Warmed)

# Model hazard as function of trimming and warming treatments
# 23 out of 60 individuals survived until end of trimming season
Mod_Surv1_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CA"))
summary(Mod_Surv1_CA)

# Get hazard ratios and calculate 95% CI for full model
summary(Mod_Surv1_CA)$coef[, 2]
exp(summary(Mod_Surv1_CA)$coef[, 1] - summary(Mod_Surv1_CA)$coef[, 3])
exp(summary(Mod_Surv1_CA)$coef[, 1] + summary(Mod_Surv1_CA)$coef[, 3])

# Fit model without interaction term; interaction term is significant, so keep it
Mod_Surv2_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CA"))
summary(Mod_Surv2_CA)
anova(Mod_Surv1_CA, Mod_Surv2_CA)

# Fit model without size term; size term is not significant, so drop it
Mod_Surv3_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y1, Species == "CA"))
summary(Mod_Surv3_CA)
anova(Mod_Surv1_CA, Mod_Surv3_CA)

# Get hazard ratios and calculate 95% CI for reduced model
summary(Mod_Surv3_CA)$coef[, 2]
exp(summary(Mod_Surv3_CA)$coef[, 1] - summary(Mod_Surv3_CA)$coef[, 3])
exp(summary(Mod_Surv3_CA)$coef[, 1] + summary(Mod_Surv3_CA)$coef[, 3])

# Examine diagnostics of reduced model; looks good
cox.zph(Mod_Surv3_CA)
plot(cox.zph(Mod_Surv3_CA))

# For CA (Year 2), create tables of number of individuals that survived
# Break data out by trimming, warming, and trimming x warming
table(Data_Temp_CA_Y2$Treatment, Data_Temp_CA_Y2$Cens)
table(Data_Temp_CA_Y2$Warmed, Data_Temp_CA_Y2$Cens)
table(Data_Temp_CA_Y2$Treatment, Data_Temp_CA_Y2$Cens, Data_Temp_CA_Y2$Warmed)

# Model hazard as function of trimming and warming treatments
# 10 out of 19 individuals survived until end of trimming season
Mod_Surv4_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed)*factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y2, Species == "CA"))
summary(Mod_Surv4_CA)

# Get hazard ratios and calculate 95% CI for full model
summary(Mod_Surv4_CA)$coef[, 2]
exp(summary(Mod_Surv4_CA)$coef[, 1] - summary(Mod_Surv4_CA)$coef[, 3])
exp(summary(Mod_Surv4_CA)$coef[, 1] + summary(Mod_Surv4_CA)$coef[, 3])

# Fit model without interaction term; interaction term is not significant, so drop it
Mod_Surv5_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t) + (1|PlotID),
                      data = subset(Data_Alt_Y2, Species == "CA"))
summary(Mod_Surv5_CA)
anova(Mod_Surv4_CA, Mod_Surv5_CA)

# Fit model without size term; size term is not significant, so drop it
Mod_Surv6_CA <- coxme(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + (1|PlotID),
                      data = subset(Data_Alt_Y2, Species == "CA"))
summary(Mod_Surv6_CA)
anova(Mod_Surv5_CA, Mod_Surv6_CA)

# Get hazard ratios and calculate 95% CI for reduced model
summary(Mod_Surv6_CA)$coef[, 2]
exp(summary(Mod_Surv6_CA)$coef[, 1] - summary(Mod_Surv6_CA)$coef[, 3])
exp(summary(Mod_Surv6_CA)$coef[, 1] + summary(Mod_Surv6_CA)$coef[, 3])

# Examine diagnostics of reduced model; looks good
cox.zph(Mod_Surv2_CA)
plot(cox.zph(Mod_Surv2_CA))

# Remove intermediate models that are no longer used
remove(Data_Temp_CN_Y1, Data_Temp_CN_Y2, Data_Temp_CA_Y1, Data_Temp_CA_Y2,
       Mod_Surv1_CN, Mod_Surv2_CN, Mod_Surv1_CA, Mod_Surv2_CA, Mod_Surv4_CA, Mod_Surv5_CA)

