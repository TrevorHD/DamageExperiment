##### Transform survival data into Kaplan-Meier for survival curves ---------------------------------------

# Generate Kaplan-Meier survival curves for comparison by trimming treatment
# Note: not using the full time horizon for CN because most of them died before winter
km_CN_NW <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Warmed == 0))
km_CN_W <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                   data = subset(Data_Alt_1, Species == "CN" & Warmed == 1))
km_CA_NW <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Warmed == 0))
km_CA_W <- survfit(Surv(ToD, Cens) ~ factor(Treatment), type = "kaplan-meier",
                   data = subset(Data_Alt, Species == "CA" & Warmed == 1))

# Generate Kaplan-Meier survival curves for comparison by warming treatment
# Note: not using the full time horizon for CN because most of them died before winter
km_CN_t1 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Treatment == 1))
km_CN_t2 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Treatment == 2))
km_CN_t3 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Treatment == 3))
km_CN_t4 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt_1, Species == "CN" & Treatment == 4))
km_CA_t1 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Treatment == 1))
km_CA_t2 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Treatment == 2))
km_CA_t3 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Treatment == 3))
km_CA_t4 <- survfit(Surv(ToD, Cens) ~ factor(Warmed), type = "kaplan-meier",
                    data = subset(Data_Alt, Species == "CA" & Treatment == 4))





##### Fit survival regressions ----------------------------------------------------------------------------

# Note: will likely need to add fixed effect instead of random since survreg doesn't support it

# Fit survival model to CN, using initial diameter as covariate
Surv1_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t),
                    data = subset(Data_Alt_1, Species == "CN" & ToD > 0))

# Same as above, but include interaction
Surv2_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment)
                    + factor(Warmed):factor(Treatment) + scale(DM_t),
                    data = subset(Data_Alt_1, Species == "CN" & ToD > 0))

# View model results
summary(Surv1_CN)
summary(Surv2_CN)
AIC(Surv1_CN)
AIC(Surv2_CN)

# In absence of warming, trimming significantly decreased survival time (excl. trt 3)
# For untrimmed individuals, warming significantly decreased survival time (accelerated life cycle?)
# For trimmed individuals, warming generally reduced the negative effects on lifespan from trimming

# Fit survival model to CA, using initial diameter as covariate
Surv1_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + scale(DM_t),
                    data = subset(Data_Alt_2, Species == "CA"))


# Same as above, but include interaction
Surv2_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) +
                      factor(Warmed):factor(Treatment) + scale(DM_t),
                    data = subset(Data_Alt_2, Species == "CA"))

# View model results
# AIC is hardly any lower when adding interaction, adds significant complexity to model
summary(Surv1_CA)
summary(Surv2_CA)
AIC(Surv1_CA)
AIC(Surv2_CA)

# In absence of warming, only the trim to ground significantly decreased survival time
# For untrimmed individuals, warming decreased survival time (accelerated life cycle?) - not significant
# For trimmed individuals, warming generally reduced the negative effects on lifespan from trimming
# Last part above not significant, except for TRT 2

# Scale less than 1 in all models; indicates that risk of death decreases over time





##### Fit regrowth models ---------------------------------------------------------------------------------

# Fit model for time-averaged regrowth (cm) to CN and CA, using initial diameter as covariate
# Then do again with interaction
RegH1_CN <- lmer(HGain ~ factor(Warmed) + factor(Treatment)
                  + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CN"))
RegH2_CN <- lmer(HGain ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                  + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CN"))
RegH1_CA <- lmer(HGain ~ factor(Warmed) + factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CA"))
RegH2_CA <- lmer(HGain ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CA"))

# View CN model results
summary(RegH1_CN)
summary(RegH2_CN)
AIC(RegH1_CN)
AIC(RegH2_CN)

# View CA model results
summary(RegH1_CA)
summary(RegH2_CA)
AIC(RegH1_CA)
AIC(RegH2_CA)

# Fit model for time-averaged stem count change to CN and CA, using initial diameter as covariate
# Then do again with interaction
RegS1_CN <- lmer(SGain ~ factor(Warmed) + factor(Treatment)
                  + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CN"))
RegS2_CN <- lmer(SGain ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                  + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CN"))
RegS1_CA <- lmer(SGain ~ factor(Warmed) + factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CA"))
RegS2_CA <- lmer(SGain ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_TA, Species == "CA"))

# View CN model results
summary(RegS1_CN)
summary(RegS2_CN)
AIC(RegS1_CN)
AIC(RegS2_CN)

# View CA model results
summary(RegS1_CA)
summary(RegS2_CA)
AIC(RegS1_CA)
AIC(RegS2_CA)

# Models generally have a lower AIC with no interaction terms, but not really by much
# As such, we can go ahead and keep interaction terms for reporting
# Plus, they help us address the question of warming effects on response to trimming

# Shapiro-Wilk test rejects normality of residuals in some cases
# However, this might be expected given sensitivity with larger sample sizes
# QQ-Plots are perfectly fine, so likely no concerns here
qqnorm(resid(RegH2_CN))
qqline(resid(RegH2_CN))
shapiro.test(resid(RegH2_CN))
qqnorm(resid(RegH2_CA))
qqline(resid(RegH2_CA))
shapiro.test(resid(RegH2_CA))
qqnorm(resid(RegS2_CN))
qqline(resid(RegS2_CN))
shapiro.test(resid(RegS2_CN))
qqnorm(resid(RegS2_CA))
qqline(resid(RegS2_CA))
shapiro.test(resid(RegS2_CA))





##### Fit reproduction models -----------------------------------------------------------------------------

# Fit model for peak number of reproductive structures for CN, using initial diameter as covariate
# Then do again with interaction
Rep1_CN <- lmer(Total ~ factor(Warmed) + factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_R, Species == "CN"))
Rep2_CN <- lmer(Total ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                 + scale(DM_t) + (1|Clust), data = subset(Data_R, Species == "CN"))

# View model results
summary(Rep1_CN)
summary(Rep2_CN)
AIC(Rep1_CN)
AIC(Rep2_CN)

# Fit model for peak number of reproductive structures for CA, using initial diameter as covariate
# Then do again with interaction
Rep1_CA <- lmer(Total ~ factor(Warmed) + factor(Treatment)
                + scale(DM_t) + (1|Clust), data = subset(Data_R, Species == "CA"))
Rep2_CA <- lmer(Total ~ factor(Warmed) + factor(Treatment) + factor(Warmed):factor(Treatment)
                + scale(DM_t) + (1|Clust), data = subset(Data_R, Species == "CA"))

# View model results
summary(Rep1_CA)
summary(Rep2_CA)
AIC(Rep1_CA)
AIC(Rep2_CA)

# Models generally had a lower AIC with no interaction terms, but not really by much
# As such, we can go ahead and keep interaction terms for reporting
# Plus, they help us address the question of warming effects on response to trimming

# Shapiro-Wilk test rejects normality of residuals
# However, this might be expected given sensitivity with larger sample sizes
# QQ-Plots are perfectly fine, so likely no concerns here
qqnorm(resid(Rep2_CN))
qqline(resid(Rep2_CN))
shapiro.test(resid(Rep2_CN))
qqnorm(resid(Rep2_CA))
qqline(resid(Rep2_CA))
shapiro.test(resid(Rep2_CA))

