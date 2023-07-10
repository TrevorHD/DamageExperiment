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

# Note: figure out what to do with the ToD>0 bit
Surv1_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CN" & ToD > 0))

# Note: figure out what to do with the ToD>0 bit
Surv2_CN <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment)
                    + factor(Warmed):factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CN" & ToD > 0))

# View model results
summary(Surv1_CN)
summary(Surv2_CN)
AIC(Surv1_CN)
AIC(Surv2_CN)

# In absence of warming, trimming significantly decreased survival time (excl. trt 3)
# For untrimmed individuals, warming significantly decreased survival time (accelerated life cycle?)
# For trimmed individuals, warming generally reduced the negative effects on lifespan from trimming

# Fit survival model to CA, using initial diameter as covariate
Surv1_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CA"))


# Same as above, but include interaction
Surv2_CA <- survreg(Surv(ToD, Cens) ~ factor(Warmed) + factor(Treatment) +
                      factor(Warmed):factor(Treatment) + DM_t,
                    data = subset(Data_Alt, Species == "CA"))

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





##### Fit growth models [WIP] -----------------------------------------------------------------------------

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