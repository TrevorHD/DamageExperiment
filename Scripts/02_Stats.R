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

# Generate Kaplan-Meier survival curves for comparison by warming treatment [WIP]
# Note: not using the full time horizon for CN because most of them died before winter
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





##### Fit growth models [WIP] -----------------------------------------------------------------------------

# Note: figure out what to do with the ToD>0 bit
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