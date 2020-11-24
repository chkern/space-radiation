library(tidyverse)
library(MASS)
library(lmtest)
library(sandwich)
library(gvlma)
library(forecast)
library(stargazer)
library(nlme)
library(drc)

load("rad.Rdata")

rad_long_1 <- dplyr::filter(rad_long, rad_long$phases == 1)
rad_long_2 <- dplyr::filter(rad_long, rad_long$phases == 2)
rad_long_3 <- dplyr::filter(rad_long, rad_long$phases == 3)

rad_clean_1 <- dplyr::filter(rad_clean, rad_clean$phases == 1)
rad_clean_2 <- dplyr::filter(rad_clean, rad_clean$phases == 2)
rad_clean_3 <- dplyr::filter(rad_clean, rad_clean$phases == 3)

com_imp_2 <- dplyr::filter(com_imp, com_imp$phases == 2)

# OLS 

m11 <- lm(log(g) ~ type, data = rad_long_1)
m12 <- lm(log(g) ~ time_h + type, data = rad_long_1)
m13 <- lm(log(g) ~ time_h*type, data = rad_long_1)

summary(m11)
summary(m12)
summary(m13)

m21 <- lm(log(g) ~ type, data = rad_long_2)
m22 <- lm(log(g) ~ time_h + type, data = rad_long_2)
m23 <- lm(log(g) ~ time_h*type, data = rad_long_2)

summary(m21)
summary(m22)
summary(m23)

m31 <- lm(log(g) ~ type, data = rad_long_3)
m32 <- lm(log(g) ~ time_h + type, data = rad_long_3)
m33 <- lm(log(g) ~ time_h*type, data = rad_long_3)

summary(m31)
summary(m32)
summary(m33)

plot(m32)
bgtest(m32)
m32gv <- gvlma(m32)
summary(m32gv)

# Wilcoxon tests

stats1 <- rad_clean_1 %>% 
  summarise(mean_exp = mean(exp_g), med_exp = median(exp_g),
            mean_ctrl = mean(ctrl_g), med_ctrl = median(ctrl_g), n())

stats3 <- rad_clean_3 %>% 
  summarise(mean_exp = mean(exp_g), med_exp = median(exp_g),
            mean_ctrl = mean(ctrl_g), med_ctrl = median(ctrl_g), n())

wt1 <- wilcox.test(rad_clean_1$ctrl_g, rad_clean_1$exp_g)
wt2 <- wilcox.test(rad_clean_2$ctrl_g, rad_clean_2$exp_g)
wt3 <- wilcox.test(rad_clean_3$ctrl_g, rad_clean_3$exp_g)

wts <- rbind(wt1$statistic, wt3$statistic)
wtp <- rbind(wt1$p.value, wt3$p.value)
stats <- rbind(stats1, stats3)
tab <- data.frame(stats, wts, wtp)
stargazer(round(tab, 3), summary = FALSE, out = "wilcox.html")

# Robust regressions
## y = g

rm11 <- rlm(g ~ type, data = rad_long_1)
rm12 <- rlm(g ~ time_h + type, data = rad_long_1)
rm13 <- rlm(g ~ time_h*type, data = rad_long_1)

ct11 <- coeftest(rm11, vcov = vcovHC)
ct12 <- coeftest(rm12, vcov = vcovHC)
ct13 <- coeftest(rm13, vcov = vcovHC)

n <- rbind(nrow(rm11$model), nrow(rm12$model), nrow(rm13$model))
aic <- round(rbind(AIC(rm11), AIC(rm12), AIC(rm13)), 3)
bic <- round(rbind(BIC(rm11), BIC(rm12), BIC(rm13)), 3)

stargazer(ct11, ct12, ct13, report = ('vcsp'), 
          add.lines = list(c("n", n), c("AIC", aic), c("BIC", bic)), 
          out = "rlm1.html")

rm21 <- rlm(g ~ type, data = rad_long_2)
rm22 <- rlm(g ~ time_h + type, data = rad_long_2)
rm23 <- rlm(g ~ time_h*type, data = rad_long_2)

coeftest(rm21, vcov = vcovHC)
coeftest(rm22, vcov = vcovHC)
coeftest(rm23, vcov = vcovHC)

rm31 <- rlm(g ~ type, data = rad_long_3)
rm32 <- rlm(g ~ time_h + type, data = rad_long_3)
rm33 <- rlm(g ~ time_h*type, data = rad_long_3)

ct31 <- coeftest(rm31, vcov = vcovHC)
ct32 <- coeftest(rm32, vcov = vcovHC)
ct33 <- coeftest(rm33, vcov = vcovHC)

n <- rbind(nrow(rm31$model), nrow(rm32$model), nrow(rm33$model))
aic <- round(rbind(AIC(rm31), AIC(rm32), AIC(rm33)), 3)
bic <- round(rbind(BIC(rm31), BIC(rm32), BIC(rm33)), 3)

stargazer(ct31, ct32, ct33, report = ('vcsp'), 
          add.lines = list(c("n", n), c("AIC", aic), c("BIC", bic)), 
          out = "rlm3.html")

## y = delta (diff in g)

rmd1 <- rlm(delta ~ relative_OD, data = com_small)
rmd2 <- rlm(delta ~ time_h, data = com_small)
rmd3 <- rlm(delta ~ time_h + relative_OD, data = com_small)

coeftest(rmd1, vcov = vcovHC)
coeftest(rmd2, vcov = vcovHC)
coeftest(rmd3, vcov = vcovHC)

rmd4 <- rlm(delta ~ relative_OD_imp1, data = com_imp)
rmd5 <- rlm(delta ~ time_h, data = com_imp)
rmd6 <- rlm(delta ~ time_h + relative_OD_imp1, data = com_imp)

coeftest(rmd4, vcov = vcovHC)
coeftest(rmd5, vcov = vcovHC)
coeftest(rmd6, vcov = vcovHC)

rmd7 <- rlm(delta ~ relative_OD_imp1, data = com_imp_2)
rmd8 <- rlm(delta ~ time_h, data = com_imp_2)
rmd9 <- rlm(delta ~ time_h + relative_OD_imp1, data = com_imp_2)

ctd7 <- coeftest(rmd7, vcov = vcovHC)
ctd8 <- coeftest(rmd8, vcov = vcovHC)
ctd9 <- coeftest(rmd9, vcov = vcovHC)

n <- rbind(nrow(rmd7$model), nrow(rmd8$model), nrow(rmd9$model))
aic <- round(rbind(AIC(rmd7), AIC(rmd8), AIC(rmd9)), 3)
bic <- round(rbind(BIC(rmd7), BIC(rmd8), BIC(rmd9)), 3)

stargazer(ctd7, ctd8, ctd9, report = ('vcsp'), 
          add.lines = list(c("n", n), c("AIC", aic), c("BIC", bic)), 
          out = "rld2.html")

rel_OD <- data.frame(relative_OD_imp1 = c(0.9, 0.925, 0.95))
predict(rmd7, newdata = rel_OD, interval = "confidence")

mtd1 <- lm(delta_total ~ relative_OD_imp1, data = com_imp)
mtd2 <- lm(delta_total ~ time_h, data = com_imp)
mtd3 <- lm(delta_total ~ time_h + relative_OD_imp1, data = com_imp)

# Compare growth data (flight vs. ground control 6)

gm6_1 <- drm(flight ~ time_h, data = growth_com6, fct = AR.3())
summary(gm6_1) # c = 1/e
plot(gm6_1)

gm6_2 <- drm(ground ~ time_h, data = growth_com6, fct = AR.3())
summary(gm6_2)
plot(gm6_2)

jgm6_1 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model (or: SSasympOrig())
               data = growth_long6) # c = exp(lrc)
summary(jgm6_1)$tTable

jgm6_2 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
               params = Asym + R0 + lrc ~ type,
               start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
               data = growth_long6)
summary(jgm6_2)$tTable

jgm6_3 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
               params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
               start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)),
               data = growth_long6)
summary(jgm6_3)$tTable

growth_long6$preds <- predict(jgm6_3, data = growth_long6)
ggplot(growth_long6) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

# Compare growth data (flight vs. ground control 7)

gm7_1 <- drm(flight ~ time_h, data = growth_com7, fct = AR.3())
summary(gm7_1)
plot(gm7_1)

gm7_2 <- drm(ground ~ time_h, data = growth_com7, fct = AR.3())
summary(gm7_2)
plot(gm7_2)

jgm7_1 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model
               data = growth_long7)
summary(jgm7_1)$tTable

jgm7_2 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
               params = Asym + R0 + lrc ~ type,
               start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
               data = growth_long7)
summary(jgm7_2)$tTable

jgm7_3 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
               params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
               start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)), 
               data = growth_long7)
summary(jgm7_3)$tTable

growth_long7$preds <- predict(jgm7_3, data = growth_long7)
ggplot(growth_long7) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

# Time series

ts_ctrl <- ts(rad_clean$ctrl_g, frequency = 786) # daily seasonality
ts_exp <- ts(rad_clean$exp_g, frequency = 786)

ts_ctrl_3 <- ts(rad_clean_3$ctrl_g, frequency = 786) # daily seasonality
ts_exp_3 <- ts(rad_clean_3$exp_g, frequency = 786)

decomp_ctrl <- decompose(ts_ctrl_3)
plot(decomp_ctrl)
decomp_exp <- decompose(ts_exp_3)
plot(decomp_exp)
