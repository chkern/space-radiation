library(tidyverse)
library(MASS)
library(lmtest)
library(stargazer)
library(nlme)
library(drc)

load("rad.Rdata")

# Compare growth data (flight vs. ground control average)

gm_f <- drm(flight ~ time_h, data = growth_com, fct = AR.3())
summary(gm_f) # c = 1/e
plot(gm_f)

gm_gm <- drm(groundm ~ time_h, data = growth_com, fct = AR.3())
summary(gm_gm)
plot(gm_gm)

jgm_m1 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model (or: SSasympOrig())
               data = growth_longm) # c = exp(lrc)
summary(jgm_m1)$tTable

jgm_m2 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
               params = Asym + R0 + lrc ~ type,
               start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
               data = growth_longm)
summary(jgm_m2)$tTable

jgm_m3 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
               params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
               start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)),
               data = growth_longm)
summary(jgm_m3)$tTable

growth_longm$preds <- predict(jgm_m3, data = growth_longm)
ggplot(growth_longm) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

# Compare growth data (flight vs. ground control 1)

jgm_g11 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model
               data = growth_long1)
summary(jgm_g11)$tTable

jgm_g12 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
               params = Asym + R0 + lrc ~ type,
               start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
               data = growth_long1)
summary(jgm_g12)$tTable

jgm_g13 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
               params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
               start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)), 
               data = growth_long1)
summary(jgm_g13)$tTable

# Compare growth data (flight vs. ground control 2)

jgm_g21 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model
                data = growth_long2)
summary(jgm_g21)$tTable

jgm_g22 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
                params = Asym + R0 + lrc ~ type,
                start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
                data = growth_long2)
summary(jgm_g22)$tTable

jgm_g23 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
                params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
                start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)), 
                data = growth_long2)
summary(jgm_g23)$tTable

# Compare growth data (flight vs. ground control 3)

jgm_g31 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Joint model
                data = growth_long3)
summary(jgm_g31)$tTable

jgm_g32 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences
                params = Asym + R0 + lrc ~ type,
                start = list(Asym = c(1, 1), R0 = c(0, 0), lrc = c(-3, -3)),
                data = growth_long3)
summary(jgm_g32)$tTable

jgm_g33 <- gnls(value ~ SSasymp(time_h, Asym, R0, lrc), # Test for differences in slope
                params = list(Asym ~ 1, R0 ~ 1, lrc ~ type),
                start = list(Asym = c(1), R0 = c(0), lrc = c(-3, -3)), 
                data = growth_long3)
summary(jgm_g33)$tTable
