library(tidyverse)
library(MASS)
library(lmtest)
library(nlme)
library(stargazer)

load("rad.Rdata")

igrowth_com <-dplyr::filter(growth_com, growth_com$time_h <= 15) # end of initial growth phase
igrowth_longm <- dplyr::filter(growth_longm, growth_longm$time_h <= 15)

growth_longm <- dplyr::filter(growth_longm, growth_longm$time_h < 50 - 4) # end of growth phase
growth_long1 <- dplyr::filter(growth_long1, growth_long1$time_h < 50 - 4)
growth_long2 <- dplyr::filter(growth_long2, growth_long2$time_h < 50 - 4)
growth_long3 <- dplyr::filter(growth_long3, growth_long3$time_h < 50 - 4)

# Initial [exponential] growth (flight vs. ground control average)

igm_f <- nls(flight ~ I(c * exp(k*time_h)), # Flight
             data = igrowth_com, 
             start = list(c = 0.05, k = 0.15))
summary(igm_f)

igm_gm <- nls(groundm ~ I(c * exp(k*time_h)), # Ground
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.15))
summary(igm_gm)

igrowth_com$preds_f <- predict(igm_f, igrowth_com)
igrowth_com$preds_g <- predict(igm_gm, igrowth_com)

ggplot(igrowth_com) +
  geom_point(aes(x = time_h, y = flight), color = "#F8766D") +
  geom_point(aes(x = time_h, y = groundm), color = "#00BFC4") +
  geom_line(aes(x = time_h, y = preds_f), color = "#F8766D") +
  geom_line(aes(x = time_h, y = preds_g), color = "#00BFC4") +
  labs(y = "relative OD" , x = "Time (in hours)")

igm_m1 <- gnls(value ~ c * exp(k*time_h), # Joint model
               data = igrowth_longm,
               start = list(c = 0.05, k = 0.15),
               control = gnlsControl(nlsTol = 0.1)) 
summary(igm_m1)$tTable

igm_m2 <- gnls(value ~ c * exp(k*time_h), # Test for differences in slope
               params = list(c ~ 1, k ~ type),
               start = list(c = 0.05, k = c(0.15, 0.15)),
               data = igrowth_longm,
               control = gnlsControl(nlsTol = 0.1))
summary(igm_m2)$tTable

igrowth_longm$preds <- predict(igm_m2, data = igrowth_longm)

ggplot(igrowth_longm) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  theme(legend.title = element_blank()) +
  labs(y = "relative OD" , x = "Time (in hours)")

# Full [logistic] growth (flight vs. ground control average)

jgm_m1 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_longm) # scal = -1/b
summary(jgm_m1)$tTable

jgm_m2 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
               params = list(Asym ~ 1, xmid ~ 1, scal ~ type),
               start = list(Asym = c(0.9), xmid = c(14), scal = c(4.5, 4.5)),
               data = growth_longm)
summary(jgm_m2)$tTable
jgm_m2t <- summary(jgm_m2)$tTable

growth_longm$preds <- predict(jgm_m2, data = growth_longm)

ggplot(growth_longm) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  theme(legend.title = element_blank()) +
  labs(y = "relative OD" , x = "Time (in hours)")

# Full [logistic] growth (flight vs. ground control 1)

jgm_g11 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_long1)
summary(jgm_g11)$tTable

jgm_g12 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ 1, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14), scal = c(4.5, 4.5)),
                data = growth_long1)
summary(jgm_g12)$tTable
jgm_g12t <- summary(jgm_g12)$tTable

# Full [logistic] growth (flight vs. ground control 2)

jgm_g21 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long2)
summary(jgm_g21)$tTable

jgm_g22 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ 1, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14), scal = c(4.5, 4.5)),
                data = growth_long2)
summary(jgm_g22)$tTable
jgm_g22t <- summary(jgm_g22)$tTable

# Full [logistic] growth (flight vs. ground control 3)

jgm_g31 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long3)
summary(jgm_g31)$tTable

jgm_g32 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ 1, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14), scal = c(4.5, 4.5)),
                data = growth_long3)
summary(jgm_g32)$tTable
jgm_g32t <- summary(jgm_g32)$tTable

n <- cbind(nrow(growth_longm), nrow(growth_long1), nrow(growth_long2), nrow(growth_long3))
aic <- round(cbind(AIC(jgm_m2), AIC(jgm_g12), AIC(jgm_g22), AIC(jgm_g32)), 3)
bic <- round(cbind(BIC(jgm_m2), BIC(jgm_g12), BIC(jgm_g22), BIC(jgm_g32)), 3)

stargazer(jgm_m2t, jgm_g12t, jgm_g22t, jgm_g32t, 
          n, aic, bic, 
          out = "jgm.html")
