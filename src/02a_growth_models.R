library(tidyverse)
library(fuzzyjoin)
library(MASS)
library(lmtest)
library(nlme)
library(stargazer)

load("rad.Rdata")

growth_temp <- growth_com %>%
  select(time_h, flight) %>%
  mutate(time_h = time_h - 4) %>%
  filter(time_h >= 0)

growth_shift <- growth_com %>%
  select(time_h, ground1, ground2, ground3, groundm) %>%
  difference_left_join(., growth_temp, by = "time_h", max_dist = 0.001) %>%
  rename(time_h = "time_h.x") %>%
  select(-time_h.y)

igrowth_com <- dplyr::filter(growth_shift, growth_shift$time_h <= 15) # exponential growth phases (4 - 19h, 0 - 15h)

growth_longm <- dplyr::filter(growth_longm, growth_longm$time_h < 50) # end of growth phase
growth_long1 <- dplyr::filter(growth_long1, growth_long1$time_h < 50)
growth_long2 <- dplyr::filter(growth_long2, growth_long2$time_h < 50)
growth_long3 <- dplyr::filter(growth_long3, growth_long3$time_h < 50)

# Initial [exponential] growth (flight vs. ground control average)

igm_f <- nls(flight ~ I(c * exp(k*time_h)), # Flight
             data = igrowth_com, 
             start = list(c = 0.05, k = 0.15))
sink("igm_f.txt")
summary(igm_f)
sink()

igm_gm <- nls(groundm ~ I(c * exp(k*time_h)), # Ground
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.15))
sink("igm_gm.txt")
summary(igm_gm)
sink()

igrowth_com$preds_f <- predict(igm_f, igrowth_com)
igrowth_com$preds_g <- predict(igm_gm, igrowth_com)

ggplot(igrowth_com) +
  geom_point(aes(x = time_h, y = flight, color = "#F8766D")) +
  geom_point(aes(x = time_h, y = groundm, color = "#00BFC4")) +
  geom_line(aes(x = time_h, y = preds_f, color = "#F8766D")) +
  geom_line(aes(x = time_h, y = preds_g, color = "#00BFC4")) +
  labs(y = "relative OD") +
  scale_x_continuous("Time (in hours, ground)", 
                     sec.axis = sec_axis(~ . + 4, name = "Time (in hours, flight)")) +
  scale_colour_manual(name = "", 
                      values = c("#F8766D" = "#F8766D", "#00BFC4" = "#00BFC4"), 
                      breaks = c("#F8766D", "#00BFC4"),
                      labels = c("Flight", "Ground \n (Average)")) +
  theme(text = element_text(size = 18))

ggsave("gp1.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control average)

jgm_m1 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_longm) # scal = -1/b
summary(jgm_m1)$tTable

jgm_m2 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
               params = list(Asym ~ 1, xmid ~ type, scal ~ type),
               start = list(Asym = c(0.95), xmid = c(15, 15), scal = c(4.5, 4.5)),
               data = growth_longm)
summary(jgm_m2)$tTable
jgm_m2t <- summary(jgm_m2)$tTable

growth_longm$preds <- predict(jgm_m2, data = growth_longm)

ggplot(growth_longm) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (Average)")) +
  theme(text = element_text(size = 18))

ggsave("gp2.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 1)

jgm_g11 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_long1)
summary(jgm_g11)$tTable

jgm_g12 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.95), xmid = c(15, 15), scal = c(4.5, 4.5)),
                data = growth_long1)
summary(jgm_g12)$tTable
jgm_g12t <- summary(jgm_g12)$tTable

growth_long1$preds <- predict(jgm_g12, data = growth_long1)

ggplot(growth_long1) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (Average)")) +
  theme(text = element_text(size = 18))

ggsave("gp3.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 2)

jgm_g21 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long2)
summary(jgm_g21)$tTable

jgm_g22 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.95), xmid = c(15, 15), scal = c(4.5, 4.5)),
                data = growth_long2)
summary(jgm_g22)$tTable
jgm_g22t <- summary(jgm_g22)$tTable

growth_long2$preds <- predict(jgm_g22, data = growth_long2)

ggplot(growth_long2) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (Average)")) +
  theme(text = element_text(size = 18))

ggsave("gp4.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 3)

jgm_g31 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long3)
summary(jgm_g31)$tTable

jgm_g32 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.95), xmid = c(15, 15), scal = c(4.5, 4.5)),
                data = growth_long3)
summary(jgm_g32)$tTable
jgm_g32t <- summary(jgm_g32)$tTable

growth_long3$preds <- predict(jgm_g32, data = growth_long3)

ggplot(growth_long3) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (Average)")) +
  theme(text = element_text(size = 18))

ggsave("gp5.png", width = 9, height = 7)

n <- cbind(nrow(growth_longm), nrow(growth_long1), nrow(growth_long2), nrow(growth_long3))
aic <- round(cbind(AIC(jgm_m2), AIC(jgm_g12), AIC(jgm_g22), AIC(jgm_g32)), 3)
bic <- round(cbind(BIC(jgm_m2), BIC(jgm_g12), BIC(jgm_g22), BIC(jgm_g32)), 3)

stargazer(jgm_m2t, jgm_g12t, jgm_g22t, jgm_g32t, 
          n, aic, bic, 
          out = "jgm.html")
