library(tidyverse)
library(fuzzyjoin)
library(MASS)
library(lmtest)
library(nlme)
library(stargazer)

load("rad.Rdata")

growth_tempg2 <- growth_com %>%
  select(time_h, ground2) %>%
  mutate(time_h = time_h - 5.5) %>%
  filter(time_h >= 0)

growth_tempg3 <- growth_com %>%
  select(time_h, ground3) %>%
  mutate(time_h = time_h - 2.5) %>%
  filter(time_h >= 0)

growth_tempf <- growth_com %>%
  select(time_h, flight) %>%
  mutate(time_h = time_h - 5) %>%
  filter(time_h >= 0)

growth_shift <- growth_com %>%
  select(time_h, ground1, groundm) %>%
  difference_left_join(., growth_tempg2, by = "time_h", max_dist = 0.01) %>%
  rename(time_h = "time_h.x") %>%
  select(-time_h.y)

growth_shift <- growth_shift %>%
  select(time_h, ground1, ground2, groundm) %>%
  difference_left_join(., growth_tempg3, by = "time_h", max_dist = 0.01) %>%
  rename(time_h = "time_h.x") %>%
  select(-time_h.y)

growth_shift <- growth_shift %>%
  select(time_h, ground1, ground2, ground3, groundm) %>%
  difference_left_join(., growth_tempf, by = "time_h", max_dist = 0.01) %>%
  rename(time_h = "time_h.x") %>%
  select(-time_h.y)

# growth_shift$means <- rowMeans(growth_shift[,2:4])

igrowth_com <- dplyr::filter(growth_shift, growth_shift$time_h <= 10) # exponential growth phases (5 - 15h flight)

igrowth_long <- igrowth_com %>% # exponential growth phases (5 - 15h flight)
  select(time_h, ground1, ground2, ground3) %>%
  pivot_longer(cols = ground1:ground3, names_to = c("type")) %>%
  drop_na()

growth_longm <- growth_shift %>% # growth phase (5 - 50h flight)
  select(time_h, flight, groundm) %>%
  filter(time_h < 45) %>%
  pivot_longer(cols = flight:groundm, names_to = c("type")) %>%
  drop_na()

growth_long1 <- growth_shift %>% # growth phase (5 - 50h flight, 0 - 45h ground1)
  select(time_h, flight, ground1) %>%
  filter(time_h < 45) %>%
  pivot_longer(cols = flight:ground1, names_to = c("type")) %>%
  drop_na()

growth_long2 <- growth_shift %>% # growth phase (5 - 50h flight, 5.5 - 50.5h ground2)
  select(time_h, flight, ground2) %>%
  filter(time_h < 45) %>%
  pivot_longer(cols = flight:ground2, names_to = c("type")) %>%
  drop_na()

growth_long3 <- growth_shift %>% # growth phase (5 - 50h flight, 2.5 - 47.5 ground3)
  select(time_h, flight, ground3) %>%
  filter(time_h < 45) %>%
  pivot_longer(cols = flight:ground3, names_to = c("type")) %>%
  drop_na()

# Initial [exponential / power] growth (flight vs. ground control average)

igm_f <- nls(flight ~ I(c * exp(k*time_h)), # Flight exponential
             data = igrowth_com, 
             start = list(c = 0.05, k = 0.1))
sink("igm_f.txt")
summary(igm_f)
sink()

igm_fp <- nls(flight ~ a * time_h^b, # Flight power
             data = igrowth_com, 
             start = list(a = 0.5, b = 2))
sink("igm_fp.txt")
summary(igm_fp)
sink()

igm_gm <- nls(groundm ~ I(c * exp(k*time_h)), # Ground exponential
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.1))
sink("igm_gm.txt")
summary(igm_gm)
sink()

igm_g1 <- nls(ground1 ~ I(c * exp(k*time_h)), # Ground 1 exponential
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.1))
sink("igm_g1.txt")
summary(igm_g1)
sink()

igm_g2 <- nls(ground2 ~ I(c * exp(k*time_h)), # Ground 2 exponential
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.1))
sink("igm_g2.txt")
summary(igm_g2)
sink()

igm_g3 <- nls(ground3 ~ I(c * exp(k*time_h)), # Ground 3 exponential
              data = igrowth_com, 
              start = list(c = 0.05, k = 0.1))
sink("igm_g3.txt")
summary(igm_g3)
sink()

igm_gmp <- nls(groundm ~ a * time_h^b, # Ground power
              data = igrowth_com, 
              start = list(a = 0.5, b = 2))
sink("igm_gmp.txt")
summary(igm_gmp)
sink()

igm_gmc <- nls(value ~ I(c * exp(k*time_h)), # Ground combined exponential
              data = igrowth_long, 
              start = list(c = 0.05, k = 0.1))
sink("igm_gmc.txt")
summary(igm_gmc)
confint(igm_gmc)
sink()

igm_gmcp <- nls(value ~ a * time_h^b, # Ground combined power
               data = igrowth_long, 
               start = list(a = 0.5, b = 2))
sink("igm_gmcp.txt")
summary(igm_gmcp)
confint(igm_gmcp)
sink()

igrowth_com$preds_f <- predict(igm_f, igrowth_com)
igrowth_com$preds_g <- predict(igm_gm, igrowth_com)
igrowth_com$preds_gc <- predict(igm_gmc, igrowth_com)

ggplot(igrowth_com) +
  geom_point(aes(x = time_h, y = flight, color = "#F8766D")) +
  geom_point(aes(x = time_h, y = groundm, color = "#00BFC4")) +
  geom_line(aes(x = time_h, y = preds_f, color = "#F8766D")) +
  geom_line(aes(x = time_h, y = preds_gc, color = "#00BFC4")) +
  labs(y = "relative OD") +
  scale_x_continuous("Time (in hours, flight)", 
                     labels = c("5", "7.5", "10", "12.5", "15")) +
  scale_colour_manual(name = "", 
                      values = c("#F8766D" = "#F8766D", "#00BFC4" = "#00BFC4"), 
                      breaks = c("#F8766D", "#00BFC4"),
                      labels = c("Flight", "Ground \n (Average)")) +
  annotate(geom = "text", x = 2.5, y = 0.375, label="y = 0.017*exp(0.299*hour)", color="#F8766D", size = 6) +
  annotate(geom = "text", x = 2.5, y = 0.35, label="y = 0.033*exp(0.226*hour)", color="#00BFC4", size = 6) +
  theme(text = element_text(size = 17))

ggsave("gp1a.png", width = 9, height = 7)

igrowth_com$preds_fp <- predict(igm_fp, igrowth_com)
igrowth_com$preds_gp <- predict(igm_gmp, igrowth_com)
igrowth_com$preds_gpc <- predict(igm_gmcp, igrowth_com)

ggplot(igrowth_com) +
  geom_point(aes(x = time_h, y = flight, color = "#F8766D")) +
  geom_point(aes(x = time_h, y = groundm, color = "#00BFC4")) +
  geom_line(aes(x = time_h, y = preds_fp, color = "#F8766D")) +
  geom_line(aes(x = time_h, y = preds_gpc, color = "#00BFC4")) +
  labs(y = "relative OD") +
  scale_x_continuous("Time (in hours, flight)", 
                     labels = c("5", "7.5", "10", "12.5", "15")) +
  scale_colour_manual(name = "", 
                      values = c("#F8766D" = "#F8766D", "#00BFC4" = "#00BFC4"), 
                      breaks = c("#F8766D", "#00BFC4"),
                      labels = c("Flight", "Ground \n (Average)")) +
  annotate(geom = "text", x = 2, y = 0.375, label="y == 0.003*hour^2.004", color="#F8766D", size = 6, parse = T) +
  annotate(geom = "text", x = 2, y = 0.35, label="y == 0.013*hour^1.333", color="#00BFC4", size = 6, parse = T) +
  theme(text = element_text(size = 17))

ggsave("gp1b.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control average)

jgm_m1 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_longm) # scal = -1/b
summary(jgm_m1)$tTable

jgm_m2 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
               params = list(Asym ~ 1, xmid ~ type, scal ~ type),
               start = list(Asym = c(0.9), xmid = c(14, 14), scal = c(4.5, 4.5)),
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
  scale_x_continuous("Time (in hours, flight)", 
                     breaks = c(0, 10, 20, 30, 40),
                     labels = c("5", "15", "25", "35", "45")) +
  theme(text = element_text(size = 20))

ggsave("gp2.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 1)

jgm_g11 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
               data = growth_long1)
summary(jgm_g11)$tTable

jgm_g12 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14, 14), scal = c(4.5, 4.5)),
                data = growth_long1)
summary(jgm_g12)$tTable
jgm_g12t <- summary(jgm_g12)$tTable

growth_long1$preds <- predict(jgm_g12, data = growth_long1)

ggplot(growth_long1) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (1)")) +
  scale_x_continuous("Time (in hours, flight)", 
                     breaks = c(0, 10, 20, 30, 40),
                     labels = c("5", "15", "25", "35", "45")) +
  theme(text = element_text(size = 20))

ggsave("gp3.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 2)

jgm_g21 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long2)
summary(jgm_g21)$tTable

jgm_g22 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14, 14), scal = c(4.5, 4.5)),
                data = growth_long2)
summary(jgm_g22)$tTable
jgm_g22t <- summary(jgm_g22)$tTable

growth_long2$preds <- predict(jgm_g22, data = growth_long2)

ggplot(growth_long2) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (2)")) +
  scale_x_continuous("Time (in hours, flight)", 
                     breaks = c(0, 10, 20, 30, 40),
                     labels = c("5", "15", "25", "35", "45")) +
  theme(text = element_text(size = 20))

ggsave("gp4.png", width = 9, height = 7)

# Full [logistic] growth (flight vs. ground control 3)

jgm_g31 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Joint model
                data = growth_long3)
summary(jgm_g31)$tTable

jgm_g32 <- gnls(value ~ SSlogis(time_h, Asym, xmid, scal), # Test for differences in slope
                params = list(Asym ~ 1, xmid ~ type, scal ~ type),
                start = list(Asym = c(0.9), xmid = c(14, 14), scal = c(4.5, 4.5)),
                data = growth_long3)
summary(jgm_g32)$tTable
jgm_g32t <- summary(jgm_g32)$tTable

growth_long3$preds <- predict(jgm_g32, data = growth_long3)

ggplot(growth_long3) +
  geom_point(aes(x = time_h, y = value, color = type)) +
  geom_line(aes(x = time_h, y = preds, color = type)) +
  labs(y = "relative OD" , x = "Time (in hours)") +
  scale_color_discrete(name = "",
                       labels = c("Flight", "Ground \n (3)")) +
  scale_x_continuous("Time (in hours, flight)", 
                     breaks = c(0, 10, 20, 30, 40),
                     labels = c("5", "15", "25", "35", "45")) +
  theme(text = element_text(size = 20))

ggsave("gp5.png", width = 9, height = 7)

n <- cbind(nrow(growth_longm), nrow(growth_long1), nrow(growth_long2), nrow(growth_long3))
aic <- round(cbind(AIC(jgm_m2), AIC(jgm_g12), AIC(jgm_g22), AIC(jgm_g32)), 3)
bic <- round(cbind(BIC(jgm_m2), BIC(jgm_g12), BIC(jgm_g22), BIC(jgm_g32)), 3)

stargazer(jgm_m2t, jgm_g12t, jgm_g22t, jgm_g32t, 
          n, aic, bic, 
          out = "jgm.html")
