library(tidyverse)
library(readxl)
library(kader)
library(ggpmisc)
library(GGally)
library(fuzzyjoin)
library(imputeTS)
library(data.table)
library(WRS2)

# Read data

rad_raw <- read_excel("Supplementary File 2.xlsx", 
                      sheet = "rad")

growth_flight_raw <- read_excel("Supplementary File 2.xlsx",
                                sheet = "growth",
                                range = cell_cols("K:N"))

growth_ground_raw <- read_excel("Supplementary File 2.xlsx",
                               sheet = "growth",
                               range = cell_cols("A:J"))

summary(rad_raw)
summary(growth_flight_raw)
summary(growth_ground_raw)

# Clean

rad_clean <- rad_raw %>% 
  slice(1:(n()-4)) %>%
  rename(time_sec = "time [sec]",
         time_min = "time [min]",
         time_h = "time [h]",
         time_days = "time [days]",
         ctrl_g = "g_cont",
         ctrl_gn = "g_cont_n",
         ctrl_gt = "g_cont_t",
         ctrl_gnt = "g_cont_nt",
         ctrl_totaln = "ctrl_total_n",
         exp_g = "g_exp",
         exp_gn = "g_exp_n",
         exp_gt = "g_exp_t",
         exp_gnt = "g_exp_nt",
         exp_totaln = "exp_total_n") %>%
  mutate(delta_cr = kader:::cuberoot(delta))

growth_flight_clean <- growth_flight_raw %>% 
  slice(2:(n())) %>%
  rename(time_h = "...2",
         flight = "...4") %>% 
  select(time_h, flight) %>%
  mutate(time_h = as.numeric(time_h),
         flight = as.numeric(flight))

min(growth_flight_clean$time_h[growth_flight_clean$flight >= 0.05])
min(growth_flight_clean$time_h[growth_flight_clean$flight >= 0.95])

rad_clean <- rad_clean %>% 
  mutate(phase1 = ifelse(time_h <= 10, 1, 0),
         phase3 = ifelse(time_h >= 50, 1, 0),
         phases = case_when(phase1 == 1 ~ 1,
                            phase1 == 0 & phase3 == 0 ~ 2,
                            phase3 == 1 ~ 3))

growth_ground_clean <- growth_ground_raw %>% 
  slice(2:(n())) %>%
  rename(time_h = "...2",
         ground1 = "...6",
         ground2 = "...7",
         ground3 = "...8",
         groundm = "...9") %>% 
  select(time_h, 
         ground1, ground2, ground3, groundm) %>%
  mutate(time_h = as.numeric(time_h),
         ground1 = as.numeric(ground1),
         ground2 = as.numeric(ground2),
         ground3 = as.numeric(ground3),
         groundm = as.numeric(groundm))

# Join

com_imp <- rad_clean %>%
  difference_left_join(., growth_flight_clean, by = "time_h", max_dist = 0.05) %>% # Fuzzy match by hour
  rename(time_h = "time_h.x", 
         relative_OD = "flight") %>%
  mutate(dub = duplicated(.$relative_OD),
         relative_OD = ifelse(dub == TRUE, NA, relative_OD)) %>%
  select(-time_h.y, -dub) 

com_imp <- com_imp %>%
  mutate(relative_OD_imp1 = na_ma(relative_OD, maxgap = 50), # Impute gaps
         relative_OD_imp2 = na_interpolation(relative_OD, maxgap = 50))

com_small <- filter(com_imp, !is.na(relative_OD))

growth_com <- growth_flight_clean %>%
  mutate(time_h = time_h - 4) %>% # Start point: temp => 30 Celcius
  filter(time_h >= 0) %>%
  difference_left_join(., growth_ground_clean, by = "time_h", max_dist = 0.25) %>% # Fuzzy match by hour
  rename(time_h = "time_h.x") %>%
  select(-time_h.y) %>%
  drop_na()

# Long format

rad_long <- pivot_longer(rad_clean, 
                         cols = ctrl_g:exp_totaln,
                         names_to = c("type", ".value"), 
                         names_pattern = "(.+)_(.+)")

growth_long1 <- growth_com %>%
  select(time_h, flight, ground1) %>%
  pivot_longer(cols = flight:ground1,
               names_to = c("type"))

growth_long2 <- growth_com %>%
  select(time_h, flight, ground2) %>%
  pivot_longer(cols = flight:ground2,
               names_to = c("type"))

growth_long3 <- growth_com %>%
  select(time_h, flight, ground3) %>%
  pivot_longer(cols = flight:ground3,
               names_to = c("type"))

growth_longm <- growth_com %>%
  select(time_h, flight, groundm) %>%
  pivot_longer(cols = flight:groundm,
               names_to = c("type"))

# Save

save(rad_clean, rad_long, 
     com_small, com_imp,
     growth_com, 
     growth_long1, growth_long2, growth_long3, growth_longm,
     file = "rad.Rdata")

# Plots

rad_long %>%
  ggplot(aes(y = log(g), x = type, color = type)) +
  geom_boxplot(outlier.size = 0.1, outlier.colour = "gray", notch = TRUE) +
  stat_summary(geom="text", fun = median,
               aes(label=sprintf("%1.2f", ..y..), color = factor(type)),
               position=position_nudge(x = 0.5), size = 4) +
  facet_grid(~ phases, scales = "free") +
  coord_cartesian(ylim = c(2.25, 6.25)) +
  theme(legend.title = element_blank()) +
  xlab("")

ggsave("p01.png", width = 10, height = 7)

rad_long %>%
ggplot(aes(x = time_h, y = g, color = type)) +
  geom_line(linetype = "dashed") +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

ggsave("p02.png", width = 10, height = 7)

rad_long %>%
  ggplot(aes(x = time_h, y = log(g), color = type)) +
  geom_line(alpha = 0.5)  +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

ggsave("p03.png", width = 10, height = 7)

rad_long %>%
  ggplot(aes(x = time_h, y = log(g), color = type)) +
  geom_line(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + 
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)") +
  facet_grid(~ phases, scales = "free_x")

ggsave("p04.png", width = 10, height = 7)

rad_long %>%
  ggplot(aes(x = time_h, y = log(g), color = type)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.25) + 
  geom_vline(xintercept = 10, linetype = "dotted") + 
  geom_vline(xintercept = 60, linetype = "dotted") +
  coord_cartesian(ylim=c(3.75, 4.5)) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

ggsave("p05.png", width = 10, height = 7)

rad_clean %>%
  ggplot(aes(x = time_h, y = delta)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.7) +
  geom_vline(xintercept = 10, linetype = "dotted") + 
  geom_vline(xintercept = 60, linetype = "dotted") +
  theme(legend.title = element_blank()) +
  ylab("g(ctrl) - g(exp)") +
  xlab("Time (in hours)")

ggsave("p06.png", width = 10, height = 7)

rad_long %>%
  ggplot(aes(x = g, y = gn, color = type)) +
  geom_point(size = 1) +
  theme(legend.title = element_blank())

ggsave("p07.png", width = 10, height = 7)

rad_long %>%
  filter(time_h <= 24) %>%
  ggplot(aes(x = time_h, y = total, color = type)) +
  geom_line(linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,
               rr.digits = 3) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

ggsave("p08.png", width = 10, height = 7)

rad_long %>%
  filter(time_h >= 240) %>%
  ggplot(aes(x = time_h, y = total, color = type)) +
  geom_line(linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,
               rr.digits = 3) +
  theme(legend.title = element_blank()) +
  xlab("Time (in hours)")

ggsave("p09.png", width = 10, height = 7)

com_small %>%
  ggplot() +
  geom_line(aes(x = time_h, y = relative_OD)) +
  geom_vline(xintercept = 10, linetype = "dotted") + 
  geom_vline(xintercept = 60, linetype = "dotted") +
  ylab("Relative optical density") +
  xlab("Time (in hours)")

ggsave("p10.png", width = 10, height = 7)

com_imp %>%
  ggplot() +
  geom_point(aes(x = time_h, y = relative_OD), color = "red", size = 0.5) +
  geom_line(aes(x = time_h, y = relative_OD_imp1)) +
  ylab("Relative optical density") +
  xlab("Time (in hours)")

ggsave("p11.png", width = 10, height = 7)

com_imp %>%
  filter(phases == 2) %>%
  select(time_h, relative_OD_imp1, delta) %>%
  ggpairs(lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
          diag = list(continuous = wrap("densityDiag")),
          upper = list(continuous = wrap("cor", stars = FALSE)),
          columnLabels = c("Time (in hours)", "Relative optical density", "g(ctrl) - g(exp)"))

ggsave("p12.png", width = 10, height = 7)

com_imp %>%
  filter(time_h <= 60) %>%
  select(time_h, relative_OD_imp1, delta_total) %>%
  ggpairs(lower = list(continuous = wrap("points", alpha = 0.3, size = 0.1)),
          diag = list(continuous = wrap("densityDiag")),
          upper = list(continuous = wrap("cor", stars = FALSE)),
          columnLabels = c("Time (in hours)", "Relative optical density", "total(g(ctrl) - g(exp))"))

ggsave("p13.png", width = 10, height = 7)
