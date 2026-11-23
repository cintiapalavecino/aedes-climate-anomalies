# ====================================================================
# 01. CLIMATE ANOMALIES (Santo Tome, Santa Fe, Argentina)
# Weekly Tmin, Tmax and precipitation anomalies (2021-2024) relative to
# the 1991-2020 climatology, modelled with GAMs and quantile GAMs.
# Author: Cintia C. Palavecino
# ====================================================================

############################ Libraries #############################
library(dplyr)
library(cowplot)
library(lubridate)
library(ggplot2)
library(trend)
library(purrr)
library(writexl)
library(scales) # To format the dates on the x-axis

############################ Load Current Data #############################
Data_2021_2024 <- read.delim("data/Santo_Tome_climate_2021-2024.csv",
                             header = TRUE, sep = ";", fileEncoding = "latin1")

# Ensure date format and calculate epidemiological week
Data_2021_2024 <- Data_2021_2024 %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"),
         EW = epiweek(Date),
         Year = year(Date),
         tmax = as.numeric(ifelse(grepl("^[0-9.]+$", tmax), tmax, NA)),
         tmin = as.numeric(ifelse(grepl("^[0-9.]+$", tmin), tmin, NA)),
         prcp = as.numeric(ifelse(grepl("^[0-9.]+$", prcp), prcp, NA)))

# Create seasons (correcting labels to English)
Data_2021_2024 <- Data_2021_2024 %>%
  mutate(Season = case_when(
    (Year == 2021 & EW >= 38) | (Year == 2022 & EW <= 37) ~ "2021 - 2022 (non-epidemic)",
    (Year == 2022 & EW >= 38) | (Year == 2023 & EW <= 37) ~ "2022 - 2023 (epidemic)",
    (Year == 2023 & EW >= 38) | (Year == 2024 & EW <= 37) ~ "2023 - 2024 (epidemic)"
  ))

# Calculate current climatic means by epidemiological week and year
climate_media <- Data_2021_2024 %>%
  group_by(Year, EW) %>%
  summarise(media_tmax = mean(tmax, na.rm = TRUE),
            media_tmin = mean(tmin, na.rm = TRUE),
            media_prcp = mean(prcp, na.rm = TRUE),
            median_prcp = median(prcp, na.rm = TRUE),
            media_tprom = mean((tmax + tmin) / 2, na.rm = TRUE)) %>%
  ungroup() %>%
  filter((Year == 2021 & EW >= 38) | (Year %in% c(2022, 2023)) | (Year == 2024 & EW <= 20))

############################ Load OVI Data #############################
Data <- read.delim("data/Santo_Tome_OVI_weekly.csv",
                   header = TRUE, sep = ";", fileEncoding = "latin1") %>%
  filter(!is.na(Year))

# Create seasons
Data <- Data %>%
  mutate(Season = case_when(
    (Year == 2021 & EW >= 38) | (Year == 2022 & EW <= 37) ~ "2021 - 2022 (non-epidemic)",
    (Year == 2022 & EW >= 38) | (Year == 2023 & EW <= 37) ~ "2022 - 2023 (epidemic)",
    (Year == 2023 & EW >= 38) | (Year == 2024 & EW <= 37) ~ "2023 - 2024 (epidemic)"
  ))

# Filter Data to include only data from EW 38 of 2021 up to EW 20 of 2024
Data <- Data %>%
  filter((Year == 2021 & EW >= 38) | (Year %in% c(2022, 2023)) | (Year == 2024 & EW <= 22))

# Convert Year of Data to numeric
Data$Year <- as.numeric(as.character(Data$Year))

# Join with current climatic means
Data <- left_join(Data, climate_media, by = c("Year", "EW"))

# Calculate the average between Tmin and Tmax in the current records
medias_temp <- Data %>%
  group_by(Year, EW) %>%
  summarise(media_tprom = mean((media_tmax + media_tmin) / 2, na.rm = TRUE)) %>%
  ungroup()

# Join with medias_temp
Data <- left_join(Data, medias_temp, by = c("Year", "EW"))

############################ Load Historical Data #############################
Data_1990_2020 <- read.delim("data/Santo_Tome_climate_1991-2020.csv",
                             header = TRUE, sep = ";", fileEncoding = "latin1")

# Process historical data
climate_media_1990_2020 <- Data_1990_2020 %>%
  mutate(Date = as.Date(Date, format = "%d/%m/%Y"),
         EW = epiweek(Date),
         Year = year(Date),
         tmax = as.numeric(ifelse(grepl("^[0-9.]+$", tmax), tmax, NA)),
         tmin = as.numeric(ifelse(grepl("^[0-9.]+$", tmin), tmin, NA)),
         prcp = as.numeric(ifelse(grepl("^[0-9.]+$", prcp), prcp, NA))) %>%
  group_by(EW) %>%
  summarise(media_tmax_hist = mean(tmax, na.rm = TRUE),
            media_tmin_hist = mean(tmin, na.rm = TRUE),
            media_prcp_hist = mean(prcp, na.rm = TRUE),
            median_prcp_hist = median(prcp, na.rm = TRUE),
            media_tprom_hist = mean((tmax + tmin) / 2, na.rm = TRUE)) %>%
  ungroup()

# Join with Data
Data <- left_join(Data, climate_media_1990_2020, by = "EW")
names(Data)

# Choose which to keep and rename
Data <- Data %>%
  mutate(media_tprom = coalesce(media_tprom.x, media_tprom.y)) %>%
  dplyr::select(-any_of(c("media_tprom.x", "media_tprom.y")))

# Calculate climatic anomalies
Data <- Data %>%
  mutate(anomalies_tmax = media_tmax - media_tmax_hist,
         anomalies_tmin = media_tmin - media_tmin_hist,
         anomalies_prcp = media_prcp - media_prcp_hist,
         anomalies_prcp_median = median_prcp - median_prcp_hist,
         anomalies_tprom = media_tprom - media_tprom_hist)

library(ggplot2)
library(scales)

#anomalies tmax
Sys.setlocale("LC_TIME", "en_US.UTF-8")  

Data$Season <- dplyr::recode(Data$Season,
                             "2021 - 2022 (no epidémica)" = "2021 - 2022 (non-epidemic)",
                             "2022 - 2023 (epidémica)" = "2022 - 2023 (epidemic)",
                             "2023 - 2024 (epidémica)" = "2023 - 2024 (epidemic)")

# Filter to include up to May (month 5) and up to EW 22
Filtered_data <- Data %>%
  filter(lubridate::month(as.Date(Date, format = "%d/%m/%Y")) %in% c(10,11,12,1,2,3,4,5)) %>%
  filter((Year == 2021 & EW >= 38) | (Year %in% c(2022, 2023)) | (Year == 2024 & EW <= 22))

# Create the EW_Season variable
Filtered_data <- Filtered_data %>%
  mutate(
    EW_Season = case_when(
      EW >= 38 ~ EW - 37,
      EW < 38  ~ EW + 15
    )
  ) %>%
  group_by(Season) %>%
  mutate(EW_Season = EW_Season - min(EW_Season) + 1) %>%
  ungroup()

Filtered_data$Season <- as.factor(Filtered_data$Season)

# Check range
table(Filtered_data$EW_Season)

# Check anomalies
Filtered_data_weekly <- Filtered_data %>%
  group_by(Season, Year, EW, EW_Season) %>%
  summarise(
    anomalies_tmax = first(anomalies_tmax),
    anomalies_tmin = first(anomalies_tmin),
    anomalies_prcp = first(anomalies_prcp),
    anomalies_prcp_median = first(anomalies_prcp_median),
    anomalies_tprom = first(anomalies_tprom),
    
    .groups = "drop"
  )

#################################GAM ANOMALIES#################################
# season as factor
Filtered_data_weekly$Season <- as.factor(Filtered_data_weekly$Season)

#Download table
dir.create("output", showWarnings = FALSE)
write_xlsx(Filtered_data_weekly, path = "output/Data_anomalies.xlsx")

####GAM tmax_anomalies####
library(mgcv)
library(dplyr)

# Assuming your dataframe is named Data and has a Tmax_Anomaly column
# Fit the GAM model
gam_anomalies_tmax <- gam(anomalies_tmax ~ s(EW_Season, by = Season, k = 7) + Season, 
                          data = Filtered_data_weekly, method = "REML")

# Model summary
summary(gam_anomalies_tmax)

#install.packages("mgcViz")
library(mgcViz)
library(DHARMa)
simulateResiduals(fittedModel = gam_anomalies_tmax, plot = T, refit = F)

par(mfrow = c(2, 2))
gam.check (gam_anomalies_tmax)
par(mfrow = c(1, 1))

library(ggplot2)
library(mgcv)
library(dplyr)

# Create a new prediction dataset with all combinations of EW_Season and Season
pred_data_anomalies_tmax <- expand.grid(
  EW_Season = seq(min(Filtered_data_weekly$EW_Season), max(Filtered_data_weekly$EW_Season), length.out = 100),
  Season = unique(Filtered_data_weekly$Season)
)

# Add predictions and their confidence intervals
pred_data_anomalies_tmax$Fit <- predict(gam_anomalies_tmax, newdata = pred_data_anomalies_tmax, type = "response", se.fit = TRUE)$fit
pred_data_anomalies_tmax$EW_fit <- predict(gam_anomalies_tmax, newdata = pred_data_anomalies_tmax, type = "response", se.fit = TRUE)$se.fit

####GAM tmin_anomalies####
# Create the EW_Season variable to start at 0
# Fit the GAM model
gam_anomalies_tmin <- gam(anomalies_tmin ~ s(EW_Season, by = Season, k = 7) + Season, 
                          data = Filtered_data_weekly, method = "REML")

# Model summary
summary(gam_anomalies_tmin)

library(mgcViz)
library(DHARMa)
simulateResiduals(fittedModel = gam_anomalies_tmin, plot = T, refit = F)

par(mfrow = c(2, 2))
gam.check (gam_anomalies_tmin)
par(mfrow = c(1, 1))

library(ggplot2)
library(mgcv)
library(dplyr)

# Create a new prediction dataset with all combinations of EW_Season and Season
pred_data_anomalies_tmin <- expand.grid(
  EW_Season = seq(min(Filtered_data_weekly$EW_Season), max(Filtered_data_weekly$EW_Season), length.out = 100),
  Season = unique(Filtered_data_weekly$Season)
)

# Add predictions and their confidence intervals
pred_data_anomalies_tmin$Fit <- predict(gam_anomalies_tmin, newdata = pred_data_anomalies_tmin, type = "response", se.fit = TRUE)$fit
pred_data_anomalies_tmin$EW_fit <- predict(gam_anomalies_tmin, newdata = pred_data_anomalies_tmin, type = "response", se.fit = TRUE)$se.fit

####GAM precipitation_anomalies####
# Fit the GAM model to compare medians
library(qgam)
qgam_anomalies_prcp_median <- qgam(anomalies_prcp ~ s(EW_Season, by = Season, k = 8) + Season,
                                   data = Filtered_data_weekly,
                                   qu = 0.5)
summary(qgam_anomalies_prcp_median)

# Create a new prediction dataset with all combinations of EW_Season and Season
pred_data_anomalies_prcp_median <- expand.grid(
  EW_Season = seq(min(Filtered_data_weekly$EW_Season), max(Filtered_data_weekly$EW_Season), length.out = 100),
  Season = unique(Filtered_data_weekly$Season)
)

# Add predictions and their confidence intervals
pred_data_anomalies_prcp_median$Fit <- predict(qgam_anomalies_prcp_median, newdata = pred_data_anomalies_prcp_median, type = "response", se.fit = TRUE)$fit
pred_data_anomalies_prcp_median$EW_fit <- predict(qgam_anomalies_prcp_median, newdata = pred_data_anomalies_prcp_median, type = "response", se.fit = TRUE)$se.fit

# Installation (if necessary): You might need the 'viridis' or 'RColorBrewer' library
# install.packages("viridis") 
# library(viridis)

# ====================================================================
# COLOR-BLIND SAFE AND HIGH-CONTRAST PALETTE
# We will use a safe palette, for example, "orange" and "cyan/light blue"
# ====================================================================

# Define color-blind safe colors (e.g., from the Okabe-Ito palette or similar)
COLOR_POSITIVE <- "#E69F00" # Burnt orange/Gold (simulates Red)
COLOR_NEGATIVE <- "#0072B2" # Bright blue/Cyan (simulates Blue)

# Define shapes to reinforce distinction without color
SHAPE_POSITIVE <- 19 # Filled circle
SHAPE_NEGATIVE <- 17 # Filled triangle (high shape contrast)

# Calculate the maximum EW_Season
max_EW_Season <- max(Filtered_data_weekly$EW_Season, na.rm = TRUE)

# --------------------------------------------------------------------
# A. TMIN PLOT (Adapted)
# --------------------------------------------------------------------
plot_combined_anomalies_tmin <- ggplot() +
  geom_segment(data = Filtered_data_weekly,
               aes(x = EW_Season, xend = EW_Season,
                   y = 0, yend = anomalies_tmin,
                   color = anomalies_tmin > 0), size = 1.2) +
  geom_point(data = Filtered_data_weekly,
             aes(x = EW_Season, y = anomalies_tmin,
                 color = anomalies_tmin > 0,
                 # Reinforce distinction with different shape for Positive/negative
                 shape = anomalies_tmin > 0), size = 2) +
  geom_line(data = pred_data_anomalies_tmin,
            aes(x = EW_Season, y = Fit),
            color = "black", size = 1.2) + # GAM curve always black (high contrast)
  geom_ribbon(data = pred_data_anomalies_tmin,
              aes(x = EW_Season, ymin = Fit - 2 * EW_fit, ymax = Fit + 2 * EW_fit),
              alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  facet_wrap(~ Season, scales = "free_x",
             labeller = labeller(Season = label_value)) +
  # === ADAPTATION: Use of safe palette and shape assignment ===
  scale_color_manual(values = c("FALSE" = COLOR_NEGATIVE, "TRUE" = COLOR_POSITIVE)) +
  scale_shape_manual(values = c("FALSE" = SHAPE_NEGATIVE, "TRUE" = SHAPE_POSITIVE)) +
  # ===============================================================
scale_x_continuous(
  breaks = seq(1, max_EW_Season, by = 3),
  labels = function(x) {
    epi_week <- ifelse(x <= 13, x + 39, x - 13)
    as.character(epi_week)
  },
  limits = c(1, max_EW_Season),
  
) +
  labs(y = "Tmin anomaly (°C)", x = "") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x.top = element_text(size = 16, angle = 0, hjust = 0.5),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "none", # Distinction is now based on color AND shape
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "plain"),
    axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
    axis.line.x = element_line(color = "black", linetype = "solid"),
    plot.margin = margin(5, 5, 0, 5)
  )


print(plot_combined_anomalies_tmin)

# --------------------------------------------------------------------
# B. TMAX PLOT (Adapted)
# --------------------------------------------------------------------
plot_combined_anomalies_tmax <- ggplot() +
  geom_segment(data = Filtered_data_weekly,
               aes(x = EW_Season, xend = EW_Season,
                   y = 0, yend = anomalies_tmax,
                   color = anomalies_tmax > 0), size = 1.2) +
  geom_point(data = Filtered_data_weekly,
             aes(x = EW_Season, y = anomalies_tmax,
                 color = anomalies_tmax > 0,
                 # Reinforce distinction with shape
                 shape = anomalies_tmax > 0), size = 2) +
  geom_line(data = pred_data_anomalies_tmax,
            aes(x = EW_Season, y = Fit),
            color = "black", size = 1.2) +
  geom_ribbon(data = pred_data_anomalies_tmax,
              aes(x = EW_Season, ymin = Fit - 2 * EW_fit, ymax = Fit + 2 * EW_fit),
              alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  facet_wrap(~ Season, scales = "free_x") +
  # === ADAPTATION: Safe palette and shape assignment ===
  scale_color_manual(values = c("FALSE" = COLOR_NEGATIVE, "TRUE" = COLOR_POSITIVE)) +
  scale_shape_manual(values = c("FALSE" = SHAPE_NEGATIVE, "TRUE" = SHAPE_POSITIVE)) +
  # =======================================================
scale_x_continuous(
  breaks = seq(1, max_EW_Season, by = 3),
  labels = function(x) {
    epi_week <- ifelse(x <= 13, x + 39, x - 13)
    as.character(epi_week)
  },
  limits = c(1, max_EW_Season),
  sec.axis = sec_axis(~ .,
                      breaks = c(2.5, 6.5, 11, 15.5, 20, 24.5, 28.5, 32.5),
                      labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"),
                      name = "")
) +
  labs(y = "Tmax anomaly (°C)", x = "") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(face = "plain"),
    axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
    axis.line.x = element_line(color = "black", linetype = "solid"),
    
    plot.margin = margin(0, 5, 0, 5)
  )

print(plot_combined_anomalies_tmax)

# --------------------------------------------------------------------
# C. PRECIPITATION PLOT (Adapted)
# --------------------------------------------------------------------
plot_combined_anomalies_prcp_median <- ggplot() +
  geom_segment(data = Filtered_data_weekly,
               aes(x = EW_Season, xend = EW_Season,
                   y = 0, yend = anomalies_prcp,
                   color = anomalies_prcp > 0), size = 1.2) +
  geom_point(data = Filtered_data_weekly,
             aes(x = EW_Season, y = anomalies_prcp,
                 color = anomalies_prcp > 0,
                 # Reinforce distinction with shape
                 shape = anomalies_prcp > 0), size = 2) +
  geom_line(data = pred_data_anomalies_prcp_median,
            aes(x = EW_Season, y = Fit),
            color = "black", size = 1.2) +
  geom_ribbon(data = pred_data_anomalies_prcp_median,
              aes(x = EW_Season, ymin = Fit - 2 * EW_fit, ymax = Fit + 2 * EW_fit),
              alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  facet_wrap(~ Season, scales = "free_x",
             labeller = labeller(Season = label_value)) +
  # === ADAPTATION: Safe palette and shape assignment ===
  scale_color_manual(values = c("FALSE" = COLOR_NEGATIVE, "TRUE" = COLOR_POSITIVE)) +
  scale_shape_manual(values = c("FALSE" = SHAPE_NEGATIVE, "TRUE" = SHAPE_POSITIVE)) +
  # =======================================================
scale_x_continuous(
  breaks = seq(1, max_EW_Season, by = 3),
  labels = function(x) {
    epi_week <- ifelse(x <= 13, x + 39, x - 13)
    as.character(epi_week)
  },
  limits = c(1, max_EW_Season),
  sec.axis = sec_axis(~ .,
                      breaks = c(2.5, 6.5, 11, 15.5, 20, 24.5, 28.5, 32.5),
                      labels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"),
                      name = "")
) +
  labs(y = "Precipitation anomaly (mm)",
       x = "Epidemiological week") +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(hjust = 1, size = 16),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    
    axis.title = element_text(face = "plain"),
    axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid = element_blank(),
    panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
    axis.line.x = element_line(color = "black", linetype = "solid"),
    
    plot.margin = margin(0, 5, 5, 5)
  )

print(plot_combined_anomalies_prcp_median)

# Combine plots (no changes, since it uses the updated objects)
plot_combined_anomalies_optimized <- plot_grid(
  plot_combined_anomalies_tmin + labs(subtitle = "(a)"),
  plot_combined_anomalies_tmax + labs(subtitle = "(b)"),
  plot_combined_anomalies_prcp_median + labs(subtitle = "(c)"),
  ncol = 1,
  nrow = 3,
  align = "v",
  rel_heights = c(1, 0.85, 1)
)

plot_combined_anomalies_final <- plot_grid(
  ggdraw() + draw_label("Temperature anomalies: observed data and GAM predictions",
                        fontface = 'bold', size = 18),
  plot_combined_anomalies_optimized,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

print(plot_combined_anomalies_final)

# Save the panel
dir.create("output", showWarnings = FALSE)


ggsave("output/Fig2.jpg", 
       plot = plot_combined_anomalies_final,
       width = 19, 
       height = 12, 
       dpi = 300)
