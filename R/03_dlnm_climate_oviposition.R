# ====================================================================
# DLNM ANALYSIS: Climate Anomalies and Mosquito Oviposition
# With glmmTMB and manual predictions
# ====================================================================

# ====================================================================
# 1. LOAD LIBRARIES
# ====================================================================
library(dplyr)
library(lubridate)
library(ggplot2)
library(mgcv)
library(MASS)
library(dlnm)
library(DHARMa)
library(car)
library(splines)
library(gridExtra)
library(glmmTMB)
library(MuMIn)
library(cowplot) # For plot_grid in section 8
library(writexl)

# ====================================================================
# 2. HELPER FUNCTIONS
# ====================================================================

load_climate_data <- function(file_path, encoding = "latin1") {
  read.delim(file_path, header = TRUE, sep = ";", fileEncoding = encoding) %>%
    mutate(
      date = as.Date(Date, format = "%d/%m/%Y"),
      EW = epiweek(date),
      Year = year(date),
      tmax = as.numeric(ifelse(grepl("^[0-9.]+$", tmax), tmax, NA)),
      tmin = as.numeric(ifelse(grepl("^[0-9.]+$", tmin), tmin, NA)),
      prcp = as.numeric(ifelse(grepl("^[0-9.]+$", prcp), prcp, NA))
    )
}

create_seasons <- function(data) {
  data %>%
    mutate(Season = case_when(
      (Year == 2021 & EW >= 40) | (Year == 2022 & EW <= 13) ~ "2021-2022 (inter-epidemic)",
      (Year == 2022 & EW >= 40) | (Year == 2023 & EW <= 13) ~ "2022-2023 (epidemic)",
      (Year == 2023 & EW >= 40) | (Year == 2024 & EW <= 13) ~ "2023-2024 (epidemic)",
      (Year == 2024 & EW >= 40) | (Year == 2025 & EW <= 13) ~ "2024-2025 (inter-epidemic)"
    ))
}

# ====================================================================
# 3. DATA LOADING AND PREPROCESSING
# ====================================================================

# Load current climate data
current_climate <- load_climate_data("data/Santo_Tome_climate_2021-2024.csv") %>%
  create_seasons()

# Weekly climate means
weekly_means <- current_climate %>%
  group_by(Year, EW) %>%
  summarise(
    mean_tmax = mean(tmax, na.rm = TRUE),
    mean_tmin = mean(tmin, na.rm = TRUE),
    mean_prcp = mean(prcp, na.rm = TRUE),
    median_prcp = median(prcp, na.rm = TRUE),
    mean_tavg = mean((tmax + tmin) / 2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter((Year == 2021 & EW >= 40) | (Year %in% 2022:2024) | (Year == 2025 & EW <= 13))

# Load oviposition data
oviposition_data <- read.delim("data/Santo_Tome_OVI_weekly.csv",
                               header = TRUE, sep = ";", fileEncoding = "latin1") %>%
  filter(!is.na(Year)) %>%
  mutate(
    date_week = as.Date(Date, format = "%d/%m/%Y")
  ) %>%
  create_seasons() %>%
  filter(Season != "2024-2025 (inter-epidemic)") %>%
  filter((Year == 2021 & EW >= 40) | (Year %in% 2022:2024) | (Year == 2025 & EW <= 13))

# Join with current climate
oviposition_data <- left_join(oviposition_data, weekly_means, by = c("Year", "EW"))

# Load historical climate (1991-2020)
historical_climate <- load_climate_data("data/Santo_Tome_climate_1991-2020.csv")

historical_means <- historical_climate %>%
  group_by(EW) %>%
  summarise(
    hist_tmax = mean(tmax, na.rm = TRUE),
    hist_tmin = mean(tmin, na.rm = TRUE),
    hist_prcp = mean(prcp, na.rm = TRUE),
    hist_prcp_median = median(prcp, na.rm = TRUE),
    hist_tavg = mean((tmax + tmin) / 2, na.rm = TRUE),
    .groups = "drop"
  )

oviposition_data <- left_join(oviposition_data, historical_means, by = "EW")

# Calculate anomalies
oviposition_data <- oviposition_data %>%
  mutate(
    anom_tmax = mean_tmax - hist_tmax,
    anom_tmin = mean_tmin - hist_tmin,
    anom_prcp = mean_prcp - hist_prcp,
    anom_prcp_median = median_prcp - hist_prcp_median,
    anom_tavg = mean_tavg - hist_tavg
  )

# ====================================================================
# 4. PREPARE ANALYSIS DATASET
# ====================================================================

analysis_data <- oviposition_data %>%
  filter(month(date_week) %in% c(10:12, 1:5)) %>%
  mutate(
    EW_season = ifelse(EW >= 40, EW - 39, EW + 13),
    active_traps = ifelse(TOTAL >= 0, 1, 0),
    positive_traps = ifelse(TOTAL > 0, 1, 0)
  ) %>%
  group_by(Season) %>%
  mutate(EW_season = EW_season - min(EW_season) + 1) %>%
  ungroup()

# Summarize by season and week
summary_data <- analysis_data %>%
  group_by(Season, EW_season) %>%
  summarise(
    total_eggs = sum(TOTAL, na.rm = TRUE),
    total_active_traps = sum(active_traps, na.rm = TRUE),
    anom_tmax = first(anom_tmax),
    anom_tmin = first(anom_tmin),
    anom_prcp = first(anom_prcp),
    .groups = "drop"
  ) %>%
  filter(!is.na(Season)) %>%
  mutate(Season = factor(Season))

# ====================================================================
# 5. DLNM MODEL WITH glmmTMB
# ====================================================================

# Prepare data for DLNM
dlnm_data <- summary_data %>%
  dplyr::select(total_eggs, anom_tmin, anom_tmax, anom_prcp, 
                Season, EW_season, total_active_traps) %>%
  na.omit()

#Download table
dir.create("output", showWarnings = FALSE)
write_xlsx(dlnm_data, path = "output/Data_dlnm.xlsx")

# Check collinearity
cat("Correlation between Tmin and Tmax anomalies:", 
    cor(dlnm_data$anom_tmin, dlnm_data$anom_tmax), "\n\n")

vif_model <- lm(total_eggs ~ anom_tmin + anom_tmax + anom_prcp +
                  Season + ns(EW_season, df = 4),
                data = dlnm_data)
print("VIF Values:")
print(vif(vif_model))

# Create cross-basis matrices
cb_tmin <- crossbasis(dlnm_data$anom_tmin,
                      lag = 7,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns", df = 3))

cb_tmax <- crossbasis(dlnm_data$anom_tmax, lag = 7, 
                      argvar = list(fun = "ns", df = 3), # Use Spline
                      arglag = list(fun = "ns", df = 3))

cb_prcp <- crossbasis(dlnm_data$anom_prcp,
                      lag = 7,
                      argvar = list(fun = "lin"),
                      arglag = list(fun = "ns", df = 3))

# Fit model with glmmTMB and compare AIC
dlnm_model_1 <- glmmTMB(
  total_eggs ~ cb_tmin + cb_tmax + cb_prcp +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_2 <- glmmTMB(
  total_eggs ~ cb_tmin + cb_tmax + 
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_3 <- glmmTMB(
  total_eggs ~ cb_tmin + cb_prcp +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_4 <- glmmTMB(
  total_eggs ~ cb_tmax + cb_prcp +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_5 <- glmmTMB(
  total_eggs ~ cb_tmin +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_6 <- glmmTMB(
  total_eggs ~ cb_prcp +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

dlnm_model_7 <- glmmTMB(
  total_eggs ~ cb_tmax +
    Season + ns(EW_season, df = 4) +
    offset(log(total_active_traps + 1e-6)),
  family = nbinom1, # NB1 for underdispersion
  data = dlnm_data
)

model.sel(dlnm_model_1,dlnm_model_2,dlnm_model_3,dlnm_model_4,dlnm_model_5,dlnm_model_6,dlnm_model_7, rank = "AIC")

print("\n=== MODEL SUMMARY ===")
print(summary(dlnm_model_2))

# Model diagnostics
cat("\n=== MODEL DIAGNOSTICS ===\n")
sim_residuals <- simulateResiduals(dlnm_model_2)
plot(sim_residuals) # 

cat("\nDispersion Test:\n")
print(testDispersion(sim_residuals))

cat("\nTemporal Autocorrelation Test:\n")
sim_aggregated <- recalculateResiduals(sim_residuals, group = dlnm_data$EW_season)
print(testTemporalAutocorrelation(sim_aggregated, time = unique(dlnm_data$EW_season)))

# ====================================================================
# 6. FUNCTIONS FOR MANUAL PREDICTIONS
# ====================================================================

# Create the data frame containing all model variables
dlnm_data_clon <- cbind(
  dlnm_data, # Original variables (total_eggs, Season, EW_season, total_active_traps)
  cb_tmin,   # Cross-basis matrix for Tmin
  cb_tmax,   # Cross-basis matrix for Tmax
  cb_prcp    # Cross-basis matrix for Prcp
)

# Now, the formula inside glm.nb will not use the CB matrix names, 
# but the expanded columns (cb_tminv1.l1, etc.).

# Generate the part of the formula containing the expanded columns (cb_tminv1.l1, etc.)
cb_tmin_names <- colnames(cb_tmin)
cb_tmax_names <- colnames(cb_tmax)
cb_prcp_names <- colnames(cb_prcp)

# Create the new formula for glm.nb (which replaces cb_tmin with its columns)
formula_clon <- as.formula(
  paste("total_eggs ~ ", 
        paste(cb_tmin_names, collapse = " + "), " + ",
        paste(cb_tmax_names, collapse = " + "), " + ",
        paste(cb_prcp_names, collapse = " + "), " + ",
        "Season + ns(EW_season, df = 4) + offset(log(total_active_traps + 1e-6))")
)

# Extraction of glmmTMB Parameters
coefs_nb1 <- glmmTMB::fixef(dlnm_model_2)$cond
vcov_nb1 <- vcov(dlnm_model_2)$cond 

# Estimate theta (1/phi) from the nbinom1 dispersion parameter
log_phi <- dlnm_model_2$fit$par[1]
theta_nb <- 1 / exp(log_phi) 

# Fit the base model with glm.nb using the expanded formula
modelo_glm_clon <- MASS::glm.nb(
  formula = formula_clon,
  data = dlnm_data_clon, # We use the data frame containing ALL columns
  init.theta = theta_nb 
)

# Overwrite the base model parameters with those from glmmTMB (NB1)
modelo_glm_clon$coefficients <- coefs_nb1
modelo_glm_clon$vcov <- vcov_nb1 

# Ensure class for compatibility with crosspred
class(modelo_glm_clon) <- c("negbin", "glm", "lm")

library(dlnm) 
library(dplyr)
library(ggplot2)

# NOTE: crosspred needs the original cb_tmin matrix to know which coefficients to use.

# Generate predictions using the CLONED model 
# Tmin
pred_tmin <- crosspred(cb_tmin, model = modelo_glm_clon, by = 0.1, cen = 0, cumul = TRUE)

# Recreate the Tmin Anomalies vector 
min_tmin <- min(dlnm_data$anom_tmin) 
max_tmin <- max(dlnm_data$anom_tmin)
valores_anomalies_tmin <- seq(from = min_tmin, to = max_tmin, length.out = length(pred_tmin$allRRfit))

# Tmax
pred_tmax <- crosspred(cb_tmax, model = modelo_glm_clon, by = 0.1, cen = 0, cumul = TRUE)

# Recreate the Tmax Anomalies vector (X-axis, since allRRkey fails)
min_tmax <- min(dlnm_data$anom_tmax) 
max_tmax <- max(dlnm_data$anom_tmax)
valores_anomalies_tmax <- seq(from = min_tmax, to = max_tmax, length.out = length(pred_tmax$allRRfit))

# ====================================================================
# 7. VISUALIZATION FUNCTIONS
# ====================================================================

# Exposure–response plot (overall cumulative effect)
plot_exposure_response <- function(pred_obj, var_name, color, label_letter = NULL) {
  df <- data.frame(
    anomaly = pred_obj$predvar,
    RPEC = pred_obj$allRRfit,
    lower = pred_obj$allRRlow,
    upper = pred_obj$allRRhigh
  )
  
  ggplot(df, aes(x = anomaly, y = RPEC)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = color, alpha = 0.3) +
    geom_line(color = color, linewidth = 1.2) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    labs(
      title = paste0(if (!is.null(label_letter)) paste0(label_letter, " "), 
                     "Cumulative effect of ", var_name, " anomaly"),
      x = paste(var_name, "anomaly (°C)"),
      y = "Relative predicted\negg counts"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    scale_y_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      
      axis.title = element_text(face = "plain", size = 16),
      axis.title.x = element_text( size = 16,margin = margin(t = 5, r = 0, b = 0, l = 0)), 
      axis.title.y = element_text( size = 16,margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
      axis.line.x = element_line(color = "black", linetype = "solid"),
      panel.grid.minor = element_blank()
    )
}

# Lag–response plot (effect by lag for a given anomaly)
plot_lag_response <- function(pred_obj, var_name, anomaly_val, color, label_letter = NULL) {
  idx <- which.min(abs(pred_obj$predvar - anomaly_val))
  actual_val <- pred_obj$predvar[idx]
  
  nlags <- ncol(pred_obj$matRRfit) - 1
  lag_vec <- 0:nlags
  
  df <- data.frame(
    lag = lag_vec,
    RPEC = pred_obj$matRRfit[idx, ],
    lower = pred_obj$matRRlow[idx, ],
    upper = pred_obj$matRRhigh[idx, ]
  )
  
  ggplot(df, aes(x = lag, y = RPEC)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = color, alpha = 0.3) +
    geom_line(color = color, linewidth = 1.2) +
    geom_point(color = color, size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
    scale_x_continuous(breaks = lag_vec) +
    labs(
      subtitle = paste0(if (!is.null(label_letter)) paste0(label_letter, " "), 
                     "Lag–response for ", var_name, " = ", round(actual_val, 2), " °C"),
      x = "Lag (weeks)",
      y = "Relative predicted\negg counts"
    ) +
     scale_y_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    scale_x_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      
      axis.title = element_text(face = "plain", size = 16),
      axis.title.x = element_text( size = 16,margin = margin(t = 5, r = 0, b = 0, l = 0)), 
      axis.title.y = element_text( size = 16,margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
      axis.line.x = element_line(color = "black", linetype = "solid"),
      panel.grid.minor = element_blank()
    )
}

# 3D contour surface (anomaly × lag)
plot_3d_contour <- function(pred_obj, var_name, color_palette, label_letter = NULL) {
  nlags <- ncol(pred_obj$matRRfit) - 1
  
  df_3d <- expand.grid(
    lag = 0:nlags,
    anomaly = pred_obj$predvar
  )
  
  df_3d$RPEC <- as.vector(t(pred_obj$matRRfit))
  
  ggplot(df_3d, aes(x = lag, y = anomaly, fill = RPEC)) +
    geom_tile() +
    geom_contour(aes(z = RPEC), color = "white", alpha = 0.5, linewidth = 0.3) +
    scale_fill_gradientn(
      colors = color_palette,
      name = "Relative predicted\negg counts",
      limits = c(min(df_3d$RPEC, na.rm = TRUE), max(df_3d$RPEC, na.rm = TRUE))
    ) +
    scale_x_continuous(breaks = 0:nlags) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    scale_x_continuous(
      expand = expansion(mult = c(0, 0))
    )+
    labs(
      x = "Lag (weeks)",
      y = paste(var_name, "anomaly (°C)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      
      axis.title = element_text(face = "plain", size = 16),
      axis.title.x = element_text( size = 16,margin = margin(t = 5, r = 0, b = 0, l = 0)), 
      axis.title.y = element_text( size = 16,margin = margin(t = 0, r = 5, b = 0, l = 0)),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
      axis.line.x = element_line(color = "black", linetype = "solid"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# ====================================================================
# 8. GENERATE ALL VISUALIZATIONS
# ====================================================================

cat("\n=== GENERATING VISUALIZATIONS ===\n")

# Exposure–response plots
p1 <- plot_exposure_response(pred_tmin, "Tmin", "#2E86AB", label_letter = "(A)")
p2 <- plot_exposure_response(pred_tmax, "Tmax", "#A23B72", label_letter = "(B)")
p1
p2

# Lag–response plots (for specific anomaly values)
p4 <- plot_lag_response(pred_tmin, "Tmin", 2, "#2E86AB")
p5 <- plot_lag_response(pred_tmax, "Tmax", -2, "#A23B72")
p4
p5

# 3D contour plots
p7 <- plot_3d_contour(pred_tmin, "Tmin", c("#0A1F44", "#2E86AB", "#A6E1FA"))
p8 <- plot_3d_contour(pred_tmax, "Tmax", c("#3D0814", "#A23B72", "#F18F01"))

# Display
p1; p2
p4; p5
p7; p8 # 

# Optional: arrange them side by side
fig4 = plot_grid(p7 + labs(tag = "(a)"), p8 + labs(tag= "(b)"), p4 + labs(tag = "(c)"), p5 + labs(tag = "(d)"), ncol = 2)
fig4

# ====================================================================
# 9. KEY FINDINGS SUMMARY
# ====================================================================

cat("\n=== KEY FINDINGS ===\n")

# Function to find Relative predicted egg counts at a specific value
get_RPEC_at_value <- function(pred_obj, value, lag = 0) {
  idx <- which.min(abs(pred_obj$predvar - value))
  actual_val <- pred_obj$predvar[idx]
  RPEC <- pred_obj$matRRfit[idx, lag + 1]
  lower <- pred_obj$matRRlow[idx, lag + 1]
  upper <- pred_obj$matRRhigh[idx, lag + 1]
  
  return(list(value = actual_val, RPEC = RPEC, lower = lower, upper = upper))
}

# Example: Relative predicted egg counts for specific anomalies
tmin_res <- get_RPEC_at_value(pred_tmin, 2, lag = 0)
cat(sprintf("Tmin Anomaly = %.2f°C (lag 0): RPEC = %.3f (IC95%%: %.3f - %.3f)\n",
            tmin_res$value, tmin_res$RPEC, tmin_res$lower, tmin_res$upper))

tmax_res <- get_RPEC_at_value(pred_tmax, -2, lag = 0)
cat(sprintf("Tmax Anomaly = %.2f°C (lag 0): RPEC = %.3f (IC95%%: %.3f - %.3f)\n",
            tmax_res$value, tmax_res$RPEC, tmax_res$lower, tmax_res$upper))

# Cumulative effect
tmin_cumul <- get_RPEC_at_value(pred_tmin, 2, lag = 7)
cat(sprintf("\nCumulative effect (lag 0-7) Tmin = %.2f°C: RPEC = %.3f (IC95%%: %.3f - %.3f)\n",
            tmin_cumul$value, tmin_cumul$RPEC, tmin_cumul$lower, tmin_cumul$upper))

tmax_cumul <- get_RPEC_at_value(pred_tmax, -2, lag = 7)
cat(sprintf("Cumulative effect (lag 0-7) Tmax = %.2f°C: RPEC = %.3f (IC95%%: %.3f - %.3f)\n",
            tmax_cumul$value, tmax_cumul$RPEC, tmax_cumul$lower, tmax_cumul$upper))

cat("\n=== COMPLETE ANALYSIS ===\n")

# Save the panel
ggsave("output/Fig4.tiff",
       plot = fig4,
       width = 19,
       height = 12,
       dpi = 300)
