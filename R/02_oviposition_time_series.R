# ====================================================================
# CLIMATIC ANOMALIES AND AEDES AEGYPTI OVIPOSITION
# Time series plot of weekly egg count by season.
# ====================================================================

# ====================================================================
# 1. LOAD LIBRARIES AND DEFINE HELPER FUNCTIONS
# ====================================================================
library(dplyr)
library(lubridate)
library(ggplot2)
library(mgcv)
library(MuMIn)
library(cowplot)
library(scales)
library(fitdistrplus)
library(gridExtra)
library(purrr)
library(trend)
library(writexl)

# --- HELPER FUNCTIONS ---

# Function to create season labels (based on EW 38 start/end)
create_seasons <- function(data) {
  data %>%
    mutate(Season = case_when(
      (Year == 2021 & EW >= 38) | (Year == 2022 & EW <= 37) ~ "2021 - 2022 (non-epidemic)",
      (Year == 2022 & EW >= 38) | (Year == 2023 & EW <= 37) ~ "2022 - 2023 (epidemic)",
      (Year == 2023 & EW >= 38) | (Year == 2024 & EW <= 37) ~ "2023 - 2024 (epidemic)",
      (Year == 2024 & EW >= 38) | (Year == 2025 & EW <= 37) ~ "2024 - 2025 (non-epidemic)"
    ))
}

# ====================================================================
# 2. DATA LOADING AND INITIAL PROCESSING
# ====================================================================

# Load oviposition data
Data <- read.delim("data/Santo_Tome_OVI_weekly.csv",
                   header = TRUE, sep = ";", fileEncoding = "latin1") %>%
  filter(!is.na(Year))

# Clean and prepare the data frame for seasonal analysis
Data_filtered <- Data %>%
  # 1. Ensure Date, EW, and Year are calculated correctly for 'create_seasons'
  mutate(
    Date = as.Date(Date, format = "%d/%m/%Y"),
    EW = epiweek(Date),
    Year = year(Date)
  ) %>%
  # 2. Apply the 'create_seasons' function to define the Season column
  create_seasons() %>%
  filter(Season != "2024 - 2025 (non-epidemic)") %>%
  filter((Year == 2021 & EW >= 40) | (Year %in% 2022:2024) | (Year == 2025 & EW <= 13))%>%
  # 3. Filter for active transmission months (Oct to May)
  filter(lubridate::month(Date) %in% c(10, 11, 12, 1, 2, 3, 4, 5)) %>%
  # 4. Calculate season-specific variables
  mutate(
    Season = as.factor(Season), 
    EW_season = case_when(
      EW >= 38 ~ EW - 37,
      EW < 38  ~ EW + 15 
    ),
    Active_ovitraps = ifelse(TOTAL >= 0, 1, 0),
    positive_ovitraps = ifelse(TOTAL > 0, 1, 0)
  ) %>%
  # 5. Normalize EW_season within each Season (starts at 1)
  group_by(Season) %>%
  mutate(EW_season = EW_season - min(EW_season) + 1) %>%
  ungroup()

# ====================================================================
# 3. DATA AGGREGATION AND SUMMARY
# ====================================================================

# Summarize total eggs and active traps by Season and EW_season
# NOTE: Keeping the original variable name 'resumen_huevos_total' as requested.
resumen_huevos_total <- Data_filtered %>%
  # Filter for months (although already done above, retaining original logic for safety)
  filter(lubridate::month(Date) %in% c(10, 11, 12, 1, 2, 3, 4, 5)) %>% 
  group_by(Season, EW_season) %>%
  summarise(
    Total_eggs = sum(TOTAL, na.rm = TRUE),
    Total_active_ovitraps_SE = sum(Active_ovitraps, na.rm = TRUE), # SE means 'Season-specific' effort
    .groups = "drop"
  ) %>%
  filter(!is.na(Season) & !is.na(Total_eggs))

#Download table
dir.create("output", showWarnings = FALSE)
write_xlsx(resumen_huevos_total, path = "output/Data_oviposition.xlsx")

# -------------------------------
# PLOT GENERATION (Original Code Restored)
# -------------------------------
Fig3 <- ggplot(
  resumen_huevos_total,
  aes(x = EW_season, y = Total_eggs)
) +
  
  geom_line(color = "#2E6F40", linewidth = 1) +
  geom_point(color = "#2E6F40", size = 2) +
  
  facet_wrap(
    ~ Season,
    scales = "free_x",
    labeller = labeller(Season = label_value)
  ) +
  
  scale_x_continuous(
    breaks = seq(1, max(resumen_huevos_total$EW_season), by = 3),
    labels = function(x){
      
      epi_week <- ifelse(x <= 13, x + 39, x - 13)
      
      mes <- cut(
        x,
        breaks = c(0,4,8,13,17,22,26,30,34),
        labels = c("Oct","Nov","Dec","Jan",
                   "Feb","Mar","Apr","May")
      )
      
      paste0(epi_week, "\n", mes)
    },
    limits = c(1, max(resumen_huevos_total$EW_season)),
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Epidemiological week",
    y = "Weekly total egg counts"
  ) +
  
  theme_minimal(base_size = 18) +
  theme(axis.title = element_text(face = "plain"),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)), 
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.text = element_text(color = "black"),
        panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "#D9D9D9", linetype = "longdash",linewidth = 0.2),
        axis.line.x = element_line(color = "black", linetype = "solid"),
        strip.text = element_text(size = 18, face = "bold"),
        plot.margin = margin(10, 10, 40, 10)
  )


Fig3 # Display plot

# ====================================================================
# 4. EXPORT PLOT (OPTIONAL)
# ====================================================================
# Set working directory and save plot (Adjusted to 500 dpi for journal requirement.)
# NOTE: Replace the file path with the desired output location.
#setwd("") 

ggsave("output/Fig3.tiff", 
       plot = Fig3,
       width = 14, 
       height = 9, 
       dpi = 300)
