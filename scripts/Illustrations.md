#  Figures 
Here we show all the codes that were used to create the Figures of the results of the fermentation experiment 

## 01 Figure Vol% of all approaches (calculated by °P) with ribbons
```r
# Load nessecary Packages
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load data
data <- read.csv("~/Desktop/Plato_pH_alk.csv", sep=";", dec=",")

# Remove unnecessary columns
data <- data[, c("Ansatz", "Tag", "Vol_prozent")]

# Group replicates and rename categories
data <- data %>%
  filter(!Ansatz %in% c(7, 8)) %>% # Ignore Ansatz 7 and 8
  mutate(Ansatz = case_when(
    Ansatz %in% c(1, 11, 16) ~ "1",
    Ansatz %in% c(3, 12, 17) ~ "3",
    Ansatz %in% c(5, 13, 18) ~ "5",
    Ansatz %in% c(6, 14, 19) ~ "6",
    Ansatz %in% c(10, 15, 20) ~ "8",
    Ansatz == 9 ~ "7",
    Ansatz == 2 ~ "2",
    Ansatz == 4 ~ "4",
    TRUE ~ as.character(Ansatz)
  )) %>%
  group_by(Ansatz, Tag) %>%
  summarise(
    Vol_mean = mean(Vol_prozent, na.rm=TRUE),
    Vol_sd = sd(Vol_prozent, na.rm=TRUE)
  ) %>%
  ungroup()

# Define custom theme 
theme_custom <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),  
    axis.title = element_text(size = 30, face = "bold"),
    axis.text = element_text(size = 24),
    strip.text = element_text(size = 30, face = "bold", color = "black"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid.major = element_line(size = 0.5, color = "gray"),
    panel.grid.minor = element_line(size = 0.2, color = "gray"),
    plot.title = element_text(size = 27, face = "bold", hjust = 0.5)
  )

# Define custom color palette for the Approaches
ansatz_colors <- c("1" = "#D2B48C",   
                   "2" = "#FFDA33",   
                   "3" = "#FFA500",   
                   "4" = "#D2691E",   
                   "5" = "#6B4226",   
                   "6" = "#8e44ad",   
                   "7" = "#A9DFBF",   
                   "8" = "#2E8B57")  

# Plot Alcohol Volume over time
plot_vol <- ggplot(data, aes(x = Tag, y = Vol_mean, color = Ansatz, fill = Ansatz)) +
  geom_line(size = 1.5) +  # Increase line size here (change the number to adjust thickness)
  geom_point() +
  geom_ribbon(aes(ymin = Vol_mean - Vol_sd, ymax = Vol_mean + Vol_sd), alpha = 0.2) +
  labs(title = "Alcohol Volume over Time", x = "Time (Days)", y = "Alcohol Volume (%)", color = "Approach", fill = "Approach") +
  scale_color_manual(values = ansatz_colors) +
  scale_fill_manual(values = ansatz_colors) +
  theme_custom

# Display plot
print(plot_vol)
```
## 02 Figure pH all approaches with Ribbons
```r
#Load nessecary Packages 
library(ggplot2)
library(dplyr)
library(gridExtra)

# Load data
data <- read.csv("~/Desktop/Plato_pH_alk.csv", sep=";", dec=",")

# Remove unnecessary columns
data <- data[, c("Ansatz", "Tag", "pH")]

# Group replicates and rename categories
data <- data %>%
  filter(!Ansatz %in% c(7, 8)) %>% # Ignore Ansatz 7 and 8
  mutate(Ansatz = case_when(
    Ansatz %in% c(1, 11, 16) ~ "1",
    Ansatz %in% c(3, 12, 17) ~ "3",
    Ansatz %in% c(5, 13, 18) ~ "5",
    Ansatz %in% c(6, 14, 19) ~ "6",
    Ansatz %in% c(10, 15, 20) ~ "8",
    Ansatz == 9 ~ "7",
    Ansatz == 2 ~ "2",
    Ansatz == 4 ~ "4",
    TRUE ~ as.character(Ansatz)
  )) %>%
  group_by(Ansatz, Tag) %>%
  summarise(
    pH_mean = mean(pH, na.rm=TRUE),
    pH_sd = sd(pH, na.rm=TRUE)
  ) %>%
  ungroup()

# Define custom color palette for the Approaches
ansatz_colors <- c("1" = "#D2B48C",   # Dunkleres Beige für Ansatz 1
                   "2" = "#FFDA33",   # Gelb für Ansatz 2
                   "3" = "#FFA500",   # Orange für Ansatz 3 (new)
                   "4" = "#D2691E",   # Hellbraun für Ansatz 4 (new)
                   "5" = "#6B4226",   # Dunkelbraun für Ansatz 5 (new)
                   "6" = "#8e44ad",   # Dunkellila für Ansatz 6
                   "7" = "#A9DFBF",   # Hellgrün für Ansatz 7
                   "8" = "#2E8B57")   # Dunkelgrün für Ansatz 8

# Define custom theme 
theme_custom <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),  
    axis.title = element_text(size = 30, face = "bold"),  # Set axis title size
    axis.text = element_text(size = 24),  # Set axis text size
    strip.text = element_text(size = 30, face = "bold", color = "black"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20),
    panel.grid.major = element_line(size = 0.5, color = "gray"),
    panel.grid.minor = element_line(size = 0.2, color = "gray"),
    plot.title = element_text(size = 27, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 20),  # Adjust X axis title size
    axis.title.y = element_text(size = 20)   # Adjust Y axis title size
  )

# Plot pH over time with ribbons
plot_pH_ribbon <- ggplot(data, aes(x = Tag, y = pH_mean, color = Ansatz, fill = Ansatz)) +
  geom_line(size = 3) +  # Set line size here (adjust the number to change thickness)
  geom_point() +
  geom_ribbon(aes(ymin = pH_mean - pH_sd, ymax = pH_mean + pH_sd), alpha = 0.2) +
  labs(title = "pH over Time", x = "Time (Days)", y = "pH", color = "Approach", fill = "Approach") +
  scale_color_manual(values = ansatz_colors) +
  scale_fill_manual(values = ansatz_colors) +
  theme_custom

# Display plot
print(plot_pH_ribbon)
```
## 03 Plot all approaches individually pH, Vol%(°P), Vol%(enz.)
```r
#Load nessecary Packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(zoo)

# Load data
data <- read.csv("~/Desktop/Plato_pH_alk.csv", sep=";", dec=",")
endalkohol <- read.csv("~/Desktop/Endalkoholneuneu.csv", sep=";", dec=",")

# Aggregate endalcohol data
endalkohol_summary <- endalkohol %>%
  group_by(Ansatz) %>%
  summarise(
    Ethanol_mean = mean(Ethanol..Vol., na.rm = TRUE),
    Ethanol_sd = sd(Ethanol..Vol., na.rm = TRUE))

# Remove unnecessary columns
data <- data[, c("Ansatz", "Tag", "pH", "Plato_korr", "Vol_prozent")]

# group Approaches
data <- data %>%
  filter(!Ansatz %in% c(7, 8)) %>%
  mutate(Ansatz = case_when(
    Ansatz %in% c(1, 11, 16) ~ "1",
    Ansatz %in% c(3, 12, 17) ~ "3",
    Ansatz %in% c(5, 13, 18) ~ "5",
    Ansatz %in% c(6, 14, 19) ~ "6",
    Ansatz %in% c(10, 15, 20) ~ "8",
    Ansatz == 9 ~ "7",
    TRUE ~ as.character(Ansatz)
  )) %>%
  group_by(Ansatz, Tag) %>%
  summarise(
    pH_mean = mean(pH, na.rm=TRUE),
    pH_sd = sd(pH, na.rm=TRUE),
    Vol_mean = mean(Vol_prozent, na.rm=TRUE),
    Vol_sd = sd(Vol_prozent, na.rm=TRUE)
  ) %>%
  ungroup()

# Interpolate missing values
data <- data %>%
  group_by(Ansatz) %>%
  mutate(Vol_mean = ifelse(is.na(Vol_mean), zoo::na.approx(Vol_mean, na.rm = FALSE), Vol_mean)) %>%
  ungroup()

# Define colour palette
farben <- c("1" = "red", "3" = "green", "4" = "blue", "5" = "orchid1",
            "6" = "purple4", "7" = "beige", "8" = "turquoise")

# Axis limits
y_limits <- c(-0.5, 7)
x_breaks <- seq(1, max(data$Tag), by = 2)

# define Theme
theme_custom <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "#EAEAEA", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 22, face = "bold", color = "black"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_line(size = 0.5, color = "gray"),
    panel.grid.minor = element_line(size = 0.2, color = "gray"),
    plot.title = element_text(size = 27, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0.5, size = 18),
    axis.title.x = element_text(margin = margin(t = 10), size = 22, face = "bold"),
    axis.text.x.top = element_text(angle = -45, hjust = 1, vjust = 1, size = 22, face = "bold")
  )

# Create the plot
plot_ph_vol <- function(ansatz) {
  df_plot <- filter(data, Ansatz == ansatz)
  df_end <- filter(endalkohol_summary, Ansatz == as.numeric(ansatz))
  
  ggplot(df_plot, aes(x = Tag)) +
    geom_line(aes(y = pH_mean, linetype = "pH"), color = "blue", na.rm = TRUE) +
    geom_point(aes(y = pH_mean), color = "blue") +
    geom_errorbar(aes(y = pH_mean, ymin = pH_mean - pH_sd, ymax = pH_mean + pH_sd), width = 0.2, color = "blue") +
    
    geom_line(aes(y = Vol_mean, linetype = "Vol%"), color = "orange", na.rm = TRUE) +
    geom_point(aes(y = Vol_mean), color = "orange") +
    geom_ribbon(aes(ymin = Vol_mean - Vol_sd, ymax = Vol_mean + Vol_sd), alpha = 0.2, fill = "orange") +
    
    geom_point(data = df_end, aes(x = max(df_plot$Tag) + 1, y = Ethanol_mean), color = "orange", size = 4) +
    geom_errorbar(data = df_end, aes(x = max(df_plot$Tag) + 1, ymin = Ethanol_mean - Ethanol_sd, ymax = Ethanol_mean + Ethanol_sd), color = "orange", width = 0.2) +
    
    scale_y_continuous(name = "pH", limits = y_limits, sec.axis = sec_axis(~ ., name = "Vol%")) +
    scale_x_continuous(breaks = c(x_breaks, max(df_plot$Tag) + 1), labels = c(x_breaks, "Vol%(enz.)")) +
    scale_linetype_manual(values = c("pH" = "dotted", "Vol%" = "solid")) +
    labs(title = paste("Approach", ansatz), x = "Time (Days)", linetype = "Parameter") +
    theme_custom
}

# List of approaches
ansaetze <- unique(data$Ansatz)
ph_vol_plots <- lapply(ansaetze, plot_ph_vol)

# Show plot
grid.arrange(grobs = ph_vol_plots, ncol = 2, top = textGrob("pH and Vol% over Time per Approach", gp=gpar(fontsize=27, fontface="bold")))
```
## 04 Plot of the cell concentration of both approaches
```r
# Load data
   
file_path3 <- "~/Desktop/Zellkonzentrationsbestimmung.csv"
data_Zellkonz <- read.csv(file_path3,sep = ";",dec = ",",na.strings = "")

#Filter data 
   
filtered_data <- data_Zellkonz %>% filter(Ansatz %in% c(14, 9))
   
#Group the data by day and approach and calculate the mean and SD
    
summary_data <- filtered_data %>%
  group_by(Tag, Ansatz) %>%
  summarise(Mittelwert = mean(Zellkonz, na.rm = TRUE), 
            SD = sd(Zellkonz, na.rm = TRUE), 
            .groups = "drop")

# Create Plot 
Plot_zellkonz_ribbon <- ggplot(summary_data, aes(x = Tag, y = Mittelwert, 
                                                 color = as.factor(Ansatz), 
                                                 fill = as.factor(Ansatz), 
                                                 group = Ansatz)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = Mittelwert - SD, ymax = Mittelwert + SD), alpha = 0.3, color = NA) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = seq(min(summary_data$Tag), max(summary_data$Tag), by = 2)) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("9", "14")) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("9", "14")) +
  labs(title = "Cell Concentration Over Time", x = "Day", y = "Abundance (Mean ± SD)", 
       color = "Approach", fill = "Approach") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
     
# Print Plot
print(Plot_zellkonz_ribbon)
View(plot)
```
# Figures enzymatic tests/CO2
Here we show all the codes that were used to create the images of the results of the enzymatic tests and CO2 experiment 

## 01 Figure final alcohol content determined enzymatically
```r
# Load libraries
library(ggplot2)
library(dplyr)

# Load data
data <- read.csv("~/Desktop/Endalkneu.csv", sep=";", dec=",")

# Select relevant columns and define groups
data <- data[, c("Ansatz", "Ethanol..Vol.")] %>%
  mutate(Ansatz = case_when(
    Ansatz %in% c(10, 110, 160) ~ "1", Ansatz %in% c(30, 120, 170) ~ "3",
    Ansatz %in% c(50, 130, 180) ~ "5", Ansatz %in% c(60, 140, 190) ~ "6",
    Ansatz %in% c(100, 150, 200) ~ "8", Ansatz == 90 ~ "7",
    Ansatz == 40 ~ "4", Ansatz == 20 ~ "2", TRUE ~ NA_character_)) %>%
  filter(!is.na(Ansatz)) %>% group_by(Ansatz) %>%
  summarise(Ethanol_mean = mean(Ethanol..Vol., na.rm = TRUE),
            Ethanol_sd = ifelse(n() > 1, sd(Ethanol..Vol., na.rm = TRUE), NA)) %>%
  ungroup()

# Define groups and order
data$Group <- factor(case_when(
  data$Ansatz %in% c("1", "2", "3", "4", "5") ~ "Beer",
  data$Ansatz == "6" ~ "Grape Juice", data$Ansatz %in% c("7", "8") ~ "YPD"), 
  levels = c("YPD", "Grape Juice", "Beer"))

# Define colors
ansatz_colors <- c("1" = "#F5DEB3", "2" = "#FFDA33", "3" = "#FFA500", "4" = "#D2B48C",
                   "5" = "#8B4513", "6" = "#8e44ad", "7" = "#A9DFBF", "8" = "#27AE60")

# Create bar plot
Plot_Endalk <- ggplot(data, aes(x = Ethanol_mean, y = Ansatz, fill = Ansatz)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_errorbar(aes(xmin = Ethanol_mean - Ethanol_sd, xmax = Ethanol_mean + Ethanol_sd), 
                width = 0.2, color = "black") +
  labs(title = "Ethanol Concentration in Different Approaches", x = "Ethanol (Vol%)", y = "Approach") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"), 
        legend.position = "none") +
  scale_fill_manual(values = ansatz_colors)

print(Plot_Endalk)
```
## 02 Plot CO2 concentration
```r
# Load library
library(ggplot2)

# Import data
data <- read.csv("~/Desktop/CO2titration.csv", sep = ";", dec = ",", na.strings = "")

# Convert columns
data$CO2_g_L <- as.numeric(gsub(",", ".", data$CO2.g.L.))
data$Approach <- as.numeric(data$Approach)

# Define approach order and filter data
ordered_approaches <- c(11, 16, 14, 19, 13, 18, 12, 17, 15, 20)
data <- data[data$Approach %in% ordered_approaches, ]
data$Approach <- factor(data$Approach, levels = ordered_approaches)
data$group <- rep(1:(length(ordered_approaches) / 2), each = 2)

# Define colors
ansatz_colors <- c("11" = "#D2B48C", "16" = "#D2B48C", "14" = "#8e44ad", "19" = "#8e44ad",
                   "13" = "#6B4226", "18" = "#6B4226", "12" = "#FFA500", "17" = "#FFA500",
                   "15" = "#2E8B57", "20" = "#2E8B57")

# Create bar plot
plot <- ggplot(data, aes(x = factor(group), y = CO2_g_L, fill = factor(Approach))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
  scale_fill_manual(values = ansatz_colors, name = "Approach") +
  labs(title = "CO₂ Concentration", x = "Approaches", y = "CO₂ (g/L)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 12, face = "bold")) +
  scale_x_discrete(labels = c("Approach 11 & 16", "Approach 14 & 19", "Approach 13 & 18",
                              "Approach 12 & 17", "Approach 15 & 20")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Print plot
print(plot)
```

## 03 Plot Sugar concentration with standard deviation
```r
# Load libraries
library(ggplot2); library(dplyr); library(tidyr)

# Import data
data <- read.csv("~/Desktop/Zuckerkonz_sd_neu.csv", sep = ";", dec = ",")

# Reshape data (long format for concentration & standard deviation)
data_long <- data %>% pivot_longer(cols = starts_with("Konz"), names_to = "Sugar", values_to = "Concentration")
sd_long <- data %>% pivot_longer(cols = starts_with("sd_"), names_to = "SD_Type", values_to = "SD") %>%
  mutate(SD_Type = recode(SD_Type, "sd_Glucose" = "Konz_Glucose_g.l", "sd_Saccharose" = "Konz_Saccharose_g.l", "sd_Maltose" = "Konz_Maltose_g.l"))

# Merge concentration & SD, replace NA with 0
data_long <- data_long %>% left_join(sd_long, by = c("Ansatz", "Sugar" = "SD_Type")) %>%
  mutate(SD = replace_na(SD, 0), Sugar = recode(Sugar, "Konz_Glucose_g.l" = "Glucose", 
         "Konz_Saccharose_g.l" = "Saccharose", "Konz_Maltose_g.l" = "Maltose"))

# Define fermentation stages & approaches
data_long <- data_long %>% mutate(Group = case_when(
  Ansatz %in% c("1.1", "3.1", "5.1", "6.1", "9.1", "10.1") ~ "before fermentation",
  Ansatz %in% c("1.2", "3.2", "5.2", "6.2", "10.2") ~ "after fermentation (Mean)",
  Ansatz == "9.2" ~ "after fermentation"),
  Approach = recode(Ansatz, "1.1" = "Approach 1", "1.2" = "Approach 1", "3.1" = "Approach 3", "3.2" = "Approach 3",
                    "5.1" = "Approach 5", "5.2" = "Approach 5", "6.1" = "Approach 6", "6.2" = "Approach 6",
                    "9.1" = "Approach 7", "9.2" = "Approach 7", "10.1" = "Approach 8", "10.2" = "Approach 8"))

# Factor levels for ordering
data_long$Approach <- factor(data_long$Approach, levels = c("Approach 1", "Approach 3", "Approach 5", "Approach 6", "Approach 7", "Approach 8"))
data_long$Group <- factor(data_long$Group, levels = c("before fermentation", "after fermentation (Mean)", "after fermentation"))

# Adjust SD range
data_long <- data_long %>% mutate(ymin = pmax(Concentration - SD, 0), ymax = Concentration + SD)

# Define colors
farben <- c("before fermentation" = "#1f77b4", "after fermentation (Mean)" = "#ff7f0e", "after fermentation" = "#ff7f0e")

# Create plot
plot <- ggplot(data_long, aes(x = Sugar, y = Concentration, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8, preserve = "single"), color = "black") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width = 0.8), width = 0.2, color = "black") +
  facet_wrap(~Approach, nrow = 2) + scale_fill_manual(values = farben, name = "Fermentation Stage") +
  labs(title = "Sugar Concentrations Before and After Fermentation", x = "Sugar Type", y = "Concentration (g/L)") +
  theme_minimal() + theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18), strip.text = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))

# Print plot
print(plot)
```

# Figure alcohol tolerance 
Here we show the code used to create the Figure of the alcoholtolerance experiment 
## 01 alcohol tolerance 
```r
# Load libraries
library(ggplot2); library(dplyr)

# Import data
data_alkoholtoleranz <- read.csv("~/Desktop/Alkoholtoleranz.csv", sep = ";", dec = ",", na.strings = "")

# Convert comma to dot for numerical values
data_alkoholtoleranz$Alkkonz <- as.numeric(gsub(",", ".", data_alkoholtoleranz$Alkkonz))
data_alkoholtoleranz$Coverage <- as.numeric(gsub(",", ".", data_alkoholtoleranz$Coverage..))

# Calculate mean and error margins
data_summary <- data_alkoholtoleranz %>%
  group_by(Organismus, Alkkonz, `Zeit..h.`) %>%
  summarise(
    Mean_Coverage = mean(Coverage, na.rm = TRUE),
    SD_Coverage = sd(Coverage, na.rm = TRUE),
    SE_Coverage = SD_Coverage / sqrt(n()),
    .groups = 'drop')

# Define colors
farben <- c("Prototheca 24h" = "#6baed6", "Prototheca 72h" = "#1f77b4", 
            "S.cerevisiae 24h" = "#66c2a5", "S.cerevisiae 72h" = "#2ca02c")

# Create plot
plot_alcohol_tolerance <- ggplot(data_summary, aes(
  x = Alkkonz, y = Mean_Coverage, 
  color = interaction(Organismus, `Zeit..h.`), 
  fill = interaction(Organismus, `Zeit..h.`),
  group = interaction(Organismus, `Zeit..h.`))) +
  geom_ribbon(aes(ymin = Mean_Coverage - SE_Coverage, ymax = Mean_Coverage + SE_Coverage), 
              alpha = 0.2, linetype = "blank") +
  geom_line(linewidth = 3) +
  geom_point(size = 3) +
  labs(
    title = "Alcohol Tolerance",
    x = "Alcohol Concentration (%)",
    y = "Coverage (%)",
    color = "Organism and Time",
    fill = "Organism and Time"
  ) +
  scale_color_manual(values = farben) +
  scale_fill_manual(values = farben) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    legend.position = "top",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14, face = "italic"),
    panel.grid.major = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.minor = element_line(color = "gray90", linetype = "dotted")
  ) +
  scale_x_continuous(breaks = seq(0, max(data_summary$Alkkonz, na.rm = TRUE), by = 1)) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), expand = expansion(mult = c(0, 0.05)))

# Print plot
print(plot_alcohol_tolerance)
```
