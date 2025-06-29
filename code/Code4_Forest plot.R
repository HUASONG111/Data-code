#install.packages(c("readxl", "dplyr", "ggplot2", "patchwork"))
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)

## Please change the following paths to match your own file locations.
data_dir <- "/Users/huasong/Desktop/data/4.Forest meta data"

neg <- read_excel(file.path(data_dir, "Forest_meta_neg.xlsx"))
pos <- read_excel(file.path(data_dir, "Forest_meta_pos.xlsx"))

neg <- neg %>%
  rename(RR = `Minimum RR`, Lag = `Lags to min`) %>%
  mutate(Type = "Minimum CRR")

pos <- pos %>%
  rename(RR = `Maximum RR`, Lag = `Lags to max`) %>%
  mutate(Type = "Maximum CRR")


data <- bind_rows(pos, neg)

data$Biome <- factor(data$Biome, levels = c("Zonefour", "Zonethree", "Zonetwo", "Zoneone"))
data$Factors <- factor(data$Factors, levels = c("Phosphate", "Nitrate", "SST", "SSW"))
data$Value <- factor(data$Value, levels = c("25th", "75th"))

color_map <- c(
  "Minimum CRR" = "#4DBBD5FF",
  "Maximum CRR" = "#E64B35FF"
)
shape_map <- c(
  "Zoneone" = 16,
  "Zonetwo" = 15,
  "Zonethree" = 17,
  "Zonefour" = 18
)

plot_percentile <- function(percentile, show_y_axis = TRUE, show_strip = TRUE) {
  data_subset <- data %>% filter(Value == percentile)
  
  if (percentile %in% c("25th", "75th")) {
    x_limits <- c(0.4, 1.6)
    x_breaks <- seq(0.4, 1.6, by = 0.2)
  } else {
    xmin <- min(data_subset$LowerCI, na.rm = TRUE)
    xmax <- max(data_subset$UpperCI, na.rm = TRUE)
    padding <- (xmax - xmin) * 0.1
    x_limits <- c(xmin - padding, xmax + padding)
    x_breaks <- seq(floor(x_limits[1] * 10) / 10, ceiling(x_limits[2] * 10) / 10, by = 0.2)
  }
  
  p <- ggplot(data_subset, aes(x = RR, y = Biome, color = Type, shape = Biome)) +
    geom_point(position = position_dodge(width = 0), size = 4) +
    # geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI),
    #                position = position_dodge(width = 0), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    facet_grid(Factors ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = color_map) +
    scale_shape_manual(values = shape_map) +
    scale_x_continuous(breaks = seq(0.5, 3.5, by = 0.5)) +
    coord_cartesian(xlim = c(0, 3.5)) +
    labs(
      title = if (percentile == "25th") bquote(25^th) else bquote(75^th),
      x = "CRR", y = NULL, color = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = if (show_strip) element_text(angle = 0, face = "bold", size = 12) else element_blank(),
      axis.text.y = if (show_y_axis) element_text(size = 10) else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank(),
      axis.text.x = element_text(size = 11),
      axis.line.x = element_line(color = "black"),
      axis.ticks.x = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.justification = "top",
      legend.box = "vertical",
      legend.box.margin = margin(0, 10, 0, 0, unit = "pt"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  return(p)
}


p25 <- plot_percentile("25th", show_y_axis = TRUE, show_strip = TRUE)
p75 <- plot_percentile("75th", show_y_axis = FALSE, show_strip = FALSE)

print(p25 + p75 + plot_layout(guides = "collect"))

