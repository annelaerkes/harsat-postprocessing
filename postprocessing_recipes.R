# HARSAT post-processing recipes ----
#
# One saved HARSAT assessment, several outputs. Run from top to bottom in the
# project folder (harsat-postprocessing.Rproj). HARSAT itself is not
# needed, only the saved .rds file.

library(dplyr)
library(ggplot2)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figs/all_series", recursive = TRUE, showWarnings = FALSE)


## Change these for your own data ----
my_file      <- "data/PAH.rds"
my_series    <- "FJBA PYR Mytilus edulis SB"   # one name from names(res$assessment)
my_y_label   <- "Pyrene, ug/kg dry weight"
my_substance <- "PYR"
my_species   <- "Mytilus edulis"


## Step 1: open the object ----
res <- readRDS(my_file)

names(res)                 # the parts of the object
head(res$timeSeries)       # one row per time series
names(res$assessment)      # the name of each time series

series <- res$assessment[[my_series]]   # the results for one time series
names(series)


## Step 2: build one time series plot ----
means   <- series$annualIndex    # annual means, on the log scale
trend   <- series$pred           # fitted trend and confidence band, on the log scale
samples <- series$data           # individual samples, on the original scale
limit   <- series$summary$EQS    # threshold, on the original scale

# exp() brings log-scale values back to the original scale. Biological variables
# are usually not on the log scale: for those, remove exp().

# 2a: the annual means
plot_1 <- ggplot() +
  geom_point(data = means, aes(year, exp(index)), size = 2)
plot_1

# 2b: add the fitted trend and its confidence band
plot_2 <- plot_1 +
  geom_ribbon(data = trend, aes(year, ymin = exp(ci.lower), ymax = exp(ci.upper)), alpha = 0.25) +
  geom_line(data = trend, aes(year, exp(fit)))
plot_2

# 2c: add the individual samples
plot_3 <- plot_2 +
  geom_point(data = samples, aes(year, concentration), size = 0.75)
plot_3

# 2d: add the threshold and labels
plot_4 <- plot_3 +
  geom_hline(yintercept = limit, linetype = "dashed") +
  labs(title = my_series, x = NULL, y = my_y_label) +
  theme_classic(base_size = 15)
plot_4

ggsave("output/figs/time_series.png", plot_4, width = 6, height = 6, dpi = 300)


## Step 3: make a summary table ----
# Each time series has its own one-row summary: put them all in one table
summaries <- res$assessment |>
  lapply(function(x) x$summary) |>
  bind_rows(.id = "series")

# Add what each series is (substance, species, unit) and where it is (station)
all_series <- res$timeSeries |>
  mutate(series = rownames(res$timeSeries)) |>
  left_join(res$stations, by = "station_code") |>
  left_join(summaries, by = "series")

write.csv(all_series, "output/tables/all_series.csv", row.names = FALSE, na = "")

# One substance and species, with the columns you want
one_substance <- filter(all_series, determinand == my_substance, species == my_species)

my_table <- one_substance |>
  arrange(desc(station_latitude)) |>
  transmute(Station = station_name,
            Years = paste(firstYearAll, lastyear, sep = "-"),
            `Recent change, %/year` = recent_change,
            `p recent` = round(p_recent_change, 3),
            `Mean last year` = round(meanLY, 2),
            `Upper CL` = round(clLY, 2),
            Threshold = round(EQS, 1),
            Status = EQSbelow)
my_table

write.csv(my_table, "output/tables/one_substance.csv", row.names = FALSE, na = "")


## Step 4: make a map ----
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# The base map: change xlim (longitude) and ylim (latitude) for your own region
base_map <- ggplot(world) +
  geom_sf(colour = "white", fill = "grey") +
  coord_sf(xlim = c(8, 25), ylim = c(54, 70), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal()
base_map

# Fitted mean in the last year. Red outline: HARSAT's status is "above" the
# threshold, which it bases on the upper confidence limit, not on the mean.
map_mean <- base_map +
  geom_point(data = one_substance,
             aes(station_longitude, station_latitude, fill = meanLY, colour = EQSbelow),
             shape = 21, size = 5, stroke = 1) +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +   # dark = high
  scale_colour_manual(values = c(above = "red", below = "black"), na.value = "black") +
  labs(fill = "Mean", colour = "Status")
map_mean

ggsave("output/figs/map_mean.png", map_mean, width = 6, height = 6, dpi = 300)

# Another column gives another map: here the recent change
map_change <- base_map +
  geom_point(data = one_substance,
             aes(station_longitude, station_latitude, fill = recent_change),
             shape = 21, size = 5, stroke = 1) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(fill = "% per year")
map_change

ggsave("output/figs/map_change.png", map_change, width = 6, height = 6, dpi = 300)


## Step 5: the same plot for many time series (optional) ----
# Step 2 written as a function, so it can be repeated. Some series have no
# fitted trend or no threshold, so those parts are only added if they exist.
plot_one_series <- function(id) {
  series <- res$assessment[[id]]
  info <- res$timeSeries[id, ]

  p <- ggplot() +
    geom_point(data = series$annualIndex, aes(year, exp(index)), size = 2) +
    geom_point(data = series$data, aes(year, concentration), size = 0.75) +
    labs(title = id, x = NULL, y = paste(info$determinand, info$unit, info$basis)) +
    theme_classic(base_size = 15)

  if (!is.null(series$pred)) {
    p <- p +
      geom_ribbon(data = series$pred, aes(year, ymin = exp(ci.lower), ymax = exp(ci.upper)), alpha = 0.25) +
      geom_line(data = series$pred, aes(year, exp(fit)))
  }
  if (!is.na(series$summary$EQS)) {
    p <- p + geom_hline(yintercept = series$summary$EQS, linetype = "dashed")
  }
  p
}

for (id in one_substance$series) {
  ggsave(file.path("output/figs/all_series", paste0(id, ".png")),
         plot_one_series(id), width = 6, height = 6, dpi = 300)
}


## Want more? ----
# R/harsat_recipes.R has the full report-style figures, which also show values
# below LOQ and years left out of the trend fit. These flat plots do not.
# source("R/harsat_recipes.R")
# plot_series(res, my_series)
# plot_map(harsat_summary(res) |> filter(determinand == my_substance, species == my_species), "mean")
