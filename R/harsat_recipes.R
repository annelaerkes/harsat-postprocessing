# Functions for post-processing a saved HARSAT assessment (.rds) ----
#
# Simplified from the Swedish reporting pipeline
# (https://github.com/NRM-MOC/fmcom-annual-report: R/harsat_helpers.R and
# R/harsat_tables_figures.R). Everything is read from the HARSAT object itself,
# so there are no config or lookup files, and HARSAT does not need to be
# installed. These are the optional take-home versions of the plain code in
# postprocessing_recipes.R, with the full report-style figures.
#
# Needs: dplyr, ggplot2, patchwork, sf, rnaturalearth, rnaturalearthdata

library(dplyr)
library(ggplot2)
library(patchwork)


## Summary table ----
# One row per time series: trends, fitted mean and thresholds from
# res$assessment[[series]]$summary, with names from the rest of the object.
# Threshold columns keep HARSAT's names (e.g. EQS, EQSdiff, EQSachieved, EQSbelow).
harsat_summary <- function(res) {
  # Each series has its own one-row summary (NULL if it was not assessed)
  results <- lapply(names(res$assessment), function(id) {
    s <- res$assessment[[id]]$summary
    if (!is.null(s)) cbind(series = id, s)
  }) |>
    bind_rows() |>
    rename(n_year_all = nyall, n_year_fit = nyfit, n_year_positive = nypos,
           first_year_all = firstYearAll, first_year_fit = firstYearFit,
           last_year = lastyear, detectable_trend = dtrend,
           mean_last_year = meanLY, climit_last_year = clLY)

  # What each series is: station, substance, species, tissue, unit
  series_info <- res$timeSeries |>
    mutate(series = rownames(res$timeSeries),
           determinand_name = res$info$determinand[determinand, "common_name"],
           distribution = res$info$determinand[determinand, "distribution"],
           species_name = res$info$species[species, "common_name"]) |>
    left_join(select(res$stations, station_code, station_name,
                     station_latitude, station_longitude),
              by = "station_code") |>
    select(series, station_code, station_name, station_latitude, station_longitude,
           determinand, determinand_name, species, species_name, matrix, basis, unit,
           distribution)

  left_join(series_info, results, by = "series")
}


## Time series plot ----
# `id` is a series name, e.g. "FJBA PYR Mytilus edulis SB" (see names(res$assessment)).
#   small points: individual samples; large points: annual means
#   open blue:    below LOQ
#   grey:         years HARSAT left out of the trend fit
#   dashed line and green area: the threshold (`threshold`, e.g. "EQS" or "EAC")
plot_series <- function(res, id, threshold = "EQS", min_years = 10) {
  a <- res$assessment[[id]]
  if (is.null(a)) {
    stop("No assessment for '", id, "'. The series are in names(res$assessment).", call. = FALSE)
  }
  obs <- res$timeSeries[id, ]
  det <- res$info$determinand[obs$determinand, ]
  last_year <- res$info$max_year

  # HARSAT models lognormal series (contaminants) on the log scale and normal
  # series (e.g. biological variables) on the original scale. The fitted
  # values are transformed back to the scale of the data.
  back <- if (det$distribution %in% "normal") identity else exp

  # The colour and shape of a point tell how HARSAT treated that value
  style_points <- function(d) {
    fit_years <- unique(a$data$year)
    below_loq <- d$censoring %in% "<"
    mutate(d,
           colour = case_when(!year %in% fit_years ~ "grey65",
                              below_loq ~ "#0072B2",
                              TRUE ~ "black"),
           shape = if_else(below_loq, 21, 19))
  }
  # All samples are in fullData if HARSAT dropped some years, otherwise in data
  samples <- if (is.null(a$fullData)) a$data else a$fullData
  samples <- style_points(samples)
  means <- style_points(mutate(a$annualIndex, mean = back(index)))

  p <- ggplot()

  # Threshold, if this series has one
  refs <- a$reference.values
  limit <- if (!is.null(refs) && threshold %in% rownames(refs)) {
    back(refs[threshold, "refvalue"])
  } else {
    NA
  }
  if (!is.na(limit)) {
    p <- p +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = limit,
               fill = "lightgreen", alpha = 0.3) +
      geom_hline(yintercept = limit, linetype = "dashed")
  }

  # Fitted trend and confidence band, if a model was fitted (a$pred)
  if (!is.null(a$pred)) {
    p <- p +
      geom_ribbon(data = a$pred, aes(x = year, ymin = back(ci.lower), ymax = back(ci.upper)),
                  fill = "black", alpha = 0.25) +
      geom_line(data = a$pred, aes(x = year, y = back(fit)), alpha = 0.1)
  }

  # Axis text: biological variables have no basis or tissue
  basis <- c(W = "ww", D = "dw", L = "lw")[obs$basis]
  tissue <- c(MU = "muscle", LI = "liver", SB = "soft body", EH = "egg content")[obs$matrix]
  y_label <- if (det$biota_group %in% "Biological") {
    paste0(det$common_name, ", ", obs$unit)
  } else {
    paste(det$common_name, obs$unit, basis, tissue)
  }
  station <- res$stations$station_name[match(obs$station_code, res$stations$station_code)]
  species <- res$info$species[obs$species, "common_name"]

  p +
    geom_point(data = samples, aes(year, concentration, colour = colour, shape = shape),
               size = 0.75, fill = "white") +
    geom_point(data = means, aes(year, mean, colour = colour, shape = shape),
               size = 2, fill = "white") +
    scale_colour_identity() +
    scale_shape_identity() +
    labs(title = paste(station, species, sep = " - "), x = NULL, y = y_label) +
    theme_classic(base_size = 15) +
    theme(axis.text = element_text(colour = "black")) +
    # Whole years only on the x-axis
    scale_x_continuous(breaks = function(limits) {
      breaks <- pretty(limits)
      breaks[breaks == round(breaks)]
    }) +
    # coord_cartesian() crops the view: y-axis to 1.5 x the highest annual mean,
    # x-axis to the last year and at least `min_years` years. Limits set with a
    # scale instead would remove the parts of the band outside the axis.
    coord_cartesian(xlim = c(min(samples$year, last_year - min_years + 1), last_year),
                    ylim = c(0, max(means$mean * 1.5, na.rm = TRUE)))
}


## Map and bar chart ----
# `data`: rows of harsat_summary() with one series per station, i.e. filtered
# to one substance and species.
# what = "mean":   fitted mean in the last year. Red outline: HARSAT's status is
#                  "above" the threshold, which it bases on the upper confidence
#                  limit (climit_last_year), not on the mean itself
# what = "change": recent yearly change (blue: decrease, red: increase)
plot_map <- function(data, what = c("mean", "change"), threshold = "EQS",
                     xlim = c(8, 25), ylim = c(54, 70)) {
  what <- match.arg(what)
  if (anyDuplicated(data$station_code)) {
    stop("Several series per station. Filter to one substance, species and tissue.", call. = FALSE)
  }
  data <- mutate(data, station_name = reorder(station_name, station_latitude))

  if (what == "mean") {
    data <- mutate(data, value = mean_last_year,
                   outline = if_else(.data[[paste0(threshold, "below")]] %in% "above", "red", "black"))
    fill_scale <- scale_fill_continuous(low = "#56B1F7", high = "#132B43")
    value_label <- paste0("fitted mean ", max(data$last_year), ", ", data$unit[1])
  } else {
    data <- mutate(data, value = recent_change, outline = "black")
    fill_scale <- scale_fill_gradient2(low = scales::muted("blue"), mid = "white",
                                       high = scales::muted("red"))
    value_label <- change_label(data$distribution[1], data$unit[1])
  }

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  map <- ggplot(world) +
    geom_sf(colour = "white", fill = "grey") +
    geom_point(data = data, aes(station_longitude, station_latitude, fill = value, colour = outline),
               shape = 21, size = 5, stroke = 1, alpha = 0.8) +
    scale_colour_identity() +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = "Longitude", y = "Latitude")

  bars <- ggplot(data, aes(value, station_name, fill = value)) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.5) +
    labs(x = value_label, y = NULL)

  (map + bars) +
    plot_annotation(title = paste(data$determinand_name[1], "in", tolower(data$species_name[1]))) &
    fill_scale &
    theme_minimal() &
    theme(legend.position = "none")
}

# Trends are in % per year for lognormal series (contaminants) and in the
# variable's unit per year for normal series (biological variables)
change_label <- function(distribution, unit) {
  if (distribution %in% "normal") {
    paste0("yearly change, ", if (unit == "%") "percentage points" else unit, " per year")
  } else {
    "yearly relative change, %"
  }
}
