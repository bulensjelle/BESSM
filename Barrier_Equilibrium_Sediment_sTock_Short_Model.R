# ================================
# 0. LOAD LIBRARIES
# ================================
library(terra)
library(dplyr)
library(purrr)
library(ggplot2)

set.seed(123)

# ================================
# 1. INPUT PATHS
# ================================
dem_paths <- list(
  "2016" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2016.tif",
  "2022" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Rasters/Topography_2022.tif"
)

polygon_paths <- list(
  "Terschelling" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Terschelling2.shp",
  "Ameland"      = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Ameland2.shp",
  "Texel"        = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Texel2.shp",
  "Vlieland"     = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Vlieland2.shp",
  "Schiermonnikoog" = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Schiermonnikoog2.shp",
  "Borkum"       = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Borkum2.shp",
  "Juist"        = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Juist2.shp",
  "Norderney"    = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Nordeney2.shp",
  "Baltrum"      = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Baltrum2.shp",
  "Langeoog"     = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Langeoog2.shp",
  "Spiekeroog"   = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Spiekeroog2.shp",
  "Wangerooge"   = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Wangeroog2.shp",
  "Amrum"        = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Amrum2.shp",
  "Sylt"         = "C:/Users/6679242/OneDrive - Universiteit Utrecht/Documents/WADWAD/Data/Geodata/Accretion_vs_SLR/Islands/Sylt2.shp"
)

# ================================
# 2. LOAD DATA
# ================================
polygons <- map(polygon_paths, vect)
dems <- map(dem_paths, rast)

# Reproject polygons to DEM CRS
polygons <- map(polygons, ~project(.x, crs(dems[[1]])))

# ================================
# 3. GEO CODES AND COLORS
# ================================
geo_levels <- c("1", "2", "3", "4")
geo_labels <- c("Beach", "Dune Complex", "Polder", "Back-barrier Marshes")

geo_colors <- c(
  "Beach" = "#FFCC33",              
  "Dune Complex" = "#FFBB03",       
  "Polder" = "#39b54a",             
  "Back-barrier Marshes" = "#99CC33"
)

# ================================
# 4. CALCULATE GEO STATS (DEP/EROSION + CALIBRATION)
# ================================
calc_geo_stats <- function(poly, dem2016, dem2022) {
  
  poly$GEO <- factor(poly$GEO, levels = geo_levels, labels = geo_labels)
  
  vals <- terra::extract(c(dem2016, dem2022), poly)
  names(vals)[2:3] <- c("value_2016", "value_2022")
  vals$GEO <- poly$GEO[vals$ID]
  
  geo_stats <- vals %>%
    group_by(GEO) %>%
    summarise(
      deposition = sum(value_2022[value_2022 > value_2016] - value_2016[value_2022 > value_2016], na.rm = TRUE),
      erosion    = sum(value_2022[value_2022 < value_2016] - value_2016[value_2022 < value_2016], na.rm = TRUE),
      net_change = sum(value_2022 - value_2016, na.rm = TRUE),
      mean_2016 = mean(value_2016, na.rm = TRUE),
      mean_2022 = mean(value_2022, na.rm = TRUE),
      mean_change = mean(value_2022 - value_2016, na.rm = TRUE),
      .groups = "drop"
    )
  
  poly_df <- data.frame(GEO = poly$GEO, area = expanse(poly))
  geo_area <- poly_df %>%
    group_by(GEO) %>%
    summarise(total_area = sum(area), .groups = "drop")
  
  geo_stats <- left_join(geo_stats, geo_area, by = "GEO")
  
  # Compute mean annual flux per m² (6 years between 2016-2022)
  geo_stats <- geo_stats %>% mutate(
    flux_per_m2 = net_change / total_area / 6
  )
  
  return(geo_stats)
}
# ================================ 
# 5. SIMULATION FUNCTION PER ISLAND (STOCHASTIC, CALIBRATED) 
# ================================

simulate_island <- function(island_name, geo_stats, time_steps = 101) {
  
  # Initialize stocks (mean elevation × area)
  stocks <- geo_stats %>% 
    mutate(stock = mean_2016 * total_area) %>%
    select(GEO, stock) %>% 
    deframe()
  
  geo_areas <- geo_stats$total_area
  names(geo_areas) <- geo_stats$GEO
  
  history <- lapply(names(stocks), function(x) numeric())
  names(history) <- names(stocks)
  
  # Parameters
  big_storm_chance <- 0.05
  small_storm_chance <- 0.3
  polder_subsidence_rate <- if (island_name %in% c("Amrum", "Sylt")) {
    0
  } else {
    0.000
  }
  saltmarsh_subsidence_rate <- 0.000
  embankment_threshold <- 4000
  depoldering_value <- 2/3
  polder_original <- stocks["Polder"]
  
  overwash_factor <- 0.007
  
  # Base fluxes (calibrated)
  flux_beach_from_offshore <- 
    geo_stats$flux_per_m2[geo_stats$GEO == "Beach"] * geo_areas["Beach"]
  
  flux_beach_to_dune <- 0.002 * stocks["Beach"]
  
  for (t in 1:time_steps) {
    
    # Stochastic storms
    big_storm <- runif(1) < big_storm_chance
    small_storm <- runif(1) < small_storm_chance
    
    # Mean elevations
    beach_elev <- stocks["Beach"] / geo_areas["Beach"]
    marsh_elev <- stocks["Back-barrier Marshes"] / geo_areas["Back-barrier Marshes"]
    
    # Lagoon → marsh accretion shut-off above 2.5 m
    flux_saltmarsh_from_lagoon <- if (marsh_elev > 2.5) {
      0
    } else {
      geo_stats$flux_per_m2[geo_stats$GEO == "Back-barrier Marshes"] *
        geo_areas["Back-barrier Marshes"]
    }
    
    # Fluxes
    # Storm-driven overwash (area-weighted)
    if (big_storm) {
      
      total_overwash <- overwash_factor * stocks["Dune Complex"]
      
      A_marsh  <- geo_areas["Back-barrier Marshes"]
      A_polder <- geo_areas["Polder"]
      A_total  <- A_marsh + A_polder
      
      frac_marsh  <- A_marsh / A_total
      frac_polder <- A_polder / A_total
      
      dune_to_saltmarsh <- total_overwash * frac_marsh
      dune_to_polder    <- total_overwash * frac_polder
      
    } else {
      dune_to_saltmarsh <- 0
      dune_to_polder    <- 0
    }
    dune_to_beach <- if (small_storm) 0.002 * stocks["Dune Complex"] else 0
    polder_subsidence <- polder_subsidence_rate * stocks["Polder"]
    saltmarsh_subsidence <- saltmarsh_subsidence_rate * stocks["Back-barrier Marshes"]
    beach_erosion <- 0.005 *stocks["Beach"]
    
    # Management
    saltmarsh_to_polder <- if (stocks["Back-barrier Marshes"] > embankment_threshold) {
      2/3 * embankment_threshold
    } else 0
    
    polder_to_saltmarsh <- if (stocks["Polder"] / polder_original < depoldering_value) {
      2/3 * stocks["Polder"]
    } else 0
    
    # Update stocks
    stocks["Beach"] <- stocks["Beach"] +
      flux_beach_from_offshore - flux_beach_to_dune + dune_to_beach - beach_erosion
    
    stocks["Dune Complex"] <- stocks["Dune Complex"] +
      flux_beach_to_dune - dune_to_beach - dune_to_saltmarsh
    
    stocks["Back-barrier Marshes"] <- stocks["Back-barrier Marshes"] +
      dune_to_saltmarsh - dune_to_polder +
      polder_to_saltmarsh - saltmarsh_to_polder +
      flux_saltmarsh_from_lagoon - saltmarsh_subsidence
    
    stocks["Polder"] <- stocks["Polder"] +
      saltmarsh_to_polder + dune_to_polder -
      polder_to_saltmarsh - polder_subsidence
    
    # -------------------------------
    # BEACH NOURISHMENT FROM DUNES
    # -------------------------------
    beach_elev <- stocks["Beach"] / geo_areas["Beach"]
    
    if (beach_elev < 1) {
      
      target_beach_stock <- 1.2 * geo_areas["Beach"]
      needed_sediment <- target_beach_stock - stocks["Beach"]
      
      # Limit nourishment by available dune sediment
      available_dune <- stocks["Dune Complex"]
      transfer <- min(needed_sediment, available_dune)
      
      stocks["Beach"] <- stocks["Beach"] + transfer
      stocks["Dune Complex"] <- stocks["Dune Complex"] - transfer
    }
    
    
    # Record history
    for (k in names(stocks)) {
      history[[k]] <- c(history[[k]], stocks[[k]])
    }
    
    if (saltmarsh_to_polder > 0 || polder_to_saltmarsh > 0) {
      polder_original <- stocks["Polder"]
    }
  }
  
  # Convert to mean elevation
  history_df <- bind_rows(lapply(names(history), function(k) {
    data.frame(
      Years = 1:time_steps,
      Elevation = history[[k]] / geo_areas[k],
      GEO = k
    )
  }))
  
  history_df$island <- island_name
  return(history_df)
}

# ================================
# 6. RUN SIMULATION FOR ALL ISLANDS
# ================================
time_steps <- 278

all_histories <- imap(polygons, function(poly, island_name) {
  geo_stats <- calc_geo_stats(poly, dems[["2016"]], dems[["2022"]])
  simulate_island(island_name, geo_stats, time_steps)
})

all_histories_df <- bind_rows(all_histories)

# ================================
# 7. PLOT MEAN ELEVATION
# ================================
ggplot(all_histories_df, aes(x = Years, y = Elevation, color = GEO)) +
  geom_line(size = 1) +
  facet_wrap(~island, scales = "free_y") +
  labs(title = "Stochastic Sediment Stock Model (Calibrated with DEM data)", x = "Years", y = "Mean NAP Elevation (m)") +
  theme_minimal() +
  scale_color_manual(values = geo_colors) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
