# Load required libraries
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Time steps
time_steps <- 101

# Initial sediment stocks based on Ameland (volume based on AHN4)
stocks <- list(
  "Offshore" = 12000,
  "ETD" = 70,
  "Island Head" = 95,
  "Beach Shore" = 130,
  "Dune Complex" = 1568,
  "Tail" = 327,
  "Saltmarsh" = 422,
  "Tidal Basin" = 500,
  "Polder" = 1397
)

# Areas of the stock based on Ameland, 1 if area does not apply
stock_areas <- list(
  "Offshore" = 1,
  "ETD" = 1,
  "Island Head" = 2721624,
  "Beach Shore" = 2729134.5,
  "Dune Complex" = 12998500.25,
  "Tail" = 4713029,
  "Saltmarsh" = 7702768.5,
  "Tidal Basin" = 1,
  "Polder" = 27450395.25
)

stock_elevation <- list(
  "Offshore" = 1,
  "ETD" = 1,
  "Island Head" = 1,
  "Beach Shore" = 1,
  "Dune Complex" = 1,
  "Tail" = 1,
  "Saltmarsh" = 1,
  "Tidal Basin" = 1,
  "Polder" = 1
)

# Store history for plotting
history <- lapply(names(stocks), function(x) numeric(0))
names(history) <- names(stocks)

history_elevation <- lapply(names(stocks), function(x) numeric(0))
names(history_elevation) <- names(stocks)

# Store polder stock for management
polder_original <- stocks[["Polder"]]

# Beach stability setting
SBeach <- 1  # if 1 -> stable, 0 -> unstable
instability_factor <- 1.1
beach_to_offshore <- 0
beach_unstable_count <- 0

# Parameters
big_storm_chance <- 0.1
small_storm_chance <- 0.3
overwash_multiplier <- 2
polder_subsidence_rate <- 0.003
etd_shoaling_rate <- 0.2
embankment_threshold <- 4000
depoldering_value <- 2/3

# TSL calculation (log-spaced values)
TSL <- 10^seq(0.2, -1, length.out = time_steps)
dSLR <- 0.00

# Simulation loop
for (t in 1:time_steps) {
  if (stocks[["Dune Complex"]] > 0) {
    # Storm surge event
    big_storm <- runif(1) < big_storm_chance
    small_storm <- runif(1) < small_storm_chance
    
    # Flows
    if (stocks[["Beach Shore"]]/4 > beach_to_offshore) {
      offshore_to_etd <- 1
      offshore_to_beach <- 9
      islands_head_to_dune <- 0.01 * stocks[["Island Head"]]
      etd_to_island_head <- 0.03 * stocks[["ETD"]]
      island_head_to_beach <- 0.04 * stocks[["Island Head"]]
      beach_to_dune <- 0.03 * stocks[["Beach Shore"]]
      beach_to_tail <- 0.02 * stocks[["Beach Shore"]]
      beach_to_offshore <- if (SBeach) 0.1 * stocks[["Beach Shore"]] else offshore_to_beach * instability_factor
      dune_to_saltmarsh <- 0.001 * stocks[["Dune Complex"]] * (if (big_storm) overwash_multiplier else 1)
      dune_to_polder <- 0.001 * stocks[["Dune Complex"]] * (if (big_storm) overwash_multiplier else 1)
      dune_to_beach <- stocks[["Dune Complex"]] * (if (small_storm) 0.001 else 0)
      tail_to_offshore <- 0.01 * stocks[["Tail"]]
      tail_to_saltmarsh <- 0.001 * stocks[["Tail"]] * (if (big_storm) overwash_multiplier else 1)
      tidal_to_saltmarsh <- TSL[t]
      tidal_to_polder <- 0 * stocks[["Tidal Basin"]]
      polder_subsidence <- stocks[["Polder"]] * polder_subsidence_rate
    } else {
      beach_unstable_count <- beach_unstable_count + 1
      offshore_to_etd <- 1
      offshore_to_beach <- 9
      islands_head_to_dune <- 0.01 * stocks[["Island Head"]]
      etd_to_island_head <- 0.03 * stocks[["ETD"]]
      island_head_to_beach <- 0.04 * stocks[["Island Head"]]
      beach_to_dune <- 0.03 * stocks[["Beach Shore"]]
      beach_to_tail <- 0.02 * stocks[["Beach Shore"]]
      beach_to_offshore <- if (SBeach) 0.1 * stocks[["Beach Shore"]] else offshore_to_beach * instability_factor
      dune_to_saltmarsh <- 0.001 * stocks[["Dune Complex"]] * (if (big_storm) overwash_multiplier else 1)
      dune_to_polder <- 0.001 * stocks[["Dune Complex"]] * (if (big_storm) overwash_multiplier else 1)
      dune_to_beach <- (stocks[["Dune Complex"]] * (if (small_storm) 0.001 else 0)) + 
        abs(stocks[["Beach Shore"]]/4 - beach_to_offshore) * 4
      tail_to_offshore <- 0.01 * stocks[["Tail"]]
      tail_to_saltmarsh <- 0.001 * stocks[["Tail"]] * (if (big_storm) overwash_multiplier else 1)
      tidal_to_saltmarsh <- TSL[t]
      tidal_to_polder <- 0 * stocks[["Tidal Basin"]]
      polder_subsidence <- stocks[["Polder"]] * polder_subsidence_rate
    }
    
    # ETD shoaling event every 30 years
    etd_shoaling <- if (t %% 30 == 0) etd_shoaling_rate * stocks[["ETD"]] else 0
    
    # Embankment of Saltmarshes
    if (stocks[["Saltmarsh"]] > embankment_threshold) {
      saltmarsh_to_polder <- 2/3 * embankment_threshold
      embankment <- 1
    } else {
      saltmarsh_to_polder <- 0
      embankment <- 0
    }
    
    # De-embankment of Polder
    subsidence_threshold <- stocks[["Polder"]] / polder_original
    if (subsidence_threshold < depoldering_value) {
      polder_to_saltmarsh <- 2/3 * stocks[["Polder"]]
      de_embankment <- 1
    } else {
      polder_to_saltmarsh <- 0
      de_embankment <- 0
    }
    
    # Update stocks
    stocks[["Offshore"]] <- stocks[["Offshore"]] + tail_to_offshore - offshore_to_etd - offshore_to_beach
    stocks[["ETD"]] <- stocks[["ETD"]] + offshore_to_etd - etd_to_island_head - etd_shoaling
    stocks[["Island Head"]] <- stocks[["Island Head"]] + etd_to_island_head + etd_shoaling - 
      island_head_to_beach - islands_head_to_dune
    stocks[["Beach Shore"]] <- stocks[["Beach Shore"]] + island_head_to_beach - beach_to_dune + 
      offshore_to_beach - beach_to_offshore - beach_to_tail + dune_to_beach
    stocks[["Dune Complex"]] <- stocks[["Dune Complex"]] + beach_to_dune + islands_head_to_dune - 
      dune_to_saltmarsh - dune_to_beach - dune_to_polder
    stocks[["Tail"]] <- stocks[["Tail"]] + beach_to_tail - tail_to_offshore - tail_to_saltmarsh
    stocks[["Saltmarsh"]] <- stocks[["Saltmarsh"]] + tidal_to_saltmarsh + polder_to_saltmarsh - 
      saltmarsh_to_polder + dune_to_saltmarsh + tail_to_saltmarsh
    stocks[["Tidal Basin"]] <- stocks[["Tidal Basin"]] - tidal_to_saltmarsh - tidal_to_polder
    stocks[["Polder"]] <- stocks[["Polder"]] + tidal_to_polder - polder_to_saltmarsh - 
      polder_subsidence + saltmarsh_to_polder + dune_to_polder
    
    # Apply SLR
    if (dSLR != 0) {
      for (key in names(stocks)) {
        if (!(key %in% c("Tidal Basin", "Offshore", "ETD"))) {
          stock_elevation[[key]] <- (stocks[[key]] / stock_areas[[key]]) * 100000
          stock_elevation[[key]] <- stock_elevation[[key]] - dSLR
          stocks[[key]] <- stock_elevation[[key]] / 100000 * stock_areas[[key]]
        }
      }
    }
    
    # Record history
    for (key in names(stocks)) {
      history[[key]] <- c(history[[key]], stocks[[key]])
      history_elevation[[key]] <- c(history_elevation[[key]], stock_elevation[[key]])
    }
    
    # Update polder stock to new original value
    if (de_embankment == 1 || embankment == 1) {
      polder_original <- stocks[["Polder"]]
    }
  }
}

# Print results
cat("Beach Unstable:", beach_unstable_count, "year(s)\n")

# Create plotting data frame
plot_data <- data.frame()
for (key in names(history)) {
  if (!(key %in% c("Tidal Basin", "Offshore"))) {
    df <- data.frame(
      Years = 1:length(history[[key]]),
      Stock = history[[key]],
      Category = key
    )
    plot_data <- rbind(plot_data, df)
  }
}

# Create the plot
ggplot(plot_data, aes(x = Years, y = Stock, color = Category)) +
  geom_line(size = 1) +
  labs(
    title = "Equilibrium Sediment Stock Model",
    x = "Years",
    y = "Sediment Stock above NAP * 10^5 (cubic meters)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_brewer(palette = "Set1")

# Elevation plot
plot_data_elev <- data.frame()
for (key in names(history)) {
  if (!(key %in% c("Tidal Basin", "Offshore", "ETD"))) {
    pltarea <- rep(stock_areas[[key]], length(history[[key]]))
    pltvl <- history[[key]] * rep(100000, length(history[[key]])) / pltarea
    df <- data.frame(
      Years = 1:length(pltvl),
      Elevation = pltvl,
      Category = key
    )
    plot_data_elev <- rbind(plot_data_elev, df)
  }
}

# Create elevation plot
ggplot(plot_data_elev, aes(x = Years, y = Elevation, color = Category)) +
  geom_line(size = 1) +
  labs(
    title = "Equilibrium Sediment Stock Model",
    x = "Years",
    y = "Mean Elevation NAP (in meters)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  scale_color_brewer(palette = "Set1")