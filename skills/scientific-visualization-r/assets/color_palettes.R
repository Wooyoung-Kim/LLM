# Color Palettes for Publication-Quality Figures in R

# Okabe-Ito palette - colorblind-safe
# Optimized for all types of color vision deficiency
okabe_ito_palette <- function() {
  c(
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermilion
    "#CC79A7",  # Reddish Purple
    "#000000"   # Black
  )
}

# Wong palette - another colorblind-safe option
wong_palette <- function() {
  c(
    "#000000",  # Black
    "#E69F00",  # Orange
    "#56B4E9",  # Sky Blue
    "#009E73",  # Bluish Green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermilion
    "#CC79A7"   # Reddish Purple
  )
}

# Paul Tol's bright palette - colorblind-safe
tol_bright_palette <- function() {
  c(
    "#4477AA",  # Blue
    "#EE6677",  # Red
    "#228833",  # Green
    "#CCBB44",  # Yellow
    "#66CCEE",  # Cyan
    "#AA3377",  # Purple
    "#BBBBBB"   # Grey
  )
}

# Paul Tol's muted palette - colorblind-safe
tol_muted_palette <- function() {
  c(
    "#332288",  # Indigo
    "#88CCEE",  # Cyan
    "#44AA99",  # Teal
    "#117733",  # Green
    "#999933",  # Olive
    "#DDCC77",  # Sand
    "#CC6677",  # Rose
    "#882255",  # Wine
    "#AA4499"   # Purple
  )
}

# Viridis-inspired discrete palette
viridis_discrete <- function(n = 8) {
  viridisLite::viridis(n)
}

# Set color palette globally
set_palette <- function(palette_name = "okabe_ito") {
  palette <- switch(palette_name,
    "okabe_ito" = okabe_ito_palette(),
    "wong" = wong_palette(),
    "tol_bright" = tol_bright_palette(),
    "tol_muted" = tol_muted_palette(),
    okabe_ito_palette()  # default
  )
  
  options(ggplot2.discrete.colour = palette)
  options(ggplot2.discrete.fill = palette)
  
  invisible(palette)
}

# Nature journal color palette (from their guidelines)
nature_palette <- function() {
  c(
    "#0173B2",  # Blue
    "#DE8F05",  # Orange
    "#029E73",  # Green
    "#CC78BC",  # Purple
    "#CA9161",  # Brown
    "#949494",  # Grey
    "#ECE133",  # Yellow
    "#56B4E9"   # Light Blue
  )
}

# Science journal color palette
science_palette <- function() {
  c(
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#e377c2",  # Pink
    "#7f7f7f"   # Grey
  )
}
