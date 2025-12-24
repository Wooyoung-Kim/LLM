# Color Palettes for Scientific Visualization in R

## Why Colorblind-Friendly Palettes Matter

Approximately **8% of males** and **0.5% of females** have some form of color vision deficiency (CVD). The most common is red-green color blindness (deuteranopia/protanopia).

**Key principle:** Never rely solely on color to convey information. Use additional visual cues like:
- Line types (solid, dashed, dotted)
- Point shapes (circle, square, triangle)
- Patterns in bar fills
- Direct labels

## Recommended Discrete Palettes

### Okabe-Ito Palette (Best Overall)

Designed specifically for all types of color blindness. This is the **gold standard** for scientific figures.

```r
source('assets/color_palettes.R')
colors <- okabe_ito_palette()

# In ggplot2
ggplot(data, aes(x, y, color = group)) +
  geom_point() +
  scale_color_manual(values = okabe_ito_palette())
```

**Colors (hex):**
- `#E69F00` - Orange
- `#56B4E9` - Sky Blue
- `#009E73` - Bluish Green
- `#F0E442` - Yellow
- `#0072B2` - Blue
- `#D55E00` - Vermilion
- `#CC79A7` - Reddish Purple
- `#000000` - Black

**Also available in colorblindr package:**
```r
library(colorblindr)
scale_color_OkabeIto()
scale_fill_OkabeIto()
```

### Paul Tol's Palettes

Scientifically designed palettes with different use cases.

**Bright palette (up to 7 colors):**
```r
source('assets/color_palettes.R')
scale_color_manual(values = tol_bright_palette())
```

**Muted palette (up to 9 colors):**
```r
scale_color_manual(values = tol_muted_palette())
```

### RColorBrewer Colorblind-Safe Sets

```r
library(RColorBrewer)

# For discrete data
scale_color_brewer(palette = "Set2")  # 8 colors
scale_fill_brewer(palette = "Dark2")  # 8 colors

# View all palettes
display.brewer.all(colorblindFriendly = TRUE)
```

## Continuous Color Scales

### Perceptually Uniform (Recommended)

These colormaps have **uniform perceptual brightness**, meaning differences in data values appear as equal visual differences.

```r
library(viridis)

# For sequential data
scale_color_viridis_c(option = "viridis")  # Default
scale_color_viridis_c(option = "magma")    # Good for dark backgrounds
scale_color_viridis_c(option = "plasma")   # Warmer tones
scale_color_viridis_c(option = "cividis")  # Optimized for CVD

# For heatmaps
scale_fill_viridis_c()
```

**When to use each:**
- **viridis** (blue to yellow): General purpose, good default
- **plasma** (purple to yellow): Warmer alternative
- **magma** (black to white): For dark backgrounds
- **cividis** (blue to yellow): Best for all CVD types
- **inferno** (black to yellow): High contrast

### Diverging Scales

For data with a meaningful center (e.g., correlation, fold change).

**Colorblind-safe diverging palettes:**
```r
library(RColorBrewer)

# Good choices
scale_fill_distiller(palette = "RdBu", direction = -1)   # Red-Blue
scale_fill_distiller(palette = "PuOr", direction = -1)   # Purple-Orange
scale_fill_distiller(palette = "BrBG", direction = -1)   # Brown-Blue-Green
scale_fill_distiller(palette = "PRGn", direction = -1)   # Purple-Green

# Always set center for diverging scales
scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", 
                     midpoint = 0)
```

**Avoid:**
- Red-Green (RdYlGn) - problematic for red-green CVD
- Spectral - contains red-green transitions

## Palettes to Avoid

### Never Use

1. **Rainbow/Jet colormap** - Not perceptually uniform, problematic for CVD
   ```r
   # DON'T DO THIS
   scale_color_gradientn(colors = rainbow(7))
   ```

2. **Red-Green combinations** - ~8% of males cannot distinguish
   ```r
   # DON'T DO THIS
   scale_color_manual(values = c("red", "green"))
   ```

3. **Pure red + pure blue** - Difficult to distinguish for some CVD types

## Testing Your Figures

### 1. Grayscale Test

```r
# View your plot in grayscale
library(colorspace)
desaturate(your_plot)

# Or manually
print(your_plot)
# Then: File > Print Preview (set to black & white)
```

Your figure should still be interpretable in grayscale.

### 2. Colorblind Simulation

```r
library(colorblindr)

# Simulate different types of color blindness
cvd_grid(your_plot)

# Individual simulations
deutan(your_plot)  # Red-green (most common)
protan(your_plot)  # Red-green (another type)
tritan(your_plot)  # Blue-yellow (rare)
```

## Guidelines by Plot Type

### Line Plots

```r
# Use both color AND line type
ggplot(data, aes(x, y, color = group, linetype = group)) +
  geom_line() +
  scale_color_manual(values = okabe_ito_palette()) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))
```

### Scatter Plots

```r
# Use both color AND shape
ggplot(data, aes(x, y, color = group, shape = group)) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito_palette()) +
  scale_shape_manual(values = c(16, 17, 15, 18))
```

### Bar Plots

```r
# Use patterns for additional encoding
library(ggpattern)

ggplot(data, aes(x = category, y = value, fill = group)) +
  geom_bar_pattern(
    aes(pattern = group),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_manual(values = okabe_ito_palette()) +
  scale_pattern_manual(values = c("none", "stripe", "crosshatch"))
```

### Heatmaps

```r
# Use perceptually uniform colormaps
library(pheatmap)
library(viridis)

pheatmap(
  data_matrix,
  color = viridis(100),
  # OR for diverging data
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-3, 3, length.out = 101)
)
```

## Journal-Specific Recommendations

### Nature

Nature provides specific color recommendations:
```r
source('assets/color_palettes.R')
scale_color_manual(values = nature_palette())
```

### General Journals

Most journals accept any colorblind-friendly palette. Okabe-Ito is always safe.

## Quick Reference Table

| Use Case | Recommended Palette | R Code |
|----------|---------------------|--------|
| Up to 8 categories | Okabe-Ito | `okabe_ito_palette()` |
| 9+ categories | Split into facets or use continuous | - |
| Sequential continuous | Viridis | `scale_color_viridis_c()` |
| Diverging continuous | RdBu or PuOr | `scale_fill_distiller(palette = "RdBu")` |
| Correlation matrix | RdBu centered at 0 | `scale_fill_gradient2(...)` |
| Heatmap | Viridis or Magma | `scale_fill_viridis_c()` |

## Additional Resources

- **colorblindr**: `install.packages("colorblindr")`
- **viridis**: `install.packages("viridis")`
- **RColorBrewer**: `install.packages("RColorBrewer")`
- **colorspace**: `install.packages("colorspace")`

Online tools:
- [Coolors.co](https://coolors.co/) - Generate color palettes
- [ColorBrewer](https://colorbrewer2.org/) - Interactive palette selection
- [Viz Palette](https://projects.susielu.com/viz-palette) - Test palettes for CVD
