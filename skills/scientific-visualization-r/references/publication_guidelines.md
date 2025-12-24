# Publication Guidelines for Scientific Figures in R

## Resolution and File Format Requirements

### Raster vs Vector

**Vector formats (preferred for plots):**
- PDF (Portable Document Format)
- EPS (Encapsulated PostScript)
- SVG (Scalable Vector Graphics)
- Advantages: Infinite scaling, small file size, editable

**Raster formats (for images/photos):**
- TIFF (Tagged Image File Format)
- PNG (Portable Network Graphics)
- Never use JPEG for scientific data (lossy compression)

### DPI Requirements

- **Combination figures** (graphs + photos): 300-600 DPI
- **Line art** (graphs, plots): 600-1200 DPI or vector
- **Photographs/microscopy**: 300-600 DPI
- **Presentations/posters**: 150-200 DPI acceptable

### Saving in R

```r
# Vector format (preferred for plots)
ggsave("figure.pdf", plot = p, width = 7, height = 5, units = "in")

# High-resolution raster
ggsave("figure.png", plot = p, width = 7, height = 5, 
       units = "in", dpi = 300)

# TIFF for journals requiring raster
ggsave("figure.tiff", plot = p, width = 7, height = 5,
       units = "in", dpi = 300, compression = "lzw")
```

## Typography Guidelines

### Font Selection

**Recommended fonts:**
- Arial (most common)
- Helvetica
- Calibri
- Liberation Sans (open-source Arial alternative)

**Avoid:**
- Serif fonts (Times, Georgia) for graphs
- Decorative fonts
- System-specific fonts

### Font Sizes

**Minimum at final print size:**
- Panel labels (A, B, C): 8-12 pt (bold)
- Axis titles: 7-9 pt
- Axis tick labels: 6-8 pt
- Legend text: 6-8 pt
- Annotations: 6-8 pt

### Best Practices

- Use sentence case: "Time (hours)" not "TIME (HOURS)"
- Always include units in parentheses
- Use proper Greek letters: expression("Expression ("*mu*"M)")
- Consistent font family across all figures in manuscript

## Color Guidelines

### Colorblind-Friendly Palettes

**Always use colorblind-safe colors.**

Approximately 8% of males and 0.5% of females have some form of color vision deficiency.

**Recommended palettes:**

1. **Okabe-Ito palette** (best)
   - Distinguishable by all types of color blindness
   - Available in `colorblindr` package

2. **Paul Tol's palettes**
   - Bright, muted, and light palettes
   - Scientifically designed for accessibility

3. **Viridis family**
   - For continuous data
   - Perceptually uniform
   - `viridis`, `plasma`, `magma`, `cividis`

**Colors to avoid:**
- Red-green combinations
- Red-blue combinations for diverging data
- Rainbow/jet colormaps

### Testing

Test all figures in:
1. Grayscale (print preview)
2. Colorblind simulators (colorblindr package)
3. At final print size

## Layout and Composition

### Figure Dimensions

**Journal-specific widths:**

| Journal | Single Column | Double Column |
|---------|--------------|---------------|
| Nature  | 89 mm (3.5") | 183 mm (7.2") |
| Science | 55 mm (2.2") | 175 mm (6.9") |
| Cell    | 85 mm (3.3") | 178 mm (7.0") |
| PLOS    | 83 mm (3.3") | 173 mm (6.8") |

### Multi-Panel Figures

**Panel labels:**
- Use bold uppercase letters: **A**, **B**, **C**
  - Exception: Nature uses lowercase: **a**, **b**, **c**
- Place consistently (usually top-left)
- Size: 10-12 pt

**Spacing:**
- Adequate white space between panels
- Align panels along common edges
- Consistent styling across all panels

### White Space

- Don't cram too much into one figure
- Leave margins around plot area
- Balance visual weight

## Statistical Rigor

### Error Bars

**Always include error bars and specify type in caption:**
- Standard Deviation (SD): Spread of individual measurements
- Standard Error of Mean (SEM): Uncertainty in mean
- Confidence Intervals (CI): Range likely containing true mean
  - 95% CI most common

**In ggplot2:**
```r
# SD
stat_summary(fun.data = mean_sdl, geom = "errorbar", width = 0.2)

# SEM
stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2)

# 95% CI
stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.2)
```

### Sample Size

- Include n in figure or caption
- Show individual data points when possible
- Don't hide data behind summary statistics

### Statistical Significance

**Notation:**
- * p < 0.05
- ** p < 0.01
- *** p < 0.001
- ns = not significant

**Best practice:** Report exact p-values in caption

## Common Requirements by Journal

### Nature

- Single column: 89 mm wide
- Double column: 183 mm wide
- Font: Arial, 5-7 pt at final size
- File formats: PDF (vector) or TIFF (300-600 DPI)
- Color: RGB
- Panel labels: lowercase bold

### Science

- Single column: 55 mm wide
- Double column: 175 mm wide
- Font: Arial, 6-8 pt
- File formats: PDF, EPS, or TIFF
- Minimal styling

### Cell

- Single column: 85 mm wide
- Double column: 178 mm wide
- Font: Arial, 6-8 pt
- File formats: PDF or TIFF (300-600 DPI)
- Color: RGB

### PLOS

- Single column: 83 mm wide
- Double column: 173 mm wide
- Font: Arial, 6-8 pt
- File formats: TIFF, EPS, PDF
- Color: RGB
- Very flexible on styling

## Complete Checklist

Before submission:

- [ ] Resolution: 300+ DPI for raster, vector for plots
- [ ] File format: Correct for journal
- [ ] Dimensions: Match journal specifications
- [ ] Fonts: ≥6 pt at final size, sans-serif
- [ ] Colors: Colorblind-friendly
- [ ] Grayscale: Figure interpretable without color
- [ ] Labels: All axes labeled with units
- [ ] Error bars: Present with type specified in caption
- [ ] Sample size: n specified
- [ ] Statistics: Significance clearly marked
- [ ] Panel labels: Present and consistent
- [ ] Legend: Clear and complete
- [ ] File names: Descriptive and versioned
- [ ] White space: Adequate margins and spacing
- [ ] Consistency: Styling matches other figures in manuscript
