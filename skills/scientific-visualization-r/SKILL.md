---
name: scientific-visualization-r
description: "Create publication figures with ggplot2/cowplot/patchwork in R. Multi-panel layouts, error bars, significance markers, colorblind-safe palettes, export PDF/EPS/TIFF for journal-ready scientific plots."
---

# Scientific Visualization (R)

## Overview

Scientific visualization in R transforms data into clear, accurate figures for publication. Create journal-ready plots with multi-panel layouts, error bars, significance markers, and colorblind-safe palettes using ggplot2, cowplot, and patchwork. Export as PDF/EPS/TIFF for manuscripts.

## When to Use This Skill

This skill should be used when:
- Creating plots or visualizations for scientific manuscripts in R
- Preparing figures for journal submission (Nature, Science, Cell, PLOS, etc.)
- Ensuring figures are colorblind-friendly and accessible
- Making multi-panel figures with consistent styling
- Exporting figures at correct resolution and format
- Following specific publication guidelines
- Improving existing R figures to meet publication standards
- Creating figures that need to work in both color and grayscale

## Quick Start Guide

### Basic Publication-Quality Figure

```r
library(ggplot2)
library(cowplot)

# Apply publication theme
source('scripts/theme_publication.R')

# Create figure with appropriate size (single column = 3.5 inches)
p <- ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Weight (1000 lbs)",
    y = "Fuel Efficiency (mpg)",
    color = "Cylinders"
  ) +
  theme_publication() +
  scale_color_colorblind()

# Save in publication formats
source('scripts/save_publication_figure.R')
save_publication_figure(p, 'figure1', width = 3.5, height = 2.5, 
                        formats = c('pdf', 'png'), dpi = 300)
```

### Using Pre-configured Themes

Apply journal-specific themes using the provided helper functions:

```r
library(ggplot2)
source('scripts/theme_publication.R')

# Option 1: Use theme function directly
p <- ggplot(data, aes(x, y)) +
  geom_point() +
  theme_publication()

# Option 2: Configure for specific journal
theme_nature <- function() {
  theme_publication(base_size = 7, base_family = "Arial")
}

p <- ggplot(data, aes(x, y)) +
  geom_point() +
  theme_nature()
```

### Quick Start with Statistical Plots

For statistical comparisons, use ggplot2 with publication styling:

```r
library(ggplot2)
library(ggbeeswarm)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Create statistical comparison figure
p <- ggplot(df, aes(x = treatment, y = response, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_beeswarm(alpha = 0.4, size = 1) +
  labs(
    x = "Treatment",
    y = expression("Response ("*mu*"M)")
  ) +
  theme_publication() +
  scale_fill_manual(values = okabe_ito_palette()) +
  theme(legend.position = "none")

# Save figure
save_publication_figure(p, 'treatment_comparison', 
                        width = 3.5, height = 3, 
                        formats = c('pdf', 'png'), dpi = 300)
```

## Core Principles and Best Practices

### 1. Resolution and File Format

**Critical requirements** (detailed in `references/publication_guidelines.md`):
- **Raster images** (photos, microscopy): 300-600 DPI
- **Line art** (graphs, plots): 600-1200 DPI or vector format
- **Vector formats** (preferred): PDF, EPS, SVG
- **Raster formats**: TIFF, PNG (never JPEG for scientific data)

**Implementation:**
```r
# Use the save_publication_figure.R script for correct settings
source('scripts/save_publication_figure.R')

# Saves in multiple formats with proper DPI
save_publication_figure(p, 'myfigure', 
                        width = 7, height = 5,
                        formats = c('pdf', 'png'), dpi = 300)

# Or save for specific journal requirements
save_for_journal(p, 'figure1', journal = 'nature', 
                 figure_type = 'combination')
```

### 2. Color Selection - Colorblind Accessibility

**Always use colorblind-friendly palettes** (detailed in `references/color_palettes.md`):

**Recommended: Okabe-Ito palette** (distinguishable by all types of color blindness):
```r
# Option 1: Use assets/color_palettes.R
source('assets/color_palettes.R')

# For discrete colors
p + scale_color_manual(values = okabe_ito_palette())
p + scale_fill_manual(values = okabe_ito_palette())

# Option 2: Use colorblindr package
library(colorblindr)
p + scale_color_OkabeIto()
p + scale_fill_OkabeIto()
```

**For heatmaps/continuous data:**
- Use perceptually uniform colormaps: `viridis`, `magma`, `plasma`, `cividis`
- Avoid red-green diverging maps (use `PuOr`, `RdBu`, `BrBG` instead)
- Never use the `rainbow` palette

```r
library(viridis)

# For continuous data
p + scale_color_viridis_c()
p + scale_fill_viridis_c()

# For diverging data
library(RColorBrewer)
p + scale_fill_distiller(palette = "RdBu", direction = -1)
```

**Always test figures in grayscale** to ensure interpretability.

### 3. Typography and Text

**Font guidelines** (detailed in `references/publication_guidelines.md`):
- Sans-serif fonts: Arial, Helvetica, Calibri
- Minimum sizes at **final print size**:
  - Axis labels: 7-9 pt
  - Tick labels: 6-8 pt
  - Panel labels: 8-12 pt (bold)
- Sentence case for labels: "Time (hours)" not "TIME (HOURS)"
- Always include units in parentheses

**Implementation:**
```r
# Set fonts in theme
theme_publication <- function(base_size = 8, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) +
  theme(
    # Text elements
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    # ... other theme elements
  )
}
```

### 4. Figure Dimensions

**Journal-specific widths** (detailed in `references/journal_requirements.md`):
- **Nature**: Single 89 mm (3.5"), Double 183 mm (7.2")
- **Science**: Single 55 mm (2.2"), Double 175 mm (6.9")
- **Cell**: Single 85 mm (3.3"), Double 178 mm (7.0")

**Set figure dimensions:**
```r
# Create figure with specific dimensions
# ggsave respects inches by default
ggsave('figure.pdf', plot = p, 
       width = 3.5, height = 2.5, units = "in")

# For Nature single column
ggsave('nature_fig.pdf', plot = p,
       width = 89, height = 60, units = "mm")
```

### 5. Multi-Panel Figures

**Best practices:**
- Label panels with bold letters: **A**, **B**, **C** (uppercase for most journals, lowercase for Nature)
- Maintain consistent styling across all panels
- Align panels along edges where possible
- Use adequate white space between panels

**Example implementation using patchwork:**
```r
library(patchwork)

# Create individual panels
p1 <- ggplot(...) + theme_publication() + ggtitle("A")
p2 <- ggplot(...) + theme_publication() + ggtitle("B")
p3 <- ggplot(...) + theme_publication() + ggtitle("C")
p4 <- ggplot(...) + theme_publication() + ggtitle("D")

# Combine with patchwork
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 10))

# Save combined figure
ggsave('multi_panel.pdf', combined, width = 7, height = 6)
```

**Alternative using cowplot:**
```r
library(cowplot)

# Create individual panels
p1 <- ggplot(...) + theme_publication()
p2 <- ggplot(...) + theme_publication()
p3 <- ggplot(...) + theme_publication()
p4 <- ggplot(...) + theme_publication()

# Combine with cowplot
combined <- plot_grid(p1, p2, p3, p4, 
                      labels = c('A', 'B', 'C', 'D'),
                      label_size = 10,
                      label_fontface = 'bold',
                      ncol = 2)

# Save
save_plot('multi_panel.pdf', combined, base_width = 7, base_height = 6)
```

## Common Tasks

### Task 1: Create a Publication-Ready Line Plot

See `references/ggplot2_examples.md` Example 1 for complete code.

**Key steps:**
1. Apply publication theme
2. Set appropriate figure size for target journal
3. Use colorblind-friendly colors
4. Add error bars with correct representation (SEM, SD, or CI)
5. Label axes with units
6. Save in vector format

**Example with ribbons for confidence intervals:**
```r
library(ggplot2)
library(dplyr)

# Calculate summary statistics
summary_df <- df %>%
  group_by(time, treatment) %>%
  summarise(
    mean = mean(value),
    se = sd(value) / sqrt(n()),
    .groups = 'drop'
  )

p <- ggplot(summary_df, aes(x = time, y = mean, 
                            color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), 
              alpha = 0.2, color = NA) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  labs(
    x = "Time (hours)",
    y = "Measurement (AU)"
  ) +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette()) +
  scale_fill_manual(values = okabe_ito_palette())
```

### Task 2: Create a Multi-Panel Figure

See `references/ggplot2_examples.md` Example 2 for complete code.

**Key steps:**
1. Use `patchwork` or `cowplot` for flexible layout
2. Ensure consistent styling across panels
3. Add bold panel labels (A, B, C, etc.)
4. Align related panels
5. Verify all text is readable at final size

### Task 3: Create a Heatmap with Proper Colormap

See `references/ggplot2_examples.md` Example 4 for complete code.

**Key steps:**
1. Use perceptually uniform colormap (`viridis`, `plasma`, `cividis`)
2. Include labeled colorbar
3. For diverging data, use colorblind-safe diverging map
4. Set appropriate center value for diverging maps
5. Test appearance in grayscale

**Example correlation matrix:**
```r
library(ggplot2)
library(reshape2)

# Calculate correlation
corr_matrix <- cor(df)
corr_df <- melt(corr_matrix)

p <- ggplot(corr_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       limits = c(-1, 1), name = "Correlation") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  coord_fixed()
```

### Task 4: Prepare Figure for Specific Journal

**Workflow:**
1. Check journal requirements: `references/journal_requirements.md`
2. Configure ggplot2 for journal:
   ```r
   source('scripts/theme_publication.R')
   theme_set(theme_nature())
   ```
3. Create figure (will auto-size correctly)
4. Export with journal specifications:
   ```r
   save_for_journal(p, 'figure1', journal = 'nature', 
                    figure_type = 'line_art')
   ```

### Task 5: Add Statistical Significance Markers

**Using ggsignif package:**
```r
library(ggplot2)
library(ggsignif)

p <- ggplot(df, aes(x = treatment, y = response)) +
  geom_boxplot() +
  geom_signif(
    comparisons = list(c("Control", "Treatment1"), 
                       c("Control", "Treatment2")),
    map_signif_level = TRUE,
    test = "t.test"
  ) +
  theme_publication()
```

**Manual significance bars:**
```r
# Add custom significance brackets
p <- p +
  annotate("segment", x = 1, xend = 2, 
           y = max_y * 1.1, yend = max_y * 1.1) +
  annotate("text", x = 1.5, y = max_y * 1.12, 
           label = "***", size = 4)
```

### Task 6: Create Colorblind-Friendly Visualizations

**Strategy:**
1. Use approved palettes from `assets/color_palettes.R`
2. Add redundant encoding (line types, shapes, patterns)
3. Test with colorblind simulator
4. Ensure grayscale compatibility

**Example:**
```r
source('assets/color_palettes.R')

# Add redundant encoding beyond color
p <- ggplot(df, aes(x = time, y = value, 
                    color = treatment, 
                    linetype = treatment, 
                    shape = treatment)) +
  geom_line(size = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = okabe_ito_palette()) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_shape_manual(values = c(16, 17, 15)) +
  theme_publication()
```

## Statistical Rigor

**Always include:**
- Error bars (SD, SEM, or CI - specify which in caption)
- Sample size (n) in figure or caption
- Statistical significance markers (*, **, ***)
- Individual data points when possible (not just summary statistics)

**Example with statistics:**
```r
library(ggplot2)
library(ggbeeswarm)

p <- ggplot(df, aes(x = treatment, y = response)) +
  # Show individual points
  geom_beeswarm(alpha = 0.4, size = 1.5) +
  # Add mean and error bars
  stat_summary(fun = mean, geom = "point", size = 3, color = "red") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, color = "red") +
  theme_publication()
```

## Working with Different R Packages

### ggplot2
- Most flexible and consistent grammar of graphics
- Best for publication-quality figures
- Use provided themes for consistent formatting
- See `references/ggplot2_examples.md` for extensive examples

### cowplot
- Built on top of ggplot2 for publication-ready themes
- Excellent for multi-panel figures
- Provides `plot_grid()` for combining plots
- Use `save_plot()` for proper sizing

```r
library(cowplot)

# Apply cowplot theme (clean, publication-ready)
theme_set(theme_cowplot())

# Combine plots
combined <- plot_grid(p1, p2, p3, p4, 
                      labels = c('A', 'B', 'C', 'D'),
                      ncol = 2)
```

### patchwork
- Modern approach to combining ggplot2 plots
- Intuitive syntax with `+`, `|`, and `/` operators
- Automatic alignment and sizing

```r
library(patchwork)

# Simple combination
combined <- p1 + p2 + p3 + p4 +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A')
```

### ggpubr
- Publication-ready themes and utilities
- Statistical comparison functions built-in
- Easy significance markers

```r
library(ggpubr)

p <- ggboxplot(df, x = "treatment", y = "response",
               fill = "treatment", palette = "jco") +
  stat_compare_means(comparisons = my_comparisons) +
  theme_publication()
```

### ComplexHeatmap
- Advanced heatmaps with clustering
- Publication-quality defaults
- Excellent for genomics data

```r
library(ComplexHeatmap)

Heatmap(matrix_data,
        name = "Expression",
        col = viridis(100),
        cluster_rows = TRUE,
        cluster_columns = TRUE)
```

## Resources

### References Directory

**Load these as needed for detailed information:**

- **`publication_guidelines.md`**: Comprehensive best practices
  - Resolution and file format requirements
  - Typography guidelines
  - Layout and composition rules
  - Statistical rigor requirements
  - Complete publication checklist

- **`color_palettes.md`**: Color usage guide
  - Colorblind-friendly palette specifications with hex values
  - Sequential and diverging colormap recommendations
  - Testing procedures for accessibility
  - Domain-specific palettes (genomics, microscopy)

- **`journal_requirements.md`**: Journal-specific specifications
  - Technical requirements by publisher
  - File format and DPI specifications
  - Figure dimension requirements
  - Quick reference table

- **`ggplot2_examples.md`**: Practical code examples
  - 10 complete working examples
  - Line plots, bar plots, heatmaps, multi-panel figures
  - Journal-specific figure examples
  - Tips for each package (ggplot2, cowplot, patchwork)

### Scripts Directory

**Use these helper scripts for automation:**

- **`save_publication_figure.R`**: Export utilities
  - `save_publication_figure()`: Save in multiple formats with correct DPI
  - `save_for_journal()`: Use journal-specific requirements automatically
  - `check_figure_size()`: Verify dimensions meet journal specs

- **`theme_publication.R`**: Pre-configured themes
  - `theme_publication()`: General publication theme
  - `theme_nature()`, `theme_science()`, `theme_cell()`: Journal-specific themes
  - Run `source('scripts/theme_publication.R')` in your scripts

### Assets Directory

**Use these files in figures:**

- **`color_palettes.R`**: Importable color definitions
  - All recommended palettes as R vectors
  - `okabe_ito_palette()` helper function
  - Can be sourced directly into scripts

## Workflow Summary

**Recommended workflow for creating publication figures:**

1. **Plan**: Determine target journal, figure type, and content
2. **Configure**: Apply appropriate theme for journal
   ```r
   source('scripts/theme_publication.R')
   theme_set(theme_nature())
   ```
3. **Create**: Build figure with proper labels, colors, statistics
4. **Verify**: Check size, fonts, colors, accessibility
5. **Export**: Save in required formats
   ```r
   source('scripts/save_publication_figure.R')
   save_for_journal(p, 'figure1', 'nature', 'combination')
   ```
6. **Review**: View at final size in manuscript context

## Common Pitfalls to Avoid

1. **Font too small**: Text unreadable when printed at final size
2. **JPEG format**: Never use JPEG for graphs/plots (creates artifacts)
3. **Red-green colors**: ~8% of males cannot distinguish
4. **Low resolution**: Pixelated figures in publication
5. **Missing units**: Always label axes with units
6. **3D effects**: Distorts perception, avoid completely
7. **Chart junk**: Remove unnecessary gridlines, decorations
8. **Truncated axes**: Start bar charts at zero unless scientifically justified
9. **Inconsistent styling**: Different fonts/colors across figures in same manuscript
10. **No error bars**: Always show uncertainty

## Final Checklist

Before submitting figures, verify:

- [ ] Resolution meets journal requirements (300+ DPI)
- [ ] File format is correct (vector for plots, TIFF for images)
- [ ] Figure size matches journal specifications
- [ ] All text readable at final size (≥6 pt)
- [ ] Colors are colorblind-friendly
- [ ] Figure works in grayscale
- [ ] All axes labeled with units
- [ ] Error bars present with definition in caption
- [ ] Panel labels present and consistent
- [ ] No chart junk or 3D effects
- [ ] Fonts consistent across all figures
- [ ] Statistical significance clearly marked
- [ ] Legend is clear and complete

Use this skill to ensure scientific figures meet the highest publication standards while remaining accessible to all readers.
