# Scientific Visualization in R

An R-based skill for creating publication-quality scientific figures using ggplot2, cowplot, and patchwork.

## Overview

This skill provides tools, themes, and guidelines for creating journal-ready figures in R that meet publication standards for:
- Nature, Science, Cell, PLOS, and other major journals
- Colorblind accessibility
- Proper resolution and file formats
- Consistent typography and styling

## Quick Start

```r
# Load the publication theme
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Create a plot
library(ggplot2)
p <- ggplot(mtcars, aes(x = wt, y = mpg, color = factor(cyl))) +
  geom_point(size = 2) +
  labs(x = "Weight (1000 lbs)", y = "Fuel Efficiency (mpg)") +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette())

# Save for publication
source('scripts/save_publication_figure.R')
save_publication_figure(p, 'figure1', width = 3.5, height = 2.5)
```

## Directory Structure

```
scientific-visualization-r/
├── SKILL.md                          # Main skill documentation
├── README.md                         # This file
├── scripts/                          # Helper R scripts
│   ├── theme_publication.R          # ggplot2 themes for journals
│   └── save_publication_figure.R    # Export functions
├── assets/                           # Reusable resources
│   └── color_palettes.R             # Colorblind-safe palettes
└── references/                       # Detailed guides
    ├── publication_guidelines.md    # Best practices
    ├── color_palettes.md           # Color usage guide
    ├── journal_requirements.md     # Journal specs
    └── ggplot2_examples.md         # Complete examples
```

## Key Features

### Publication-Ready Themes

Journal-specific themes with appropriate font sizes, spacing, and styling:

```r
source('scripts/theme_publication.R')

# General publication theme
theme_publication()

# Journal-specific themes
theme_nature()    # Nature journal specifications
theme_science()   # Science journal specifications
theme_cell()      # Cell Press specifications
theme_plos()      # PLOS journals
```

### Colorblind-Safe Palettes

All recommended palettes are safe for colorblind readers:

```r
source('assets/color_palettes.R')

# Okabe-Ito palette (best overall)
scale_color_manual(values = okabe_ito_palette())

# Other options
scale_color_manual(values = tol_bright_palette())
scale_color_manual(values = nature_palette())
```

### Easy Export Functions

Automatically save figures in multiple formats with correct settings:

```r
source('scripts/save_publication_figure.R')

# Save in multiple formats
save_publication_figure(plot, 'figure1', 
                        width = 7, height = 5,
                        formats = c('pdf', 'png', 'tiff'),
                        dpi = 300)

# Or use journal-specific settings
save_for_journal(plot, 'figure1', 
                 journal = 'nature',
                 figure_type = 'single')
```

## Common Use Cases

### Multi-Panel Figures with patchwork

```r
library(patchwork)

# Combine plots
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))
```

### Statistical Comparisons

```r
library(ggpubr)

ggplot(df, aes(x = group, y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("A", "B"))) +
  theme_publication()
```

### Heatmaps

```r
library(ComplexHeatmap)
library(viridis)

Heatmap(matrix, col = viridis(100), 
        name = "Expression")
```

## Requirements

### Required R Packages

```r
# Core visualization
install.packages("ggplot2")
install.packages("cowplot")
install.packages("patchwork")

# Colors
install.packages("viridis")
install.packages("RColorBrewer")
install.packages("colorblindr")  # devtools::install_github("clauswilke/colorblindr")

# Statistical plots
install.packages("ggpubr")
install.packages("ggbeeswarm")

# Specialized
install.packages("ComplexHeatmap")  # BiocManager::install("ComplexHeatmap")
```

## Documentation

- **[SKILL.md](SKILL.md)** - Complete skill documentation with examples
- **[Publication Guidelines](references/publication_guidelines.md)** - Best practices for scientific figures
- **[Color Palettes Guide](references/color_palettes.md)** - Colorblind-safe color selection
- **[Journal Requirements](references/journal_requirements.md)** - Specifications by journal
- **[ggplot2 Examples](references/ggplot2_examples.md)** - 10+ complete working examples

## Related Skills

This R version complements the Python-based `scientific-visualization` skill. Use whichever matches your workflow:

- **Python version**: matplotlib, seaborn, plotly
- **R version** (this): ggplot2, cowplot, patchwork

## License

This skill is designed for scientific research and publication. All code and resources are freely usable for academic and commercial purposes.

## Contributing

Contributions welcome! Areas for expansion:
- Additional journal-specific themes
- More example plots
- Domain-specific palettes (genomics, neuroscience, etc.)
- Accessibility enhancements

## References

- Okabe, M., & Ito, K. (2008). Color Universal Design. J*FLY Data Depository for Drosophila researchers.
- Tol, P. (2021). Colour Schemes. Technical Note SRON/EPS/TN/09-002.
- Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer.
