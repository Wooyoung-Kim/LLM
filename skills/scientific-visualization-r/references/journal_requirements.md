# Journal Requirements for Scientific Figures

## Quick Reference Table

| Journal | Single Column | Double Column | Font Size | DPI (Raster) | Formats |
|---------|--------------|---------------|-----------|--------------|---------|
| Nature | 89 mm (3.5") | 183 mm (7.2") | 5-7 pt | 300-600 | PDF, TIFF |
| Science | 55 mm (2.2") | 175 mm (6.9") | 6-8 pt | 300 | PDF, EPS |
| Cell | 85 mm (3.3") | 178 mm (7.0") | 6-8 pt | 300-600 | PDF, TIFF |
| PLOS | 83 mm (3.3") | 173 mm (6.8") | 6-8 pt | 300-600 | TIFF, EPS, PDF |

## Nature Family

### Dimensions
- **Single column**: 89 mm (3.5 inches)
- **Double column**: 183 mm (7.2 inches)
- **Maximum height**: 247 mm (9.7 inches)

### File Format
- **Preferred**: PDF (vector) or TIFF (raster)
- **DPI**: 300-600 for combination figures, 600-1200 for line art
- **Color mode**: RGB

### Typography
- **Font**: Arial, Helvetica, or similar sans-serif
- **Minimum size**: 5-7 pt at final size
- **Panel labels**: Lowercase bold (a, b, c)

### R Code
```r
source('scripts/save_publication_figure.R')
save_for_journal(plot, "figure1", journal = "nature", 
                 figure_type = "single")
```

## Science

### Dimensions
- **Single column**: 55 mm (2.2 inches)
- **Double column**: 175 mm (6.9 inches)
- **Maximum height**: 225 mm (8.9 inches)

### File Format
- **Preferred**: PDF, EPS
- **Acceptable**: TIFF (300 DPI minimum)
- **Color mode**: RGB or CMYK

### Typography
- **Font**: Arial or Helvetica
- **Minimum size**: 6-8 pt at final size
- **Panel labels**: Uppercase bold (A, B, C)

### R Code
```r
save_for_journal(plot, "figure1", journal = "science",
                 figure_type = "single")
```

## Cell Press

### Dimensions
- **Single column**: 85 mm (3.3 inches)
- **Double column**: 178 mm (7.0 inches)
- **Maximum height**: 234 mm (9.2 inches)

### File Format
- **Preferred**: PDF (vector)
- **Acceptable**: TIFF (300-600 DPI)
- **Color mode**: RGB

### Typography
- **Font**: Arial, 6-8 pt minimum
- **Panel labels**: Uppercase bold (A, B, C)

### R Code
```r
save_for_journal(plot, "figure1", journal = "cell",
                 figure_type = "double")
```

## PLOS (Public Library of Science)

### Dimensions
- **Single column**: 83 mm (3.27 inches)
- **Double column**: 173 mm (6.83 inches)
- **Maximum height**: 223 mm (8.78 inches)

### File Format
- **Preferred**: TIFF (uncompressed or LZW compression)
- **Acceptable**: EPS, PDF
- **DPI**: 300-600 for combination, 600-1200 for line art

### Typography
- **Font**: Arial, Helvetica, or similar
- **Minimum size**: 6-8 pt
- **Panel labels**: Uppercase bold (A, B, C)

### R Code
```r
save_for_journal(plot, "figure1", journal = "plos",
                 figure_type = "double")
```

## General Best Practices

### Resolution Guidelines

| Figure Type | Minimum DPI | Recommended Format |
|-------------|-------------|-------------------|
| Line art (graphs, plots) | 600-1200 | PDF (vector) |
| Halftone (photos) | 300 | TIFF |
| Combination (graphs + photos) | 300-600 | PDF or TIFF |

### Color Mode

- **RGB**: Most common for online/digital publication
- **CMYK**: Sometimes required for print journals (check guidelines)

In R, save as RGB by default:
```r
ggsave("figure.pdf", device = cairo_pdf)  # Ensures RGB
```

### Font Embedding

Always embed fonts in PDFs:
```r
# Using cairo for better font handling
ggsave("figure.pdf", device = cairo_pdf, 
       width = 7, height = 5)
```

## Specific Journal Notes

### Nature

- Very strict about dimensions
- Prefers lowercase panel labels (unique among major journals)
- RGB color mode required
- High-quality TIFF acceptable

### Science

- Smallest single-column width (55 mm)
- Plan for very compact figures
- Text must be legible at small size
- EPS or PDF preferred

### Cell

- Moderate sizing requirements
- Very flexible on file formats
- RGB for online, CMYK acceptable for print

### PLOS

- Open access, more flexible
- Strongly prefers TIFF
- Must provide figures in original file format
- All figures published online in color at no charge

## Checking Compliance

```r
source('scripts/save_publication_figure.R')

# Check if your figure meets requirements
check_figure_size(plot, journal = "nature", figure_type = "single")
```

## Common Pitfalls

1. **Wrong dimensions**: Always check journal guidelines
2. **Font too small**: Test at final print size (zoom out to 100%)
3. **Low resolution**: Vector for plots, ≥300 DPI for images
4. **JPEG format**: Never use JPEG for data visualization
5. **RGB vs CMYK**: Most journals want RGB for online publication

## Tips for Multi-Journal Submission

If you're unsure where to submit:

1. **Design for the smallest**: Use Science's dimensions (55 mm)
2. **Use vector formats**: Easily scalable
3. **Keep fonts ≥7 pt**: Readable in all journals
4. **Use colorblind-safe palettes**: Required by most journals
5. **Save in multiple formats**: PDF + high-res PNG/TIFF

## Using save_for_journal() Function

The helper function automatically applies correct settings:

```r
source('scripts/save_publication_figure.R')

# Automatically sets dimensions, DPI, and formats
save_for_journal(
  plot = my_plot,
  filename = "figure1",
  journal = "nature",      # or "science", "cell", "plos"
  figure_type = "single",  # or "double", "combination"
  path = "./figures"
)
```

This will create files with appropriate:
- Width and height
- DPI
- File formats
- Naming convention
