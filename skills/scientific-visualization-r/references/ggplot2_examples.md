# ggplot2 Examples for Publication-Quality Figures

## Example 1: Basic Line Plot with Error Bars

```r
library(ggplot2)
library(dplyr)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Prepare data with summary statistics
summary_data <- df %>%
  group_by(time, treatment) %>%
  summarise(
    mean_response = mean(response),
    se = sd(response) / sqrt(n()),
    .groups = 'drop'
  )

# Create plot
p <- ggplot(summary_data, aes(x = time, y = mean_response, 
                               color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = mean_response - se, ymax = mean_response + se),
              alpha = 0.2, color = NA) +
  geom_line(size = 0.8) +
  geom_point(size = 2, shape = 21, fill = "white") +
  labs(
    x = "Time (hours)",
    y = "Response (AU)",
    color = "Treatment",
    fill = "Treatment"
  ) +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette()) +
  scale_fill_manual(values = okabe_ito_palette())

# Save
source('scripts/save_publication_figure.R')
save_publication_figure(p, 'figure1_lineplot', 
                        width = 3.5, height = 2.5,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 2: Bar Plot with Individual Points

```r
library(ggplot2)
library(ggbeeswarm)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

p <- ggplot(df, aes(x = treatment, y = response, fill = treatment)) +
  # Bar shows mean
  stat_summary(fun = mean, geom = "bar", alpha = 0.7, color = "black", size = 0.3) +
  # Error bars show SEM
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.5) +
  # Individual points
  geom_beeswarm(alpha = 0.4, size = 1.5, color = "black") +
  labs(
    x = "Treatment",
    y = expression("Response ("*mu*"M)")
  ) +
  theme_publication() +
  scale_fill_manual(values = okabe_ito_palette()) +
  theme(legend.position = "none")

save_publication_figure(p, 'figure2_barplot',
                        width = 3.5, height = 3,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 3: Multi-Panel Figure with patchwork

```r
library(ggplot2)
library(patchwork)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Panel A: Scatter plot
p1 <- ggplot(df, aes(x = dose, y = response, color = cell_line)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, size = 0.8) +
  labs(x = "Dose (nM)", y = "Response (%)", color = "Cell Line") +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette())

# Panel B: Box plot
p2 <- ggplot(df, aes(x = cell_line, y = response, fill = cell_line)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  labs(x = "Cell Line", y = "Response (%)") +
  theme_publication() +
  scale_fill_manual(values = okabe_ito_palette()) +
  theme(legend.position = "none")

# Panel C: Time series
p3 <- ggplot(timeseries, aes(x = time, y = signal, color = treatment)) +
  geom_line(size = 0.8) +
  geom_point(size = 1.5) +
  labs(x = "Time (min)", y = "Signal (AU)", color = "Treatment") +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette())

# Panel D: Histogram
p4 <- ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  labs(x = "Value", y = "Count", fill = "Group") +
  theme_publication() +
  scale_fill_manual(values = okabe_ito_palette())

# Combine panels
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold', size = 10))

# Save
save_publication_figure(combined, 'figure3_multipanel',
                        width = 7, height = 6,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 4: Heatmap with Hierarchical Clustering

```r
library(ComplexHeatmap)
library(circlize)
library(viridis)

# Prepare data matrix
mat <- as.matrix(expression_data)

# Create color function
col_fun <- colorRamp2(
  c(min(mat), 0, max(mat)),
  c("#2166AC", "white", "#B2182B")
)

# Create heatmap
ht <- Heatmap(
  mat,
  name = "Expression",
  col = col_fun,
  
  # Clustering
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  
  # Styling
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6),
  
  # Annotations
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  )
)

# Save as PDF
pdf("figure4_heatmap.pdf", width = 7, height = 5)
draw(ht)
dev.off()

# Save as PNG
png("figure4_heatmap.png", width = 7, height = 5, 
    units = "in", res = 300)
draw(ht)
dev.off()
```

## Example 5: Correlation Plot with ggplot2

```r
library(ggplot2)
library(reshape2)
library(RColorBrewer)
source('scripts/theme_publication.R')

# Calculate correlation
corr_matrix <- cor(df[, numeric_columns])

# Get lower triangle only
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  cormat
}
lower_tri <- get_lower_tri(corr_matrix)

# Melt for ggplot2
corr_df <- melt(lower_tri, na.rm = TRUE)

# Create heatmap
p <- ggplot(corr_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), 
            size = 2.5, color = "black") +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = -1,
    limits = c(-1, 1),
    name = "Pearson\nCorrelation"
  ) +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

save_publication_figure(p, 'figure5_correlation',
                        width = 5, height = 4,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 6: Violin Plot with Statistics

```r
library(ggplot2)
library(ggpubr)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Define comparisons
my_comparisons <- list(
  c("Control", "Treatment1"),
  c("Control", "Treatment2"),
  c("Treatment1", "Treatment2")
)

p <- ggplot(df, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test",
                     label = "p.signif") +
  labs(
    x = "Group",
    y = "Value (AU)"
  ) +
  theme_publication() +
  scale_fill_manual(values = okabe_ito_palette()) +
  theme(legend.position = "none")

save_publication_figure(p, 'figure6_violin',
                        width = 4, height = 4,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 7: Faceted Plot

```r
library(ggplot2)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

p <- ggplot(df, aes(x = time, y = response, color = treatment)) +
  geom_line(size = 0.8) +
  geom_point(size = 1.5) +
  facet_wrap(~ cell_line, ncol = 3) +
  labs(
    x = "Time (hours)",
    y = "Response (AU)",
    color = "Treatment"
  ) +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette()) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 8)
  )

save_publication_figure(p, 'figure7_faceted',
                        width = 7, height = 4,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 8: Dose-Response Curve

```r
library(ggplot2)
library(drc)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Fit dose-response model
model <- drm(response ~ dose, data = df, fct = LL.4())

# Create prediction data
pred_data <- data.frame(
  dose = seq(min(df$dose), max(df$dose), length.out = 100)
)
pred_data$response <- predict(model, newdata = pred_data)

# Plot
p <- ggplot() +
  # Raw data points
  geom_point(data = df, aes(x = dose, y = response),
             alpha = 0.5, size = 2) +
  # Fitted curve
  geom_line(data = pred_data, aes(x = dose, y = response),
            color = okabe_ito_palette()[1], size = 1) +
  scale_x_log10() +
  labs(
    x = "Dose (nM)",
    y = "Response (%)"
  ) +
  theme_publication() +
  annotation_logticks(sides = "b")

save_publication_figure(p, 'figure8_dose_response',
                        width = 3.5, height = 2.5,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 9: PCA Plot

```r
library(ggplot2)
library(ggrepel)
source('scripts/theme_publication.R')
source('assets/color_palettes.R')

# Perform PCA
pca_result <- prcomp(df[, numeric_columns], scale. = TRUE)
pca_data <- data.frame(pca_result$x, group = df$group)

# Calculate variance explained
var_exp <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = rownames(pca_data))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(size = 2.5, max.overlaps = 10) +
  labs(
    x = sprintf("PC1 (%s%%)", var_exp[1]),
    y = sprintf("PC2 (%s%%)", var_exp[2]),
    color = "Group"
  ) +
  theme_publication() +
  scale_color_manual(values = okabe_ito_palette())

save_publication_figure(p, 'figure9_pca',
                        width = 4, height = 3.5,
                        formats = c('pdf', 'png'), dpi = 300)
```

## Example 10: Survival Curve

```r
library(survival)
library(survminer)
source('assets/color_palettes.R')

# Fit survival model
fit <- survfit(Surv(time, status) ~ group, data = survival_df)

# Create plot with survminer
p <- ggsurvplot(
  fit,
  data = survival_df,
  
  # Styling
  palette = okabe_ito_palette(),
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  
  # Labels
  xlab = "Time (months)",
  ylab = "Survival Probability",
  legend.title = "Group",
  legend.labs = c("Control", "Treatment"),
  
  # Font sizes
  font.main = 10,
  font.x = 9,
  font.y = 9,
  font.tickslab = 7,
  font.legend = 8
)

# Save
pdf("figure10_survival.pdf", width = 5, height = 6)
print(p)
dev.off()

png("figure10_survival.png", width = 5, height = 6, 
    units = "in", res = 300)
print(p)
dev.off()
```
