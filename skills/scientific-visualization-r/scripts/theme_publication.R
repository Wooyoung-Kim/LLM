# Publication-ready theme for ggplot2
# Based on common requirements for scientific journals

library(ggplot2)

#' Publication-ready theme
#'
#' A clean, minimal theme suitable for scientific publications
#'
#' @param base_size Base font size (default: 8 pt)
#' @param base_family Base font family (default: "Arial")
#' @param base_line_size Base line size
#' @param base_rect_size Base rectangle size
#'
#' @return A ggplot2 theme object
#' @export
theme_publication <- function(base_size = 8,
                               base_family = "Arial",
                               base_line_size = 0.5,
                               base_rect_size = 0.5) {
  
  theme_bw(base_size = base_size,
           base_family = base_family,
           base_line_size = base_line_size,
           base_rect_size = base_rect_size) +
    theme(
      # Text elements
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      axis.title = element_text(size = base_size + 1, face = "plain"),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      
      # Remove unnecessary elements
      panel.grid.major = element_line(color = "grey90", size = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      
      # Legend
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(4, "mm"),
      
      # Strip (facet labels)
      strip.background = element_rect(fill = "grey90", color = "black", size = 0.5),
      
      # Margins
      plot.margin = margin(5, 5, 5, 5)
    )
}

#' Nature journal theme
#'
#' Theme following Nature journal specifications
#'
#' @return A ggplot2 theme object
#' @export
theme_nature <- function() {
  theme_publication(base_size = 7, base_family = "Arial") +
    theme(
      # Nature prefers minimal gridlines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Nature style panel
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      
      # Nature uses lowercase panel labels
      plot.tag = element_text(face = "bold", size = 10)
    )
}

#' Science journal theme
#'
#' Theme following Science journal specifications
#'
#' @return A ggplot2 theme object
#' @export
theme_science <- function() {
  theme_publication(base_size = 6, base_family = "Arial") +
    theme(
      # Science prefers very clean plots
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    )
}

#' Cell journal theme
#'
#' Theme following Cell journal specifications
#'
#' @return A ggplot2 theme object
#' @export
theme_cell <- function() {
  theme_publication(base_size = 7, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(color = "grey95", size = 0.25),
      panel.grid.minor = element_blank()
    )
}

#' PLOS journal theme
#'
#' Theme following PLOS journal specifications
#'
#' @return A ggplot2 theme object
#' @export
theme_plos <- function() {
  theme_publication(base_size = 8, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(color = "grey90", size = 0.25),
      panel.grid.minor = element_blank()
    )
}

#' Minimal theme for presentations
#'
#' Larger fonts suitable for posters and presentations
#'
#' @param base_size Base font size (default: 14 pt)
#' @return A ggplot2 theme object
#' @export
theme_presentation <- function(base_size = 14) {
  theme_publication(base_size = base_size) +
    theme(
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.title = element_text(size = base_size + 1),
      plot.title = element_text(size = base_size + 4, face = "bold")
    )
}
