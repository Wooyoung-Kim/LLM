# Save Publication-Quality Figures in R

library(ggplot2)

#' Save a figure in multiple formats for publication
#'
#' @param plot A ggplot2 object
#' @param filename Base filename (without extension)
#' @param width Width in inches
#' @param height Height in inches
#' @param formats Vector of formats to save (default: c("pdf", "png"))
#' @param dpi DPI for raster formats (default: 300)
#' @param path Output directory (default: current directory)
#'
#' @export
save_publication_figure <- function(plot,
                                     filename,
                                     width = 7,
                                     height = 5,
                                     formats = c("pdf", "png"),
                                     dpi = 300,
                                     path = ".") {
  
  # Ensure plot is valid
  if (!inherits(plot, "gg") && !inherits(plot, "patchwork")) {
    stop("plot must be a ggplot2 or patchwork object")
  }
  
  # Save in each requested format
  for (format in formats) {
    output_file <- file.path(path, paste0(filename, ".", format))
    
    message(sprintf("Saving %s...", output_file))
    
    switch(format,
      "pdf" = {
        ggsave(output_file, plot = plot,
               width = width, height = height,
               units = "in", device = "pdf",
               useDingbats = FALSE)  # For compatibility
      },
      "png" = {
        ggsave(output_file, plot = plot,
               width = width, height = height,
               units = "in", device = "png",
               dpi = dpi, bg = "white")
      },
      "tiff" = {
        ggsave(output_file, plot = plot,
               width = width, height = height,
               units = "in", device = "tiff",
               dpi = dpi, compression = "lzw")
      },
      "eps" = {
        ggsave(output_file, plot = plot,
               width = width, height = height,
               units = "in", device = "eps")
      },
      "svg" = {
        ggsave(output_file, plot = plot,
               width = width, height = height,
               units = "in", device = "svg")
      },
      {
        warning(sprintf("Unknown format: %s", format))
      }
    )
  }
  
  invisible(NULL)
}

#' Save figure according to journal specifications
#'
#' @param plot A ggplot2 object
#' @param filename Base filename (without extension)
#' @param journal Journal name ("nature", "science", "cell", "plos")
#' @param figure_type Figure type ("single", "double", "combination")
#' @param path Output directory (default: current directory)
#'
#' @export
save_for_journal <- function(plot,
                              filename,
                              journal = "nature",
                              figure_type = "single",
                              path = ".") {
  
  # Journal-specific specifications
  specs <- list(
    nature = list(
      single = list(width = 89/25.4, height = 89/25.4, dpi = 300),
      double = list(width = 183/25.4, height = 183/25.4, dpi = 300),
      combination = list(width = 89/25.4, height = 120/25.4, dpi = 600),
      formats = c("pdf", "tiff")
    ),
    science = list(
      single = list(width = 55/25.4, height = 55/25.4, dpi = 300),
      double = list(width = 175/25.4, height = 175/25.4, dpi = 300),
      combination = list(width = 55/25.4, height = 80/25.4, dpi = 600),
      formats = c("pdf", "png")
    ),
    cell = list(
      single = list(width = 85/25.4, height = 85/25.4, dpi = 300),
      double = list(width = 178/25.4, height = 178/25.4, dpi = 300),
      combination = list(width = 85/25.4, height = 115/25.4, dpi = 600),
      formats = c("pdf", "tiff")
    ),
    plos = list(
      single = list(width = 83/25.4, height = 83/25.4, dpi = 300),
      double = list(width = 173/25.4, height = 173/25.4, dpi = 300),
      combination = list(width = 83/25.4, height = 115/25.4, dpi = 600),
      formats = c("pdf", "png", "tiff")
    )
  )
  
  # Get specifications
  journal <- tolower(journal)
  if (!journal %in% names(specs)) {
    stop(sprintf("Unknown journal: %s. Choose from: %s",
                 journal, paste(names(specs), collapse = ", ")))
  }
  
  if (!figure_type %in% c("single", "double", "combination")) {
    stop("figure_type must be 'single', 'double', or 'combination'")
  }
  
  spec <- specs[[journal]][[figure_type]]
  formats <- specs[[journal]]$formats
  
  # Save figure
  save_publication_figure(
    plot = plot,
    filename = filename,
    width = spec$width,
    height = spec$height,
    formats = formats,
    dpi = spec$dpi,
    path = path
  )
  
  message(sprintf("\nSaved %s for %s journal (%s column)",
                  filename, journal, figure_type))
}

#' Check if figure size meets journal requirements
#'
#' @param plot A ggplot2 object or NULL (checks last plot)
#' @param journal Journal name
#' @param figure_type Figure type ("single", "double")
#'
#' @export
check_figure_size <- function(plot = NULL, 
                               journal = "nature",
                               figure_type = "single") {
  
  if (is.null(plot)) {
    plot <- ggplot2::last_plot()
  }
  
  # Journal requirements (in mm)
  requirements <- list(
    nature = list(single = 89, double = 183),
    science = list(single = 55, double = 175),
    cell = list(single = 85, double = 178),
    plos = list(single = 83, double = 173)
  )
  
  journal <- tolower(journal)
  if (!journal %in% names(requirements)) {
    stop(sprintf("Unknown journal: %s", journal))
  }
  
  required_width_mm <- requirements[[journal]][[figure_type]]
  required_width_in <- required_width_mm / 25.4
  
  # Build plot to get dimensions
  built <- ggplot2::ggplot_build(plot)
  
  message(sprintf("\n%s %s column requirements:", 
                  toupper(journal), figure_type))
  message(sprintf("  Required width: %.1f mm (%.2f inches)",
                  required_width_mm, required_width_in))
  message("\nSet figure width using:")
  message(sprintf("  ggsave(..., width = %.2f, height = your_height, units = 'in')",
                  required_width_in))
}
