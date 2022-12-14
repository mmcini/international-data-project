# Libs #############################################################################################

library(rnaturalearthdata)
library(rnaturalearth)
library(RColorBrewer)
library(conflicted)
library(tidyverse)
library(patchwork)
library(prospectr)
library(ggspatial)
library(extrafont)
library(corrplot)
library(tidytext)
library(stringr)
library(geodata)
library(ggrepel)
library(ggtern)
library(ggpubr)
library(readxl)
library(raster)
library(caret)
library(tune)
library(sp)
library(sf)

sf_use_s2(FALSE) # turns off s2 processing, doesnt check geometry

## avoids masking conflincs
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("theme_bw", "ggplot2")
conflict_prefer("ggsave", "ggplot2")
conflict_prefer("aes", "ggplot2")
conflict_prefer("annotate", "ggplot2")

# Functions ########################################################################################

## Functions not related to plots
## Calculates CV
cv <- function(x, na.rm) {
  return(sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
}

## Extracts all variables in the importance table so they can be counted
extract_pxrf_vars <- function(variables) {
  all_elements <- c()
  for (i in seq_len(length(variables))) {
    ## tokenizes the strings in the variables column and appends to a single vector
    tokens <- unlist(strsplit(str_remove_all(variables[i], ","), "\\s+"))
    all_elements <- c(all_elements, tokens)
  }
  return(all_elements)
}

## Normalizes data based in min-max
range_preprocess <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

## Automatic descriptive stats
descriptive_stats <- function(data, group_by = NULL) {
functions <- c(max, min, mean, median, sd, cv, IQR)
names <- c("variables", "max", "min", "mean", "median", "sd", "cv", "IQR")
  if (!is.null(group_by)) {
    names <- c(group_by, names)
    stats <- data %>%
             pivot_longer(-c(group_by), names_to = "variables", values_to = "values") %>%
             group_by_at(c(group_by, "variables")) %>%
             summarize(across(everything(), functions, na.rm = T)) %>%
             rename_with(~ names)
  } else {
    stats <- data %>%
             pivot_longer(everything(), names_to = "variables", values_to = "values") %>%
             group_by(variables) %>%
             summarize(across(everything(), functions, na.rm = T)) %>%
             rename_with(~ names)
  }
return(stats)
}

# Plot layouts #####################################################################################
## Lengthy and repetitive plot formatting and functions are defined here

### Country colors
country_colors <- brewer.pal(5, "Set1")
names(country_colors) <- c("Mozambique", "Brazil", "France", "India", "US")

## histograms
hist_layout <- list(theme_bw(), ylab("Count"),
                    theme(text = element_text(family = "Times New Roman"),
                          legend.title = element_blank()))

create_histograms <- function(data = NULL, variables = NULL, group_by = NULL) {
  if (is.null(variables)) {
    variables <- as.character(colnames(data))
  } else {
    variables <- data %>%
                select(all_of(variables)) %>%
                colnames()
  }
  plots <- list()
  for (i in seq_len(length(variables))) {
    plots[[i]] <- ggplot(data, aes_string(x = variables[i], fill = group_by)) +
                  geom_histogram() + hist_layout
  }
  return(plots)
}

## Boxplots
boxplot_layout <- list(theme_bw(),
                       theme(text = element_text(family = "Times New Roman"),
                             legend.title = element_blank()),
                       xlab(""), ylab(""))
create_boxplots <- function(data, variables, by_country = T,
                            ncol = NULL, nrow = NULL, ordered = F, levels = NULL) {
  if (by_country) {
    boxplot_data <- data %>%
                    select(all_of(variables)) %>%
                    pivot_longer(-c(country), names_to = "variables", values_to = "values")
    if (ordered) {
      boxplot_data <- boxplot_data %>%
                      mutate(variables = factor(variables,
                                                levels = levels,
                                                ordered = T))
    }
    plot <- ggplot(boxplot_data, aes(y = variables, x = values, fill = country)) +
            geom_boxplot() +
            scale_fill_manual(drop = T, limits = levels(boxplot_data$country),
                              values = country_colors) +
            facet_wrap( ~ variables, scales = "free", ncol = ncol, nrow = nrow) +
            boxplot_layout
  } else {
    boxplot_data <- data %>%
                    select(all_of(variables)) %>%
                    pivot_longer(everything(), names_to = "variables", values_to = "values")
    plot <- ggplot(boxplot_data, aes(y = variables, x = values)) +
            geom_boxplot() +
            facet_wrap( ~ variables, scales = "free") +
            boxplot_layout
  }
  return(plot)
}

## Vis-NIR
visnir_layout <- list(theme_bw(), ylab("Reflectance"), xlab("Wavelength (nm)"),
                      scale_x_continuous(breaks = seq(0, 2500, 250)),
                      scale_color_manual(limits = c("Brazil", "France", "India", "US"),
                                         values = country_colors),
                      theme(text = element_text(family = "Times New Roman"),
                            legend.title = element_blank()))
## Function to add labels to Vis-NIR plots
features_pos <- c(420, 480, 650, 1415, 1930, 2205)
features_names <- c("Gt", "Gt", "Hm", "2:1/1:1", "2:1/Water", "Kt")
positions <- rep(0, 6)
visnir_labels <- function(plot) {
  for (i in seq_len(length(features_names))) {
    plot <- plot +
      geom_vline(linetype = 2, xintercept = features_pos[i]) +
      annotate("text", label = features_names[i],
               x = features_pos[i] - 35, y = positions[i], angle = 90,
               hjust = 0,
               size = 3,
               family = "Times New Roman")
  }
  return(plot)
}

visnir_cor_labels <- function(plot) {
  for (i in seq_len(length(features_names))) {
    plot <- plot +
      geom_vline(linetype = 2, xintercept = features_pos[i]) +
      annotate("text", label = features_names[i],
               x = features_pos[i] - 35, y = positions[i] + 0.6, angle = 90,
               hjust = 0,
               size = 3,
               family = "Times New Roman")
  }
  return(plot)
}

cor_visnir_layout <- list(xlab(""), ylab(""),
                          scale_x_continuous(breaks = seq(0, 2500, 250)),
                          scale_fill_distiller(limits = c(-1, 1), type = "div",
                                               palette = "RdBu", aesthetics = "fill"),
                          theme_bw(),
                          theme(text = element_text(family = "Times New Roman"),
                                panel.grid = element_blank(),
                                legend.title = element_blank()))

create_cor_visnir <- function(data = NULL, bands = NULL, property = NULL, group_by = NULL) {
  plots <- list()
  if (!is.null(group_by)) {
    groups <- unique(data[group_by])
    n <- seq_len(nrow(groups))
    for (i in n) {
      group <- groups[i, , drop = T]
      cor_matrix <- data %>%
        filter(get(group_by) == group) %>%
        select(property, bands) %>%
        drop_na() %>%
        cor() %>%
        as_tibble() %>%
        select(property)
      x_axis <- seq(350, 2500, 10)
      plots[[i]] <- ggplot(cor_matrix[-1, ], aes(x = x_axis, y = "")) +
                    geom_tile(aes_string(fill = property)) + ylab(group) +
                    cor_visnir_layout
      if (i == max(n)) {
      plots[[i]] <- plots[[i]] + xlab("Wavelength (cm)")
      }
    }
  } else {
    cor_matrix <- data %>%
      select(property, bands) %>%
      drop_na() %>%
      cor() %>%
      as_tibble() %>%
      select(property)
    x_axis <- seq(350, 2500, 10)
    plots[[1]] <- ggplot(cor_matrix[-1, ], aes(x = x_axis, y = "")) +
                  geom_tile(aes_string(fill = property)) + ylab("All countries") +
                  ggtitle(property) +
                  cor_visnir_layout
  }
  return(plots)
}

annotate_valid_scores <- function(data, r2, rmse, y_coord) {
  ## data must have pred and obs columns for relative positioning
  annotate("text", label = c(r2, rmse), x = min(data$obs),
           y = c(y_coord, y_coord * 0.90),
           vjust = 0, hjust = 0, family = "Times New Roman")
}

annotate_coord <- function(y_limit, x_limit, coord_scale) {
  # gets the the highest value from x and y ranges to use as annotates y coord
  y_coord <- min(y_limit) + (max(y_limit) - min(y_limit)) * coord_scale
  x_coord <- min(x_limit) + (max(x_limit) - min(x_limit)) * coord_scale
  coords <- c(y_coord, x_coord)
  return(max(coords))
}

mean_from_folds <- function(data, type = "RMSE", folds_column = "Resample") {
  folds <- unique(data[[folds_column]])
  values <- c()
  count <- 1
  if (type == "RMSE") {
    for (fold in folds) {
       set <- data %>%
              filter(Resample == fold)
       values[count] <- caret::RMSE(pred = set$pred, obs = set$obs)
       count <- count + 1
    }
  } else if (type == "R2") {
    for (fold in folds) {
       set <- data %>%
              filter(Resample == fold)
       values[count] <- caret::R2(pred = set$pred, obs = set$obs)
       count <- count + 1
    }
  } else { return(print("Type must be RMSE or R2")) }
  return(mean(values))
}

validation_plot <- function(cv_set, valid_set, variable = "", dataset = "",
                            model = "", group_by = NULL, coord_scale = 0.9) {

  prediction_plot_layout <- list(geom_abline(slope = 1), coord_obs_pred(), theme_bw(),
                               theme(text = element_text(family = "Times New Roman"),
                                     legend.title = element_blank()))

  if (!is.null(dev.list())) {dev.off()} # refreshes device if not null
  data <- list(cv_set, valid_set)
  has_group <- !is.null(group_by)
  plots <- list()
  count <- 1
  for (set in data) {
    if (count == 1) { # calculates the mean of each fold for cv dataset (count = 1)
      rmse <- mean_from_folds(set, type = "RMSE")
      r2 <- mean_from_folds(set, type = "R2")
    } else {
      rmse <- RMSE(pred = set$pred, obs = set$obs)
      r2 <- caret::R2(pred = set$pred, obs = set$obs)
    }
    rmse_text <- paste("RMSE: ", round(rmse, 2))
    r2_text <- paste("R2: ", round(r2, 2))
    if (has_group) { # adds colors to points if grouped
      set$country <- droplevels(set$country)
      colors <- scale_color_manual(drop = T, limits = levels(set$country), values = country_colors)
      colored_points <- list(geom_point(aes_string(color = group_by)), colors)
      }
    plots[[count]] <- ggplot(set, aes(x = obs, y = pred)) +
                             {if (has_group) {colored_points}
                              else {geom_point(color = "gray")}} +
                            xlab(paste("Observed", variable, "(%)")) +
                            ylab(paste("Predicted", variable, "(%)")) +
                            geom_smooth(aes(color = null), method = "lm",se = F,
                                        linetype = "dashed", col = "black") +
                            prediction_plot_layout
    coord <- annotate_coord(layer_scales(plots[[1]])$y$get_limits(),
                            layer_scales(plots[[1]])$x$get_limits(),
                            coord_scale)
    plots[[count]] <- plots[[count]] +
                      annotate_valid_scores(set, r2_text, rmse_text, coord)
    count <- count + 1
  }
  plots[[1]] <- plots[[1]] + ggtitle(paste(dataset, "-", model, "(80%)"))
  plots[[2]] <- plots[[2]] + ggtitle(paste(dataset, "-", model, "(20%)"))
  return(plots)
}

importance_plot_layout <- list(geom_bar(stat = "identity", width = 0.2),coord_flip(),
                               xlab("Variables"), ylab("Importance (%)"),
                               theme_bw(), theme(text = element_text("Times New Roman")))
