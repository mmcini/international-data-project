# Libs #############################################################################################

library(RColorBrewer)
library(tidyverse)
library(prospectr)
library(extrafont)
library(corrplot)
library(stringr)
library(ggpubr)
library(readxl)
library(caret)

# Merging Vis-NIR data/fixing format　##############################################################

## Some data from Brazil has 2151 rows
## but most data was binned in intervals of 10 (350, 360 ... 2500)
unbinned_oc_data <- read_excel("../data/_unbinned_oc_data.xlsx", na = "NA")
unbinned_oc_spectra <- unbinned_oc_data %>%
                       select(c("355":"2500"))

## Creating bins between 360 and 2498 and adding columns 350 and 2500
## This is the format of the rest of the data
binned_oc_spectra <- binning(unbinned_oc_spectra, bin.size = 10) %>%
                     as_tibble()
binned_oc_spectra <- binned_oc_spectra[-215] %>% # removing band 2498 and replacing for 2500
                     add_column("350" = unbinned_oc_data$`350`, .before = "360") %>%
                     add_column("2500" = unbinned_oc_data$`2500`)

write_csv(binned_oc_spectra, "../data/_binned_oc_data.csv")

# Functions ########################################################################################

## Functions not related to plots
## Calculates CV
cv <- function(x, na.rm) {
  return(sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm))
}
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
boxplot_layout <- list(theme_bw(),theme(text = element_text(family = "Times New Roman")),
                       xlab(""), ylab(""))
create_boxplots <- function(data, variables, by_country = T) {
  if (by_country) {
    boxplot_data <- data %>%
                    select(all_of(variables)) %>%
                    pivot_longer(-c(country), names_to = "variables", values_to = "values")
    plot <- ggplot(boxplot_data, aes(y = variables, x = values)) +
            geom_boxplot() +
            facet_grid(country ~ variables, scales = "free") +
            boxplot_layout
  } else {
    boxplot_data <- data %>%
                    select(all_of(variables)) %>%
                    pivot_longer(everything(), names_to = "variables", values_to = "values")
    plot <- ggplot(boxplot_data, aes(y = variables, x = values)) +
            geom_boxplot() +
            facet_grid( ~ variables, scales = "free") +
            boxplot_layout
  }
  return(plot)
}
## Vis-NIR
visnir_layout <- list(theme_bw(), ylab("Reflectance factor"), xlab("Wavelength (nm)"),
                      scale_x_continuous(breaks = seq(0, 2500, 250)),
                      theme(text = element_text(family = "Times New Roman"),
                            legend.title = element_blank()))
## Function to add labels to Vis-NIR plots
features <- c(420, 480, 650, 1415, 1930, 2205, 2265, 2350, 2385)
feature_names <- c(as.character(features[1:7]), "2350, 2385", "")
positions <- rep(0, 9)
visnir_labels <- function(plot) {
  for (i in seq_len(length(features))) {
    plot <- plot +
      geom_vline(linetype = 2, xintercept = features[i]) +
      annotate("text", label = feature_names[i],
               x = features[i] - 35, y = positions[i], angle = 90,
               hjust = 0,
               size = 3,
               family = "Times New Roman")
  }
  return(plot)
}
cor_visnir_layout <- list(xlab("Wavelength (nm)"),
                          scale_x_continuous(breaks = seq(350, 2500, 150)),
                          scale_fill_distiller(limits = c(-1, 1), type = "div",
                                               palette = "RdBu", aesthetics = "fill"),
                          theme_bw(),
                          theme(text = element_text(family = "Times New Roman"),
                                plot.title = element_text(hjust = 0.5),
                                panel.grid = element_blank(),
                                axis.text.x = element_text(angle = 90),
                                legend.title = element_blank()))
create_cor_visnir <- function(data = NULL, bands = NULL, property = NULL, group_by = NULL) {
  plots <- list()
  if (!is.null(group_by)) {
    groups <- unique(data[group_by])
    for (i in seq_len(nrow(groups))) {
      group <- groups[i, , drop = T]
      cor_matrix <- data %>%
        filter(get(group_by) == group) %>%
        select(property, bands) %>%
        drop_na() %>%
        cor() %>%
        as.tibble() %>%
        select(property)
      x_axis <- seq(350, 2500, 10)
      plots[[i]] <- ggplot(cor_matrix[-1, ], aes(x = x_axis, y = "")) +
                    geom_tile(aes_string(fill = property)) + ylab(group) +
                    cor_visnir_layout
    }
  } else {
    cor_matrix <- data %>%
      select(property, bands) %>%
      drop_na() %>%
      cor() %>%
      as.tibble() %>%
      select(property)
    x_axis <- seq(350, 2500, 10)
    plots[[1]] <- ggplot(cor_matrix[-1, ], aes(x = x_axis, y = "")) +
                  geom_tile(aes_string(fill = property)) + ylab("All countries") +
                  cor_visnir_layout
  }
  return(plots)
}
prediction_plot_layout <- list(xlab("Observed OC (%)"), ylab("Predicted OC (%)"),
                               geom_point(), geom_abline(slope = 1), theme_bw(),
                               theme(text = element_text(family = "Times New Roman"),
                                     legend.title = element_blank()))
annotate_valid_scores <- function(data, r2, rmse) {
  ## data must have pred and obs columns for relative positioning
  annotate("text", label = c(r2, rmse), x = min(data$obs),
           y = c(max(data$pred), max(data$pred) - 0.3),
           vjust = 0, hjust = 0, family = "Times New Roman")
}
importance_plot_layout <- list(geom_bar(stat = "identity", width = 0.2),coord_flip(),
                               xlab("Variables"), ylab("Importance (%)"),
                               theme_bw(), theme(text = element_text("Times New Roman")))

# Descriptive stats　###############################################################################

raw_oc_data <- read_excel("../data/oc_data.xlsx", na = "NA")

## All countries
desc_stats_allcountries <- raw_oc_data %>%
                           select("OC", "K":"Pb") %>%
                           descriptive_stats()

## By country
desc_stats_bycountries <- raw_oc_data %>%
                          select("country", "OC", "K":"Pb") %>%
                          descriptive_stats(group_by = "country")

## OC plots
ggplot(raw_oc_data, aes(x = OC, fill = country)) +
  geom_histogram() + hist_layout

## Boxplots
create_boxplots(raw_oc_data, c("Fe", "Ca", "K"), by_country = F)
create_boxplots(raw_oc_data, c("country","Ti"), by_country = T)
create_boxplots(raw_oc_data, c("country","OC"), by_country = T)

## PXRF plots
pxrf_data <- raw_oc_data %>%
             select(country, c(K:Pb))
pxrf_hists <- create_histograms(pxrf_data, c(2:17), "country")
ggarrange(plotlist = pxrf_hists, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")

## Vis-NIR plots
## No continuum removal
visnir_data <- raw_oc_data %>%
               select(country, c("350":"2500")) %>%
               filter(country != "Africa") %>% # Africa has no Vis-NIR data
               group_by(country) %>%
               summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
               pivot_longer(c("350":"2500"), names_to = "wavelength", values_to = "reflectance") %>%
               mutate(wavelength = as.numeric(wavelength))
visnir_plot <- ggplot(visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
               geom_line() + visnir_layout
visnir_labels(visnir_plot)

## With continuum removal
cr_visnir_data <- raw_oc_data %>%
                  select(country, c("350":"2500")) %>%
                  filter(country != "Africa") %>% # Africa has no Vis-NIR data
                  group_by(country) %>%
                  summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
                  select(c("350":"2500")) %>%
                  continuumRemoval() %>%
                  as_tibble()
colnames(cr_visnir_data) <- seq(350, 2500, 10)
cr_visnir_data <- add_column(cr_visnir_data, country = c("Brazil", "France", "India", "US"),
                             .before = "350") %>%
                             pivot_longer(c("350":"2500"), names_to = "wavelength",
                                          values_to = "reflectance") %>%
                             mutate(wavelength = as.numeric(wavelength))
cr_visnir_plot <- ggplot(cr_visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
                  geom_line() + visnir_layout
visnir_labels(cr_visnir_plot)

## Savitzky-Golay + first derivative
deriv_visnir_data <- raw_oc_data %>%
  select(country, c("350":"2500")) %>%
  filter(country != "Africa") %>% # Africa has no Vis-NIR data
  group_by(country) %>%
  summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
  select(c("350":"2500")) %>%
  savitzkyGolay(m = 1, p = 3, w = 11) %>%
  as_tibble()
deriv_visnir_data <- add_column(deriv_visnir_data, country = c("Brazil", "France", "India", "US"),
                               .before = "400") %>%
                     pivot_longer(c("400":"2450"), names_to = "wavelength",
                                   values_to = "reflectance") %>%
                     mutate(wavelength = as.numeric(wavelength))
ggplot(deriv_visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
               geom_line() + visnir_layout

## Correlations
## PXRF correlations
cor_pxrf <- raw_oc_data %>%
            select(OC, c(K:Pb)) %>%
            drop_na() %>%
            cor()
cor_pxrf_pvalue <- raw_oc_data %>%
                   select(OC, c(K:Pb)) %>%
                   drop_na() %>%
                   cor.mtest(conf.level = 0.95) # calculates p-value
par(family = "Times New Roman")
corrplot(cor_pxrf, method = "color", order = "hclust", p.mat = cor_pxrf_pvalue$p, addrect = 2,
         tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
         title = "PXRF vs. OC", mar = c(0, 0, 1.5, 0))

## Vis-NIR correlations
## All countries
cor_visnir_plot_allcountries <- create_cor_visnir(raw_oc_data,
                                                  bands = c(27:242), property = "OC")

## By country
cor_visnir_bycountry <- raw_oc_data %>%
                        filter(country != "Africa") # Africa has no Vis-NIR data
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(27:242),
                                                property = "OC", group_by = "country")

## All Vis-NIR correlation plots
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")

# Modeling　########################################################################################

## All countries ###################################################################################

## Vis-NIR models ##################################################################################
visnir_data_allcountries <- raw_oc_data %>%
                            select(OC, "350":"2500") %>%
                            drop_na()
countries <- raw_oc_data %>% # keeping countries to plot later
            select(country, OC, "350":"2500") %>%
            drop_na() %>%
            select(country)

## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("zv", "center", "scale")
set.seed(100)
visnir_pls_model <- train(OC ~ ., data = visnir_data_allcountries, method = "pls",
                          preProcess = preprocess_pls,trControl = control_pls)

## PLS results
visnir_pls_model # best model with comps = 3
visnir_pls_pred <- visnir_pls_model$pred %>%
                   filter(ncomp == 3) %>%
                   arrange(rowIndex) %>%
                   add_column(country = countries$country)
rmse_pls <- paste("RMSE: ", round(min(visnir_pls_model$results$RMSE), 2))
r2_pls <- paste("R2: ", round(max(visnir_pls_model$results$Rsquared), 2))
ggplot(visnir_pls_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(visnir_pls_pred, r2_pls, rmse_pls)

## PLS coefficients
pls_coefs <- visnir_pls_model$finalModel$coefficients %>%
             as_tibble() %>%
             rename("Component 1" = contains("1"),
                    "Component 2" = contains("2"),
                    "Component 3" = contains("3")) %>%
             add_column(bands = seq(350, 2500, 10)) %>%
             pivot_longer(-c("bands"), names_to = "comps", values_to = "values")
ggplot(pls_coefs, aes(x = bands, y = values, color = comps, linetype = comps)) +
       geom_line() + visnir_layout + ylab("Coefficients")

## PLS importance
visnir_pls_importance <- varImp(visnir_pls_model)$importance %>%
                         rownames_to_column("variables") %>%
                         slice_max(n = 20, order_by = Overall) %>% # 10 biggest values
                         arrange(Overall) %>%
                         mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                         mutate(variables = factor(variables,levels = variables))
ggplot(visnir_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout

## RF with PCA preprocessing
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T,
                           preProcOptions = list(pcaComp = 10)) # use 10 first PCs
preprocess_rf <- c("zv", "center", "scale", "pca")
set.seed(100)
visnir_rf_model <- train(OC ~ ., data = visnir_data_allcountries, method = "rf",
                          preProcess = preprocess_rf,trControl = control_rf)
visnir_rf_model
## RF with PCA preprocessing results
visnir_rf_pred <- visnir_rf_model$pred %>%
                  filter(mtry == 2) %>%
                  add_column(country = countries$country)
rmse_rf <- paste("RMSE: ", round(min(visnir_rf_model$results$RMSE), 2))
r2_rf <- paste("R2: ", round(max(visnir_rf_model$results$Rsquared), 2))
ggplot(visnir_rf_pred, aes(x = obs, y = pred, color = country)) +
      prediction_plot_layout +
      annotate_valid_scores(visnir_rf_pred, r2_rf, rmse_rf)

visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
                        rownames_to_column("variables") %>%
                        arrange(Overall) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("zv", "center", "scale")
set.seed(100)
visnir_cubist_model <- train(OC ~ ., data = visnir_data_allcountries, method = "cubist",
                          preProcess = preprocess_cubist,trControl = control_cubist)

## Cubist results
visnir_cubist_pred <- visnir_cubist_model$pred %>%
                  filter(committees == 20, neighbors == 0) %>% # best scores
                  add_column(country = countries$country)
rmse_cubist <- paste("RMSE: ", round(min(visnir_cubist_model$results$RMSE), 2))
r2_cubist <- paste("R2: ", round(max(visnir_cubist_model$results$Rsquared), 2))
ggplot(visnir_cubist_pred, aes(x = obs, y = pred, color = country)) +
      prediction_plot_layout +
      annotate_valid_scores(visnir_cubist_pred, r2_cubist, rmse_cubist)

## Cubist importance
visnir_cubist_importance <- varImp(visnir_cubist_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout
