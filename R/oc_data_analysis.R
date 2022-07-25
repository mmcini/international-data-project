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
unbinned_oc_data <- read_excel("data/_unbinned_oc_data.xlsx", na = "NA")
unbinned_oc_spectra <- unbinned_oc_data %>%
                       select(c("355":"2500"))

## Creating bins between 360 and 2498 and adding columns 350 and 2500
## This is the format of the rest of the data
binned_oc_spectra <- binning(unbinned_oc_spectra, bin.size = 10) %>%
                     as_tibble()
binned_oc_spectra <- binned_oc_spectra[-215] %>% # removing band 2498 and replacing for 2500
                     add_column("350" = unbinned_oc_data$`350`, .before = "360") %>%
                     add_column("2500" = unbinned_oc_data$`2500`)

write_csv(binned_oc_spectra, "data/_binned_oc_data.csv")

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
        as_tibble() %>%
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
      as_tibble() %>%
      select(property)
    x_axis <- seq(350, 2500, 10)
    plots[[1]] <- ggplot(cor_matrix[-1, ], aes(x = x_axis, y = "")) +
                  geom_tile(aes_string(fill = property)) + ylab("All countries") +
                  cor_visnir_layout
  }
  return(plots)
}
## Prediction plots
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

raw_oc_data <- read_excel("data/oc_data.xlsx", na = "NA")

## All countries
desc_stats_allcountries <- raw_oc_data %>%
                           select("OC", "K":"Pb") %>%
                           descriptive_stats()
write_excel_csv(desc_stats_allcountries, "tables/descriptive_stats_all_countries.csv")

## By country
desc_stats_bycountries <- raw_oc_data %>%
                          select("country", "OC", "K":"Pb") %>%
                          descriptive_stats(group_by = "country")
write_excel_csv(desc_stats_bycountries, "tables/descriptive_stats_by_country.csv")

## OC plots
ggplot(raw_oc_data, aes(x = OC, fill = country)) +
  geom_histogram() + hist_layout + ggtitle("OC Data - All Countries")
ggsave("figures/pxrf_oc_histogram.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Boxplots
create_boxplots(raw_oc_data, c("Fe", "Ca", "K"), by_country = F)
create_boxplots(raw_oc_data, c("country","Ti"), by_country = T)
create_boxplots(raw_oc_data, c("country","OC"), by_country = T)

## PXRF plots
pxrf_data <- raw_oc_data %>%
             select(country, c(K:Pb))
pxrf_hists <- create_histograms(pxrf_data, c(2:17), "country")
ggarrange(plotlist = pxrf_hists, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")
ggsave("figures/pxrf_variables_histograms.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

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
               geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra")
visnir_labels(visnir_plot)
ggsave("figures/visnir_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

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
                  geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra - Continuum Removal")
visnir_labels(cr_visnir_plot)
ggsave("figures/visnir_cr_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

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
               geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra - First Derivative")
ggsave("figures/visnir_deriv_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

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
png("figures/pxrf_corrplot_all_countries.png", units = "mm", res = 300,
    width = 200, height = 200, bg = "white")
par(family = "Times New Roman")
corrplot::corrplot(cor_pxrf, method = "color", order = "hclust",
                   p.mat = cor_pxrf_pvalue$p, addrect = 2,
                   tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
                   insig = "label_sig", title = "PXRF vs. OC", mar = c(0, 0, 1.5, 0))
dev.off()

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
ggsave("figures/visnir_corrplot_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

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
                          preProcess = preprocess_pls, trControl = control_pls)

## PLS results
visnir_pls_pred <- visnir_pls_model$pred %>%
                   filter(ncomp == 3) %>% # best model with comps = 3
                   arrange(rowIndex) %>%
                   add_column(country = countries$country)
rmse_pls <- min(visnir_pls_model$results$RMSE)
rmse_pls_text <- paste("RMSE: ", round(rmse_pls, 2))
r2_pls <- max(visnir_pls_model$results$Rsquared)
r2_pls_text <- paste("R2: ", round(r2_pls, 2))
ggplot(visnir_pls_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(visnir_pls_pred, r2_pls_text, rmse_pls_text) +
       ggtitle("Vis-NIR - PLS")
ggsave("figures/visnir_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## PLS coefficients
pls_coefs <- visnir_pls_model$finalModel$coefficients %>%
             as_tibble() %>%
             rename("Component 1" = contains("1"),
                    "Component 2" = contains("2"),
                    "Component 3" = contains("3")) %>%
             add_column(bands = seq(350, 2500, 10)) %>%
             pivot_longer(-c("bands"), names_to = "comps", values_to = "values")
ggplot(pls_coefs, aes(x = bands, y = values, color = comps, linetype = comps)) +
       geom_line() + visnir_layout + ylab("Coefficients") +
       ggtitle("Vis-NIR - PLS Coefficients")
ggsave("figures/visnir_pls_coefs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## PLS importance
visnir_pls_importance <- varImp(visnir_pls_model)$importance %>%
                         rownames_to_column("variables") %>%
                         slice_max(n = 20, order_by = Overall) %>% # 10 biggest values
                         arrange(Overall) %>%
                         mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                         mutate(variables = factor(variables,levels = variables))
ggplot(visnir_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - PLS Variable Importance")
ggsave("figures/visnir_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("zv", "center", "scale")
set.seed(100)
visnir_rf_model <- train(OC ~ ., data = visnir_data_allcountries, method = "rf",
                          preProcess = preprocess_rf, trControl = control_rf)

## RF
visnir_rf_pred <- visnir_rf_model$pred %>%
                  filter(mtry == 109) %>%
                  add_column(country = countries$country)
rmse_rf <- min(visnir_rf_model$results$RMSE)
rmse_rf_text <- paste("RMSE: ", round(rmse_rf, 2))
r2_rf <- max(visnir_rf_model$results$Rsquared)
r2_rf_text <- paste("R2: ", round(r2_rf, 2))
ggplot(visnir_rf_pred, aes(x = obs, y = pred, color = country)) +
      prediction_plot_layout +
      annotate_valid_scores(visnir_rf_pred, r2_rf_text, rmse_rf_text) +
      ggtitle("Vis-NIR - RF")
ggsave("figures/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 20, order_by = Overall) %>% # 20 biggest values
                        arrange(Overall) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
      ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("zv", "center", "scale")
set.seed(100)
visnir_cubist_model <- train(OC ~ ., data = visnir_data_allcountries, method = "cubist",
                          preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results
visnir_cubist_pred <- visnir_cubist_model$pred %>%
                  filter(committees == 20, neighbors == 0) %>% # best scores
                  add_column(country = countries$country)
rmse_cubist <- min(visnir_cubist_model$results$RMSE)
rmse_cubist_text <- paste("RMSE: ", round(rmse_cubist, 2))
r2_cubist <- max(visnir_cubist_model$results$Rsquared)
r2_cubist_text <- paste("R2: ", round(r2_cubist, 2))
ggplot(visnir_cubist_pred, aes(x = obs, y = pred, color = country)) +
      prediction_plot_layout +
      annotate_valid_scores(visnir_cubist_pred, r2_cubist_text, rmse_cubist_text) +
      ggtitle("Vis-NIR - Cubist")
ggsave("figures/visnir_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist importance
visnir_cubist_importance <- varImp(visnir_cubist_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - Cubist Variable Importance")
ggsave("figures/visnir_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), n = numeric(), Model = character(),
                       RMSE = numeric(), R2 = numeric()) %>%
                add_row(Dataset = "Vis-NIR", n = nrow(visnir_data_allcountries),
                        Model = "PLS", RMSE = rmse_pls, R2 = r2_pls) %>%
                add_row(Dataset = "Vis-NIR", n = nrow(visnir_data_allcountries),
                        Model = "RF", RMSE = rmse_rf, R2 = r2_rf) %>%
                add_row(Dataset = "Vis-NIR", n = nrow(visnir_data_allcountries),
                        Model = "Cubist", RMSE = rmse_cubist, R2 = r2_cubist)

## PXRF models ##################################################################################
pxrf_data_allcountries <- raw_oc_data %>%
                          select(OC, "K":"Pb") %>%
                          mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
                          drop_na()
countries <- raw_oc_data %>% # keeping countries to plot later
             select(country, OC, "K":"Pb") %>%
             mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # trat NAs as 0
             drop_na() %>%
             select(country)

## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("zv", "center", "scale")
set.seed(100)
pxrf_pls_model <- train(OC ~ ., data = pxrf_data_allcountries, method = "pls",
                          preProcess = preprocess_pls, trControl = control_pls)

## PLS results
pxrf_pls_pred <- pxrf_pls_model$pred %>%
                 filter(ncomp == 3) %>% # best model with comps = 3
                 arrange(rowIndex) %>%
                 add_column(country = countries$country)
rmse_pls <- min(pxrf_pls_model$results$RMSE)
rmse_pls_text <- paste("RMSE: ", round(rmse_pls, 2))
r2_pls <- max(pxrf_pls_model$results$Rsquared)
r2_pls_text <- paste("R2: ", round(r2_pls, 2))
ggplot(pxrf_pls_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pxrf_pls_pred, r2_pls_text, rmse_pls_text) +
       ggtitle("PXRF - PLS")
ggsave("figures/pxrf_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## PLS importance
pxrf_pls_importance <- varImp(pxrf_pls_model)$importance %>%
                       rownames_to_column("variables") %>%
                       slice_max(n = 20, order_by = Overall) %>% # 10 biggest values
                       arrange(Overall) %>%
                       mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - PLS Variable Importance")
ggsave("figures/pxrf_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("zv", "center", "scale")
set.seed(100)
pxrf_rf_model <- train(OC ~ ., data = pxrf_data_allcountries, method = "rf",
                       preProcess = preprocess_rf, trControl = control_rf)

## RF
pxrf_rf_pred <- pxrf_rf_model$pred %>%
                filter(mtry == 9) %>% # parameter with best results
                add_column(country = countries$country)
rmse_rf <- min(pxrf_rf_model$results$RMSE)
rmse_rf_text <- paste("RMSE: ", round(rmse_rf, 2))
r2_rf <- max(pxrf_rf_model$results$Rsquared)
r2_rf_text <- paste("R2: ", round(r2_rf, 2))
ggplot(pxrf_rf_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pxrf_rf_pred, r2_rf_text, rmse_rf_text) +
       ggtitle("PXRF - RF")
ggsave("figures/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
                      rownames_to_column("variables") %>%
                      arrange(Overall) %>%
                      mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - RF Variable Importance")
ggsave("figures/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("zv", "center", "scale")
set.seed(100)
pxrf_cubist_model <- train(OC ~ ., data = pxrf_data_allcountries, method = "cubist",
                           preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results
pxrf_cubist_pred <- pxrf_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 9) %>% # best scores
                    add_column(country = countries$country)
rmse_cubist <- min(pxrf_cubist_model$results$RMSE)
rmse_cubist_text <- paste("RMSE: ", round(rmse_cubist, 2))
r2_cubist <- max(pxrf_cubist_model$results$Rsquared)
r2_cubist_text <- paste("R2: ", round(r2_cubist, 2))
ggplot(pxrf_cubist_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pxrf_cubist_pred, r2_cubist_text, rmse_cubist_text) +
       ggtitle("PXRF - Cubist")
ggsave("figures/pxrf_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist importance
pxrf_cubist_importance <- varImp(pxrf_cubist_model)$importance %>%
                          rownames_to_column("variables") %>%
                          slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                          arrange(Overall) %>%
                          mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - Cubist Variable Importance")
ggsave("figures/pxrf_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
                add_row(Dataset = "PXRF", n = nrow(pxrf_data_allcountries),
                        Model = "PLS", RMSE = rmse_pls, R2 = r2_pls) %>%
                add_row(Dataset = "PXRF", n = nrow(pxrf_data_allcountries),
                        Model = "RF", RMSE = rmse_rf, R2 = r2_rf) %>%
                add_row(Dataset = "PXRF", n = nrow(pxrf_data_allcountries),
                        Model = "Cubist", RMSE = rmse_cubist, R2 = r2_cubist)

## PXRF + Vis-NIR (pv) models ######################################################################
pv_data_allcountries <- raw_oc_data %>%
                        select(OC, "K":"Pb", "350":"2500") %>%
                        mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
                        drop_na()
countries <- raw_oc_data %>% # keeping countries to plot later
             select(country, OC, "K":"Pb", "350":"2500") %>%
             mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
             drop_na() %>%
             select(country)

## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("zv", "center", "scale")
set.seed(100)
pv_pls_model <- train(OC ~ ., data = pv_data_allcountries, method = "pls",
                      preProcess = preprocess_pls, trControl = control_pls)

## PLS results
pv_pls_pred <- pv_pls_model$pred %>%
               filter(ncomp == 3) %>% # best model with comps = 3
               arrange(rowIndex) %>%
               add_column(country = countries$country)
rmse_pls <- min(pv_pls_model$results$RMSE)
rmse_pls_text <- paste("RMSE: ", round(rmse_pls, 2))
r2_pls <- max(pv_pls_model$results$Rsquared)
r2_pls_text <- paste("R2: ", round(r2_pls, 2))
ggplot(pv_pls_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pv_pls_pred, r2_pls_text, rmse_pls_text) +
       ggtitle("PXRF + Vis-NIR - PLS")
ggsave("figures/pv_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## PLS importance
pv_pls_importance <- varImp(pv_pls_model)$importance %>%
                     rownames_to_column("variables") %>%
                     slice_max(n = 20, order_by = Overall) %>% # 20 biggest values
                     arrange(Overall) %>%
                     mutate(variables = factor(variables,levels = variables))
ggplot(pv_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - PLS Variable Importance")
ggsave("figures/pv_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("zv", "center", "scale")
set.seed(100)
pv_rf_model <- train(OC ~ ., data = pv_data_allcountries, method = "rf",
                     preProcess = preprocess_rf, trControl = control_rf)

## RF
pv_rf_pred <- pv_rf_model$pred %>%
              filter(mtry == 2) %>% # parameter with best results
              add_column(country = countries$country)
rmse_rf <- min(pv_rf_model$results$RMSE)
rmse_rf_text <- paste("RMSE: ", round(rmse_rf, 2))
r2_rf <- max(pv_rf_model$results$Rsquared)
r2_rf_text <- paste("R2: ", round(r2_rf, 2))
ggplot(pv_rf_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pv_rf_pred, r2_rf_text, rmse_rf_text) +
       ggtitle("PXRF + Vis-NIR - RF")
ggsave("figures/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
                    rownames_to_column("variables") %>%
                    slice_max(n = 20, order_by = Overall) %>% # 20 biggest values
                    arrange(Overall) %>%
                    mutate(variables = factor(variables,levels = variables))
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("zv", "center", "scale")
set.seed(100)
pv_cubist_model <- train(OC ~ ., data = pv_data_allcountries, method = "cubist",
                         preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results
pv_cubist_pred <- pv_cubist_model$pred %>%
                  filter(committees == 20, neighbors == 0) %>% # best scores
                  add_column(country = countries$country)
rmse_cubist <- min(pv_cubist_model$results$RMSE)
rmse_cubist_text <- paste("RMSE: ", round(rmse_cubist, 2))
r2_cubist <- max(pv_cubist_model$results$Rsquared)
r2_cubist_text <- paste("R2: ", round(r2_cubist, 2))
ggplot(pv_cubist_pred, aes(x = obs, y = pred, color = country)) +
       prediction_plot_layout +
       annotate_valid_scores(pv_cubist_pred, r2_cubist_text, rmse_cubist_text) +
       ggtitle("PXRF + Vis-NIR - Cubist")
ggsave("figures/pv_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist importance
pv_cubist_importance <- varImp(pv_cubist_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                        arrange(Overall) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(pv_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - Cubist Variable Importance")
ggsave("figures/pv_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
                add_row(Dataset = "PXRF + Vis-NIR", n = nrow(pv_data_allcountries),
                        Model = "PLS", RMSE = rmse_pls, R2 = r2_pls) %>%
                add_row(Dataset = "PXRF + Vis-NIR", n = nrow(pv_data_allcountries),
                        Model = "RF", RMSE = rmse_rf, R2 = r2_rf) %>%
                add_row(Dataset = "PXRF + Vis-NIR", n = nrow(pv_data_allcountries),
                        Model = "Cubist", RMSE = rmse_cubist, R2 = r2_cubist)
write_excel_csv(model_scores, "tables/model_scores.csv")
