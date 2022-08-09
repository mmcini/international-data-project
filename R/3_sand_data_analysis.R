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
prediction_plot_layout <- list(xlab("Observed Sand (%)"), ylab("Predicted Sand (%)"),
                               geom_point(), geom_abline(slope = 1), theme_bw(),
                               theme(text = element_text(family = "Times New Roman"),
                                     legend.title = element_blank()))
annotate_valid_scores <- function(data, r2, rmse) {
  ## data must have pred and obs columns for relative positioning
  annotate("text", label = c(r2, rmse), x = min(data$obs),
           y = c(max(data$pred), max(data$pred) - 1.3),
           vjust = 0, hjust = 0, family = "Times New Roman")
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
validation_plot <- function(cv_set, valid_set, dataset = "", model = "") {
  data <- list(cv_set, valid_set)
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
    plots[[count]] <- ggplot(set, aes(x = obs, y = pred, color = country)) +
                            prediction_plot_layout +
                            annotate_valid_scores(set, r2_text, rmse_text)
    count <- count + 1
  }
  plots[[1]] <- plots[[1]] + ggtitle(paste(dataset, "-", model, "(Cross-Validation, 80%)"))
  plots[[2]] <- plots[[2]] + ggtitle(paste(dataset, "-", model, "(Hold-out Validation, 20%)"))
  return(plots)
}
importance_plot_layout <- list(geom_bar(stat = "identity", width = 0.2),coord_flip(),
                               xlab("Variables"), ylab("Importance (%)"),
                               theme_bw(), theme(text = element_text("Times New Roman")))

# Descriptive stats　###############################################################################

raw_sand_data <- read_excel("data/oc_texture_data.xlsx", na = "NA")

## All countries
desc_stats_allcountries <- raw_sand_data %>%
                           select("sand", "K":"Pb") %>%
                           descriptive_stats()
write_excel_csv(desc_stats_allcountries, "tables/sand/descriptive_stats_all_countries.csv")

## By country
desc_stats_bycountries <- raw_sand_data %>%
                          select("country", "sand", "K":"Pb") %>%
                          descriptive_stats(group_by = "country")
write_excel_csv(desc_stats_bycountries, "tables/sand/descriptive_stats_by_country.csv")

## sand plots
ggplot(raw_sand_data, aes(x = sand, fill = country)) +
  geom_histogram() + hist_layout + ggtitle("sand Data - All Countries")
ggsave("figures/sand/pxrf_oc_histogram.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Boxplots
create_boxplots(raw_sand_data, c("Fe", "Ca", "K"), by_country = F)
create_boxplots(raw_sand_data, c("country","Ti"), by_country = T)
create_boxplots(raw_sand_data, c("country","sand"), by_country = T)

## PXRF plots
pxrf_data <- raw_sand_data %>%
             select(country, c(K:Pb))
pxrf_hists <- create_histograms(pxrf_data, c(2:17), "country")
ggarrange(plotlist = pxrf_hists, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")
ggsave("figures/sand/pxrf_variables_histograms.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Vis-NIR plots
## No continuum removal
visnir_data <- raw_sand_data %>%
               select(country, c("350":"2500")) %>%
               filter(country != "Africa") %>% # Africa has no Vis-NIR data
               group_by(country) %>%
               summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
               pivot_longer(c("350":"2500"), names_to = "wavelength", values_to = "reflectance") %>%
               mutate(wavelength = as.numeric(wavelength))
visnir_plot <- ggplot(visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
               geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra")
visnir_labels(visnir_plot)
ggsave("figures/sand/visnir_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## With continuum removal
cr_visnir_data <- raw_sand_data %>%
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
ggsave("figures/sand/visnir_cr_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Savitzky-Golay + first derivative
deriv_visnir_data <- raw_sand_data %>%
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
ggsave("figures/sand/visnir_deriv_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Correlations
## PXRF correlations
cor_pxrf <- raw_sand_data %>%
            select(sand, c(K:Pb)) %>%
            drop_na() %>%
            cor()
cor_pxrf_pvalue <- raw_sand_data %>%
                   select(sand, c(K:Pb)) %>%
                   drop_na() %>%
                   cor.mtest(conf.level = 0.95) # calculates p-value
png("figures/sand/pxrf_corrplot_all_countries.png", units = "mm", res = 300,
    width = 200, height = 200, bg = "white")
par(family = "Times New Roman")
corrplot::corrplot(cor_pxrf, method = "color", order = "hclust",
                   p.mat = cor_pxrf_pvalue$p, addrect = 2,
                   tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
                   insig = "label_sig", title = "PXRF vs. sand", mar = c(0, 0, 1.5, 0))
dev.off()

## Vis-NIR correlations
## All countries
cor_visnir_plot_allcountries <- create_cor_visnir(raw_sand_data,
                                                  bands = c(27:242), property = "sand")

## By country
cor_visnir_bycountry <- raw_sand_data %>%
                        filter(country != "Africa") # Africa has no Vis-NIR data
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(27:242),
                                                property = "sand", group_by = "country")

## All Vis-NIR correlation plots
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")
ggsave("figures/sand/visnir_corrplot_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

# Modeling　########################################################################################

## All countries ###################################################################################

## Vis-NIR models ##################################################################################
## Training and validation data
visnir_data_allcountries <- raw_sand_data %>%
                            select(country, sand, "350":"2500") %>%
                            drop_na() %>%
                            # substituting numbers by legal variable names
                            rename_with(~ str_replace(., "^[0-9]+$", paste0("band_",. , "_nm")))

## Partitioning
set.seed(100)
partition_index <- createDataPartition(visnir_data_allcountries$sand, p = 0.8, list = F)
visnir_train_data <- visnir_data_allcountries %>%
                     slice(partition_index)
visnir_valid_data <- visnir_data_allcountries %>%
                     slice(-partition_index)

## Models
## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("nzv", "center", "scale")
set.seed(100)
visnir_pls_model <- train(sand ~ ., data = visnir_train_data[-1], method = "pls",
                          preProcess = preprocess_pls, trControl = control_pls)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
visnir_pls_cv <- visnir_pls_model$pred %>%
                 filter(ncomp == 3) %>% # best model with comps = 3
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_pls_valid <- predict(visnir_pls_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$sand) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_pls_plots <- validation_plot(visnir_pls_cv, visnir_pls_valid, "Vis-NIR", "PLS")
ggarrange(plotlist = visnir_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/visnir_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

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
ggsave("figures/sand/visnir_pls_coefs.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## PLS importance
visnir_pls_importance <- varImp(visnir_pls_model)$importance %>%
                         rownames_to_column("variables") %>%
                         slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                         arrange(Overall) %>%
                         mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                         mutate(variables = factor(variables,levels = variables))
ggplot(visnir_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - PLS Variable Importance")
ggsave("figures/sand/visnir_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("nzv", "center", "scale")
set.seed(100)
visnir_rf_model <- train(sand ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess_rf, trControl = control_rf)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
                filter(mtry == 109) %>% # best model with mtry = 109
                arrange(rowIndex) %>%
                add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = visnir_valid_data$sand) %>%
                   add_column(country = visnir_valid_data$country) %>%
                   rename(pred = value)
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid, "Vis-NIR", "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
      ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/sand/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("nzv", "center", "scale")
set.seed(100)
visnir_cubist_model <- train(sand ~ ., data = visnir_train_data[-1], method = "cubist",
                             preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
visnir_cubist_cv <- visnir_cubist_model$pred %>%
                    filter(committees == 10, neighbors == 0) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = visnir_train_data$country)
visnir_cubist_valid <- predict(visnir_cubist_model, newdata = visnir_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = visnir_valid_data$sand) %>%
                       add_column(country = visnir_valid_data$country) %>%
                       rename(pred = value)
visnir_cubist_plots <- validation_plot(visnir_cubist_cv, visnir_cubist_valid, "Vis-NIR", "Cubist")
ggarrange(plotlist = visnir_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/visnir_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
visnir_cubist_importance <- varImp(visnir_cubist_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - Cubist Variable Importance")
ggsave("figures/sand/visnir_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
                add_row(Dataset = "Vis-NIR", Model = "PLS",
                        n_cv = nrow(visnir_pls_cv), n_valid = nrow(visnir_pls_valid),
                        RMSE_cv = RMSE(visnir_pls_cv$pred, visnir_pls_cv$obs),
                        R2_cv = caret::R2(visnir_pls_cv$pred, visnir_pls_cv$obs),
                        RMSE_valid = RMSE(visnir_pls_valid$pred, visnir_pls_valid$obs),
                        R2_valid = caret::R2(visnir_pls_valid$pred, visnir_pls_valid$obs)) %>%
                add_row(Dataset = "Vis-NIR", Model = "RF",
                        n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
                        RMSE_cv = RMSE(visnir_rf_cv$pred, visnir_rf_cv$obs),
                        R2_cv = caret::R2(visnir_rf_cv$pred, visnir_rf_cv$obs),
                        RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
                        R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs)) %>%
                add_row(Dataset = "Vis-NIR", Model = "Cubist",
                        n_cv = nrow(visnir_cubist_cv), n_valid = nrow(visnir_cubist_valid),
                        RMSE_cv = RMSE(visnir_cubist_cv$pred, visnir_cubist_cv$obs),
                        R2_cv = caret::R2(visnir_cubist_cv$pred, visnir_cubist_cv$obs),
                        RMSE_valid = RMSE(visnir_cubist_valid$pred, visnir_cubist_valid$obs),
                        R2_valid = caret::R2(visnir_cubist_valid$pred, visnir_cubist_valid$obs))

## PXRF models #####################################################################################
pxrf_data_allcountries <- raw_sand_data %>%
                          select(country, sand, "K":"Pb") %>%
                          mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
                          drop_na()

## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_allcountries[-c(1, 2)], y = pxrf_data_allcountries[["sand"]],
                        sizes = c(1:16), rfeControl = control_rfe)
predictors(pxrf_feat_select)
## Selected variables: "Ca" "Sr" "Zr" "K" "Zn" "Mn" "Ti" "Fe" "Rb" "Cu" "Pb" "V" "Cr" "Ni" "Ba" "As"
pxrf_data_allcountries <- pxrf_data_allcountries %>%
                          select(country, sand, predictors(pxrf_feat_select))

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_allcountries$sand, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_allcountries %>%
                   slice(partition_index)
pxrf_valid_data <- pxrf_data_allcountries %>%
                   slice(-partition_index)

## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("nzv", "center", "scale")
set.seed(100)
pxrf_pls_model <- train(sand ~ ., data = pxrf_train_data[-1], method = "pls",
                        preProcess = preprocess_pls, trControl = control_pls)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_pls_cv <- pxrf_pls_model$pred %>%
               filter(ncomp == 3) %>% # best model with comps = 3
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_pls_valid <- predict(pxrf_pls_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$sand) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_pls_plots <- validation_plot(pxrf_pls_cv, pxrf_pls_valid, "PXRF", "PLS")
ggarrange(plotlist = pxrf_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pxrf_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PLS importance
pxrf_pls_importance <- varImp(pxrf_pls_model)$importance %>%
                       rownames_to_column("variables") %>%
                       slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                       arrange(Overall) %>%
                       mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - PLS Variable Importance")
ggsave("figures/sand/pxrf_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("nzv", "center", "scale")
set.seed(100)
pxrf_rf_model <- train(sand ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess_rf, trControl = control_rf)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
              filter(mtry == 9) %>% # best model with mtry = 9
              arrange(rowIndex) %>%
              add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
                 as_tibble() %>%
                 add_column(obs= pxrf_valid_data$sand) %>%
                 add_column(country = pxrf_valid_data$country) %>%
                 rename(pred = value)
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid, "PXRF", "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
                      rownames_to_column("variables") %>%
                      slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                      arrange(Overall) %>%
                      mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - RF Variable Importance")
ggsave("figures/sand/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("nzv", "center", "scale")
set.seed(100)
pxrf_cubist_model <- train(sand ~ ., data = pxrf_train_data[-1], method = "cubist",
                           preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_cubist_cv <- pxrf_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 9) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = pxrf_train_data$country)
pxrf_cubist_valid <- predict(pxrf_cubist_model, newdata = pxrf_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs= pxrf_valid_data$sand) %>%
                       add_column(country = pxrf_valid_data$country) %>%
                       rename(pred = value)
pxrf_cubist_plots <- validation_plot(pxrf_cubist_cv, pxrf_cubist_valid, "PXRF", "Cubist")
ggarrange(plotlist = pxrf_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pxrf_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
pxrf_cubist_importance <- varImp(pxrf_cubist_model)$importance %>%
                          rownames_to_column("variables") %>%
                          slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                          arrange(Overall) %>%
                          mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - Cubist Variable Importance")
ggsave("figures/sand/pxrf_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
                add_row(Dataset = "PXRF", Model = "PLS",
                        n_cv = nrow(pxrf_pls_cv), n_valid = nrow(pxrf_pls_valid),
                        RMSE_cv = RMSE(pxrf_pls_cv$pred, pxrf_pls_cv$obs),
                        R2_cv = caret::R2(pxrf_pls_cv$pred, pxrf_pls_cv$obs),
                        RMSE_valid = RMSE(pxrf_pls_valid$pred, pxrf_pls_valid$obs),
                        R2_valid = caret::R2(pxrf_pls_valid$pred, pxrf_pls_valid$obs)) %>%
                add_row(Dataset = "PXRF", Model = "RF",
                        n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
                        RMSE_cv = RMSE(pxrf_rf_cv$pred, pxrf_rf_cv$obs),
                        R2_cv = caret::R2(pxrf_rf_cv$pred, pxrf_rf_cv$obs),
                        RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
                        R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs)) %>%
                add_row(Dataset = "PXRF", Model = "Cubist",
                        n_cv = nrow(pxrf_cubist_cv), n_valid = nrow(pxrf_cubist_valid),
                        RMSE_cv = RMSE(pxrf_cubist_cv$pred, pxrf_cubist_cv$obs),
                        R2_cv = caret::R2(pxrf_cubist_cv$pred, pxrf_cubist_cv$obs),
                        RMSE_valid = RMSE(pxrf_cubist_valid$pred, pxrf_cubist_valid$obs),
                        R2_valid = caret::R2(pxrf_cubist_valid$pred, pxrf_cubist_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
pv_data_allcountries <- raw_sand_data %>%
                        select(country, sand, "K":"Pb", "350":"2500") %>%
                        # substituting numbers by legal variable names
                        rename_with(~ str_replace(., "^[0-9]+$", paste0("band_",. , "_nm"))) %>%
                        mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
                        drop_na()

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pv_data_allcountries$sand, p = 0.8, list = F)
pv_train_data <- pv_data_allcountries %>%
                 slice(partition_index)
pv_valid_data <- pv_data_allcountries %>%
                 slice(-partition_index)

## PLS
control_pls <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_pls <- c("nzv", "center", "scale")
set.seed(100)
pv_pls_model <- train(sand ~ ., data = pv_train_data[-1], method = "pls",
                      preProcess = preprocess_pls, trControl = control_pls)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pv_pls_cv <- pv_pls_model$pred %>%
             filter(ncomp == 3) %>% # best model with comps = 3
             arrange(rowIndex) %>%
             add_column(country = pv_train_data$country)
pv_pls_valid <- predict(pv_pls_model, newdata = pv_valid_data) %>%
                as_tibble() %>%
                add_column(obs = pv_valid_data$sand) %>%
                add_column(country = pv_valid_data$country) %>%
                rename(pred = value)
pv_pls_plots <- validation_plot(pv_pls_cv, pv_pls_valid, "PXRF + Vis-NIR", "PLS")
ggarrange(plotlist = pv_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pv_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PLS importance
pv_pls_importance <- varImp(pv_pls_model)$importance %>%
                     rownames_to_column("variables") %>%
                     slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                     arrange(Overall) %>%
                     mutate(variables = str_replace(variables, "^band_", "")) %>%
                     mutate(variables = str_replace(variables, "_nm$", "")) %>%
                     mutate(variables = factor(variables,levels = variables))
ggplot(pv_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - PLS Variable Importance")
ggsave("figures/sand/pv_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
control_rf <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_rf <- c("nzv", "center", "scale")
set.seed(100)
pv_rf_model <- train(sand ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess_rf, trControl = control_rf)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
            filter(mtry == 117) %>% # best model with mtry = 117
            arrange(rowIndex) %>%
            add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
               as_tibble() %>%
               add_column(obs = pv_valid_data$sand) %>%
               add_column(country = pv_valid_data$country) %>%
               rename(pred = value)
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid, "PXRF + Vis-NIR", "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
                    rownames_to_column("variables") %>%
                    slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                    arrange(Overall) %>%
                    mutate(variables = str_replace(variables, "^band_", "")) %>%
                    mutate(variables = str_replace(variables, "_nm$", "")) %>%
                    mutate(variables = factor(variables,levels = variables))
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/sand/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("nzv", "center", "scale")
set.seed(100)
pv_cubist_model <- train(sand ~ ., data = pv_train_data[-1], method = "cubist",
                         preProcess = preprocess_cubist, trControl = control_cubist)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pv_cubist_cv <- pv_cubist_model$pred %>%
                filter(committees == 20, neighbors == 0) %>% # filtering best model
                arrange(rowIndex) %>%
                add_column(country = pv_train_data$country)
pv_cubist_valid <- predict(pv_cubist_model, newdata = pv_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = pv_valid_data$sand) %>%
                   add_column(country = pv_valid_data$country) %>%
                   rename(pred = value)
pv_cubist_plots <- validation_plot(pv_cubist_cv, pv_cubist_valid, "PXRF + Vis-NIR", "Cubist")
ggarrange(plotlist = pv_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/sand/pv_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
pv_cubist_importance <- varImp(pv_cubist_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_replace(variables, "^band_", "")) %>%
                        mutate(variables = str_replace(variables, "_nm$", "")) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(pv_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - Cubist Variable Importance")
ggsave("figures/sand/pv_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
                add_row(Dataset = "PXRF + Vis-NIR", Model = "PLS",
                        n_cv = nrow(pv_pls_cv), n_valid = nrow(pv_pls_valid),
                        RMSE_cv = RMSE(pv_pls_cv$pred, pv_pls_cv$obs),
                        R2_cv = caret::R2(pv_pls_cv$pred, pv_pls_cv$obs),
                        RMSE_valid = RMSE(pv_pls_valid$pred, pv_pls_valid$obs),
                        R2_valid = caret::R2(pv_pls_valid$pred, pv_pls_valid$obs)) %>%
                add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
                        n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
                        RMSE_cv = RMSE(pv_rf_cv$pred, pv_rf_cv$obs),
                        R2_cv = caret::R2(pv_rf_cv$pred, pv_rf_cv$obs),
                        RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
                        R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs)) %>%
                add_row(Dataset = "PXRF + Vis-NIR", Model = "Cubist",
                        n_cv = nrow(pv_cubist_cv), n_valid = nrow(pv_cubist_valid),
                        RMSE_cv = RMSE(pv_cubist_cv$pred, pv_cubist_cv$obs),
                        R2_cv = caret::R2(pv_cubist_cv$pred, pv_cubist_cv$obs),
                        RMSE_valid = RMSE(pv_cubist_valid$pred, pv_cubist_valid$obs),
                        R2_valid = caret::R2(pv_cubist_valid$pred, pv_cubist_valid$obs))
write_excel_csv(model_scores, "tables/sand/model_scores.csv")
