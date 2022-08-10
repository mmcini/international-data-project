# Libs and functions ###############################################################################

source("R/_functions.R")

# Modeling　########################################################################################

raw_data <- read_excel("data/oc_texture_data.xlsx", na = "NA")

## Control and preprocessing parameters
control <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess <- c("nzv", "center", "scale")

## Datasets
## Vis-NIR
visnir_data <- raw_data %>%
               select(country, clay, "350":"2500") %>%
               drop_na() %>%
               # substituting numbers by legal variable names
               rename_with(~ str_replace(., "^[0-9]+$", paste0("band_", ., "_nm")))

## PXRF
pxrf_data <- raw_data %>%
             select(country, clay, "K":"Pb") %>%
             mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
             drop_na()

## PXRF + Vis-NIR
pv_data <- raw_data %>%
           select(country, clay, "K":"Pb", "350":"2500") %>%
           # substituting numbers by legal variable names
           rename_with(~ str_replace(., "^[0-9]+$", paste0("band_", ., "_nm"))) %>%
           mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
           drop_na()

# All countries ###################################################################################

## Vis-NIR models ##################################################################################
## Partitioning
set.seed(100)
partition_index <- createDataPartition(visnir_data$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data %>%
                     slice(partition_index)
visnir_valid_data <- visnir_data %>%
                     slice(-partition_index)

## Models
## PLS
set.seed(100)
visnir_pls_model <- train(clay ~ ., data = visnir_train_data[-1], method = "pls",
                          preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
visnir_pls_cv <- visnir_pls_model$pred %>%
                 filter(ncomp == 3) %>% # best model with comps = 3
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_pls_valid <- predict(visnir_pls_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$clay) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_pls_plots <- validation_plot(visnir_pls_cv, visnir_pls_valid,
                                    variable = "clay", dataset =  "Vis-NIR", model = "PLS")
ggarrange(plotlist = visnir_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/visnir_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/visnir_pls_coefs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/visnir_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
                filter(mtry == 216) %>% # best model with mtry = 109
                arrange(rowIndex) %>%
                add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = visnir_valid_data$clay) %>%
                   add_column(country = visnir_valid_data$country) %>%
                   rename(pred = value)
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
visnir_cubist_model <- train(clay ~ ., data = visnir_train_data[-1], method = "cubist",
                             preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
visnir_cubist_cv <- visnir_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 9) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = visnir_train_data$country)
visnir_cubist_valid <- predict(visnir_cubist_model, newdata = visnir_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = visnir_valid_data$clay) %>%
                       add_column(country = visnir_valid_data$country) %>%
                       rename(pred = value)
visnir_cubist_plots <- validation_plot(visnir_cubist_cv, visnir_cubist_valid,
                                       variable = "clay", dataset = "Vis-NIR", model = "Cubist")
ggarrange(plotlist = visnir_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/visnir_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/visnir_cubist_importance.png", dpi = 300, units = "mm",
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
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data[-c(1, 2)], y = pxrf_data[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
## [1] "Fe" "Ti" "Sr" "Zn" "K"  "Cu" "Cr" "Mn" "Rb" "Zr" "Ba" "Ca" "As" "Pb" "Ni"

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data %>%
                   select(country, clay, all_of(pxrf_selected_predictors)) %>%
                   slice(partition_index)
pxrf_valid_data <- pxrf_data %>%
                   select(country, clay, all_of(pxrf_selected_predictors)) %>%
                   slice(-partition_index)

## PLS
set.seed(100)
pxrf_pls_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "pls",
                        preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_pls_cv <- pxrf_pls_model$pred %>%
               filter(ncomp == 3) %>% # best model with comps = 3
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_pls_valid <- predict(pxrf_pls_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$clay) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_pls_plots <- validation_plot(pxrf_pls_cv, pxrf_pls_valid,
                                  variable = "clay", dataset = "PXRF", model = "PLS")
ggarrange(plotlist = pxrf_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pxrf_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PLS importance
pxrf_pls_importance <- varImp(pxrf_pls_model)$importance %>%
                       rownames_to_column("variables") %>%
                       slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                       arrange(Overall) %>%
                       mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - PLS Variable Importance")
ggsave("figures/clay/allcountries/pxrf_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
              filter(mtry == 8) %>% # best model
              arrange(rowIndex) %>%
              add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
                 as_tibble() %>%
                 add_column(obs = pxrf_valid_data$clay) %>%
                 add_column(country = pxrf_valid_data$country) %>%
                 rename(pred = value)
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
                      rownames_to_column("variables") %>%
                      slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                      arrange(Overall) %>%
                      mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/allcountries/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("nzv", "center", "scale")
set.seed(100)
pxrf_cubist_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "cubist",
                           preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_cubist_cv <- pxrf_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 9) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = pxrf_train_data$country)
pxrf_cubist_valid <- predict(pxrf_cubist_model, newdata = pxrf_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = pxrf_valid_data$clay) %>%
                       add_column(country = pxrf_valid_data$country) %>%
                       rename(pred = value)
pxrf_cubist_plots <- validation_plot(pxrf_cubist_cv, pxrf_cubist_valid,
                                     variable = "clay", dataset = "PXRF", model = "Cubist")
ggarrange(plotlist = pxrf_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pxrf_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
pxrf_cubist_importance <- varImp(pxrf_cubist_model)$importance %>%
                          rownames_to_column("variables") %>%
                          slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                          arrange(Overall) %>%
                          mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - Cubist Variable Importance")
ggsave("figures/clay/allcountries/pxrf_cubist_importance.png", dpi = 300, units = "mm",
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
## Partitioning
set.seed(100)
partition_index <- createDataPartition(pv_data$clay, p = 0.8, list = F)
pv_train_data <- pv_data %>%
                 select(country, clay, all_of(pxrf_selected_predictors),
                        "band_350_nm":"band_2500_nm") %>%
                 slice(partition_index)
pv_valid_data <- pv_data %>%
                 select(country, clay, all_of(pxrf_selected_predictors),
                        "band_350_nm":"band_2500_nm") %>%
                 slice(-partition_index)

## PLS
set.seed(100)
pv_pls_model <- train(clay ~ ., data = pv_train_data[-1], method = "pls",
                      preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pv_pls_cv <- pv_pls_model$pred %>%
             filter(ncomp == 3) %>% # best model with comps = 3
             arrange(rowIndex) %>%
             add_column(country = pv_train_data$country)
pv_pls_valid <- predict(pv_pls_model, newdata = pv_valid_data) %>%
                as_tibble() %>%
                add_column(obs = pv_valid_data$clay) %>%
                add_column(country = pv_valid_data$country) %>%
                rename(pred = value)
pv_pls_plots <- validation_plot(pv_pls_cv, pv_pls_valid,
                                variable = "clay", dataset = "PXRF + Vis-NIR", model = "PLS")
ggarrange(plotlist = pv_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pv_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/pv_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
            filter(mtry == 231) %>% # best model
            arrange(rowIndex) %>%
            add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
               as_tibble() %>%
               add_column(obs = pv_valid_data$clay) %>%
               add_column(country = pv_valid_data$country) %>%
               rename(pred = value)
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pv_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
pv_cubist_model <- train(clay ~ ., data = pv_train_data[-1], method = "cubist",
                         preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pv_cubist_cv <- pv_cubist_model$pred %>%
                filter(committees == 20, neighbors == 0) %>% # filtering best model
                arrange(rowIndex) %>%
                add_column(country = pv_train_data$country)
pv_cubist_valid <- predict(pv_cubist_model, newdata = pv_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = pv_valid_data$clay) %>%
                   add_column(country = pv_valid_data$country) %>%
                   rename(pred = value)
pv_cubist_plots <- validation_plot(pv_cubist_cv, pv_cubist_valid,
                                   variable = "clay", dataset = "PXRF + Vis-NIR", model = "Cubist")
ggarrange(plotlist = pv_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pv_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/clay/allcountries/pv_cubist_importance.png", dpi = 300, units = "mm",
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
write_excel_csv(model_scores, "tables/clay/allcountries/model_scores.csv")
