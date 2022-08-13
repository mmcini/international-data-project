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
               select(country, OC, "350":"2500") %>%
               drop_na() %>%
               # substituting numbers by legal variable names
               rename_with(~ str_replace(., "^[0-9]+$", paste0("band_", ., "_nm")))

## PXRF
pxrf_data <- raw_data %>%
             select(country, OC, "K":"Pb", "350":"2500") %>%
             mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
             drop_na() %>%
             select(country, OC, "K":"Pb")

## PXRF + Vis-NIR
pv_data <- raw_data %>%
           select(country, OC, "K":"Pb", "350":"2500") %>%
           # substituting numbers by legal variable names
           rename_with(~ str_replace(., "^[0-9]+$", paste0("band_", ., "_nm"))) %>%
           mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
           drop_na()

# All countries ###################################################################################

## Vis-NIR models ##################################################################################
## Partitioning
set.seed(100)
partition_index <- createDataPartition(visnir_data$OC, p = 0.8, list = F)
visnir_train_data <- visnir_data %>%
                     slice(partition_index)
visnir_valid_data <- visnir_data %>%
                     slice(-partition_index)

## Models
## PLS
set.seed(100)
visnir_pls_model <- train(OC ~ ., data = visnir_train_data[-1], method = "pls",
                          preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
visnir_pls_cv <- visnir_pls_model$pred %>%
                 filter(ncomp == 3) %>% # best model with comps = 3
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_pls_valid <- predict(visnir_pls_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$OC) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_pls_plots <- validation_plot(visnir_pls_cv, visnir_pls_valid,
                                    variable = "OC", dataset =  "Vis-NIR", model = "PLS")
ggarrange(plotlist = visnir_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/visnir_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/visnir_pls_coefs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/visnir_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
visnir_rf_model <- train(OC ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
                filter(mtry == 2) %>% # best model
                arrange(rowIndex) %>%
                add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = visnir_valid_data$OC) %>%
                   add_column(country = visnir_valid_data$country) %>%
                   rename(pred = value)
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "OC", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
visnir_cubist_model <- train(OC ~ ., data = visnir_train_data[-1], method = "cubist",
                             preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
visnir_cubist_cv <- visnir_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 0) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = visnir_train_data$country)
visnir_cubist_valid <- predict(visnir_cubist_model, newdata = visnir_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = visnir_valid_data$OC) %>%
                       add_column(country = visnir_valid_data$country) %>%
                       rename(pred = value)
visnir_cubist_plots <- validation_plot(visnir_cubist_cv, visnir_cubist_valid,
                                       variable = "OC", dataset = "Vis-NIR", model = "Cubist")
ggarrange(plotlist = visnir_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/visnir_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/visnir_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
visnir_xgb_model <- train(OC ~ ., data = visnir_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
visnir_xgb_cv <- visnir_xgb_model$pred %>%
                 filter(nrounds == 50, max_depth == 3,
                        eta == 0.3,gamma == 0,colsample_bytree == 0.6,
                        min_child_weight == 1, subsample == 0.75) %>% # filtering best model
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_xgb_valid <- predict(visnir_xgb_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$OC) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_xgb_plots <- validation_plot(visnir_xgb_cv, visnir_xgb_valid,
                                    variable = "OC", dataset = "Vis-NIR", model = "XGB")
ggarrange(plotlist = visnir_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/visnir_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
visnir_xgb_importance <- varImp(visnir_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - XGB Variable Importance")
ggsave("figures/OC/allcountries/visnir_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(visnir_cubist_valid$pred, visnir_cubist_valid$obs)) %>%
                add_row(Dataset = "Vis-NIR", Model = "XGB",
                        n_cv = nrow(visnir_xgb_cv), n_valid = nrow(visnir_xgb_valid),
                        RMSE_cv = RMSE(visnir_xgb_cv$pred, visnir_xgb_cv$obs),
                        R2_cv = caret::R2(visnir_xgb_cv$pred, visnir_xgb_cv$obs),
                        RMSE_valid = RMSE(visnir_xgb_valid$pred, visnir_xgb_valid$obs),
                        R2_valid = caret::R2(visnir_xgb_valid$pred, visnir_xgb_valid$obs))

## PXRF models #####################################################################################
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data[-c(1, 2)], y = pxrf_data[["OC"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
### [1] "Ca" "Sr" "K"  "Cu" "Zr" "Mn" "Fe" "Cr" "Rb" "Ti" "V"  "Zn" "Ni" "Pb"

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data$OC, p = 0.8, list = F)
pxrf_train_data <- pxrf_data %>%
                   select(country, OC, all_of(pxrf_selected_predictors)) %>%
                   slice(partition_index)
pxrf_valid_data <- pxrf_data %>%
                   select(country, OC, all_of(pxrf_selected_predictors)) %>%
                   slice(-partition_index)

## PLS
set.seed(100)
pxrf_pls_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "pls",
                        preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_pls_cv <- pxrf_pls_model$pred %>%
               filter(ncomp == 3) %>% # best model with comps = 3
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_pls_valid <- predict(pxrf_pls_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$OC) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_pls_plots <- validation_plot(pxrf_pls_cv, pxrf_pls_valid,
                                  variable = "OC", dataset = "PXRF", model = "PLS")
ggarrange(plotlist = pxrf_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pxrf_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PLS importance
pxrf_pls_importance <- varImp(pxrf_pls_model)$importance %>%
                       rownames_to_column("variables") %>%
                       slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                       arrange(Overall) %>%
                       mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - PLS Variable Importance")
ggsave("figures/OC/allcountries/pxrf_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pxrf_rf_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
              filter(mtry == 2) %>% # best model
              arrange(rowIndex) %>%
              add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
                 as_tibble() %>%
                 add_column(obs = pxrf_valid_data$OC) %>%
                 add_column(country = pxrf_valid_data$country) %>%
                 rename(pred = value)
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "OC", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
                      rownames_to_column("variables") %>%
                      slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                      arrange(Overall) %>%
                      mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - RF Variable Importance")
ggsave("figures/OC/allcountries/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
control_cubist <- trainControl(method = "cv", number = 10, savePredictions = T)
preprocess_cubist <- c("nzv", "center", "scale")
set.seed(100)
pxrf_cubist_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "cubist",
                           preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_cubist_cv <- pxrf_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 0) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = pxrf_train_data$country)
pxrf_cubist_valid <- predict(pxrf_cubist_model, newdata = pxrf_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = pxrf_valid_data$OC) %>%
                       add_column(country = pxrf_valid_data$country) %>%
                       rename(pred = value)
pxrf_cubist_plots <- validation_plot(pxrf_cubist_cv, pxrf_cubist_valid,
                                     variable = "OC", dataset = "PXRF", model = "Cubist")
ggarrange(plotlist = pxrf_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pxrf_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
pxrf_cubist_importance <- varImp(pxrf_cubist_model)$importance %>%
                          rownames_to_column("variables") %>%
                          slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                          arrange(Overall) %>%
                          mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - Cubist Variable Importance")
ggsave("figures/OC/allcountries/pxrf_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
pxrf_xgb_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_xgb_cv <- pxrf_xgb_model$pred %>%
               filter(nrounds == 50, max_depth == 3,
                      eta == 0.3,gamma == 0,colsample_bytree == 0.6,
                      min_child_weight == 1, subsample == 1) %>% # filtering best model
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_xgb_valid <- predict(pxrf_xgb_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$OC) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_xgb_plots <- validation_plot(pxrf_xgb_cv, pxrf_xgb_valid,
                                    variable = "OC", dataset = "PXRF", model = "XGB")
ggarrange(plotlist = pxrf_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pxrf_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
pxrf_xgb_importance <- varImp(pxrf_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - XGB Variable Importance")
ggsave("figures/OC/allcountries/pxrf_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(pxrf_cubist_valid$pred, pxrf_cubist_valid$obs)) %>%
                add_row(Dataset = "PXRF", Model = "XGB",
                        n_cv = nrow(pxrf_xgb_cv), n_valid = nrow(pxrf_xgb_valid),
                        RMSE_cv = RMSE(pxrf_xgb_cv$pred, pxrf_xgb_cv$obs),
                        R2_cv = caret::R2(pxrf_xgb_cv$pred, pxrf_xgb_cv$obs),
                        RMSE_valid = RMSE(pxrf_xgb_valid$pred, pxrf_xgb_valid$obs),
                        R2_valid = caret::R2(pxrf_xgb_valid$pred, pxrf_xgb_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
set.seed(100)
partition_index <- createDataPartition(pv_data$OC, p = 0.8, list = F)
pv_train_data <- pv_data %>%
                 select(country, OC, all_of(pxrf_selected_predictors),
                      "band_350_nm":"band_2500_nm") %>%
                 slice(partition_index)
pv_valid_data <- pv_data %>%
                 select(country, OC, all_of(pxrf_selected_predictors),
                      "band_350_nm":"band_2500_nm") %>%
                 slice(-partition_index)

## PLS
set.seed(100)
pv_pls_model <- train(OC ~ ., data = pv_train_data[-1], method = "pls",
                      preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pv_pls_cv <- pv_pls_model$pred %>%
             filter(ncomp == 3) %>% # best model with comps = 3
             arrange(rowIndex) %>%
             add_column(country = pv_train_data$country)
pv_pls_valid <- predict(pv_pls_model, newdata = pv_valid_data) %>%
                as_tibble() %>%
                add_column(obs = pv_valid_data$OC) %>%
                add_column(country = pv_valid_data$country) %>%
                rename(pred = value)
pv_pls_plots <- validation_plot(pv_pls_cv, pv_pls_valid,
                                variable = "OC", dataset = "PXRF + Vis-NIR", model = "PLS")
ggarrange(plotlist = pv_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pv_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/pv_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pv_rf_model <- train(OC ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
            filter(mtry == 116) %>% # best model
            arrange(rowIndex) %>%
            add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
               as_tibble() %>%
               add_column(obs = pv_valid_data$OC) %>%
               add_column(country = pv_valid_data$country) %>%
               rename(pred = value)
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "OC", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pv_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
pv_cubist_model <- train(OC ~ ., data = pv_train_data[-1], method = "cubist",
                         preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pv_cubist_cv <- pv_cubist_model$pred %>%
                filter(committees == 20, neighbors == 0) %>% # filtering best model
                arrange(rowIndex) %>%
                add_column(country = pv_train_data$country)
pv_cubist_valid <- predict(pv_cubist_model, newdata = pv_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = pv_valid_data$OC) %>%
                   add_column(country = pv_valid_data$country) %>%
                   rename(pred = value)
pv_cubist_plots <- validation_plot(pv_cubist_cv, pv_cubist_valid,
                                   variable = "OC", dataset = "PXRF + Vis-NIR", model = "Cubist")
ggarrange(plotlist = pv_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pv_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/allcountries/pv_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
pv_xgb_model <- train(OC ~ ., data = pv_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
pv_xgb_cv <- pv_xgb_model$pred %>%
               filter(nrounds == 50, max_depth == 2,
                      eta == 0.3,gamma == 0,colsample_bytree == 0.6,
                      min_child_weight == 1, subsample == 1) %>% # filtering best model
               arrange(rowIndex) %>%
               add_column(country = pv_train_data$country)
pv_xgb_valid <- predict(pv_xgb_model, newdata = pv_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pv_valid_data$OC) %>%
                  add_column(country = pv_valid_data$country) %>%
                  rename(pred = value)
pv_xgb_plots <- validation_plot(pv_xgb_cv, pv_xgb_valid,
                                    variable = "OC", dataset = "PXRF + Vis-NIR", model = "XGB")
ggarrange(plotlist = pv_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/allcountries/pv_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
pv_xgb_importance <- varImp(pv_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_replace(variables, "^band_", "")) %>%
                        mutate(variables = str_replace(variables, "_nm$", "")) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(pv_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - XGB Variable Importance")
ggsave("figures/OC/allcountries/pv_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(pv_cubist_valid$pred, pv_cubist_valid$obs)) %>%
                add_row(Dataset = "PXRF + Vis-NIR", Model = "XGB",
                        n_cv = nrow(pv_xgb_cv), n_valid = nrow(pv_xgb_valid),
                        RMSE_cv = RMSE(pv_xgb_cv$pred, pv_xgb_cv$obs),
                        R2_cv = caret::R2(pv_xgb_cv$pred, pv_xgb_cv$obs),
                        RMSE_valid = RMSE(pv_xgb_valid$pred, pv_xgb_valid$obs),
                        R2_valid = caret::R2(pv_xgb_valid$pred, pv_xgb_valid$obs))
write_excel_csv(model_scores, "tables/OC/allcountries/model_scores.csv")

# US ##############################################################################################

## Vis-NIR models ##################################################################################
## Partitioning
visnir_data_us <- visnir_data %>%
                  filter(country == "US")
set.seed(100)
partition_index <- createDataPartition(visnir_data$OC, p = 0.8, list = F)
visnir_train_data <- visnir_data_us %>%
                     slice(partition_index)
visnir_valid_data <- visnir_data_us %>%
                     slice(-partition_index)

## Models
## PLS
set.seed(100)
visnir_pls_model <- train(OC ~ ., data = visnir_train_data[-1], method = "pls",
                          preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
visnir_pls_cv <- visnir_pls_model$pred %>%
                 filter(ncomp == 3) %>% # best model with comps = 3
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_pls_valid <- predict(visnir_pls_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$OC) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_pls_plots <- validation_plot(visnir_pls_cv, visnir_pls_valid,
                                    variable = "OC", dataset =  "Vis-NIR", model = "PLS")
ggarrange(plotlist = visnir_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/visnir_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/visnir_pls_coefs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/visnir_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
visnir_rf_model <- train(OC ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
                filter(mtry == 109) %>% # best model
                arrange(rowIndex) %>%
                add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = visnir_valid_data$OC) %>%
                   add_column(country = visnir_valid_data$country) %>%
                   rename(pred = value)
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "OC", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
visnir_cubist_model <- train(OC ~ ., data = visnir_train_data[-1], method = "cubist",
                             preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
visnir_cubist_cv <- visnir_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 0) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = visnir_train_data$country)
visnir_cubist_valid <- predict(visnir_cubist_model, newdata = visnir_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = visnir_valid_data$OC) %>%
                       add_column(country = visnir_valid_data$country) %>%
                       rename(pred = value)
visnir_cubist_plots <- validation_plot(visnir_cubist_cv, visnir_cubist_valid,
                                       variable = "OC", dataset = "Vis-NIR", model = "Cubist")
ggarrange(plotlist = visnir_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/visnir_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/visnir_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
visnir_xgb_model <- train(OC ~ ., data = visnir_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
visnir_xgb_cv <- visnir_xgb_model$pred %>%
                 filter(nrounds == 150, max_depth == 1,
                        eta == 0.4,gamma == 0,colsample_bytree == 0.8,
                        min_child_weight == 1, subsample == 0.5) %>% # filtering best model
                 arrange(rowIndex) %>%
                 add_column(country = visnir_train_data$country)
visnir_xgb_valid <- predict(visnir_xgb_model, newdata = visnir_valid_data) %>%
                    as_tibble() %>%
                    add_column(obs = visnir_valid_data$OC) %>%
                    add_column(country = visnir_valid_data$country) %>%
                    rename(pred = value)
visnir_xgb_plots <- validation_plot(visnir_xgb_cv, visnir_xgb_valid,
                                    variable = "OC", dataset = "Vis-NIR", model = "XGB")
ggarrange(plotlist = visnir_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/visnir_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
visnir_xgb_importance <- varImp(visnir_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(visnir_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("Vis-NIR - XGB Variable Importance")
ggsave("figures/OC/US/visnir_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(visnir_cubist_valid$pred, visnir_cubist_valid$obs)) %>%
                add_row(Dataset = "Vis-NIR", Model = "XGB",
                        n_cv = nrow(visnir_xgb_cv), n_valid = nrow(visnir_xgb_valid),
                        RMSE_cv = RMSE(visnir_xgb_cv$pred, visnir_xgb_cv$obs),
                        R2_cv = caret::R2(visnir_xgb_cv$pred, visnir_xgb_cv$obs),
                        RMSE_valid = RMSE(visnir_xgb_valid$pred, visnir_xgb_valid$obs),
                        R2_valid = caret::R2(visnir_xgb_valid$pred, visnir_xgb_valid$obs))

## PXRF models #####################################################################################
pxrf_data_us <- pxrf_data %>%
               filter(country == "US")
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_us[-c(1, 2)], y = pxrf_data_us[["OC"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
### [1] "Cu" "Sr" "Mn" "Pb" "Cr" "Fe" "Rb" "As" "Ca" "Ti" "Zr"

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_us$OC, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_us %>%
                   select(country, OC, all_of(pxrf_selected_predictors)) %>%
                   slice(partition_index)
pxrf_valid_data <- pxrf_data_us %>%
                   select(country, OC, all_of(pxrf_selected_predictors)) %>%
                   slice(-partition_index)

## PLS
set.seed(100)
pxrf_pls_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "pls",
                        preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_pls_cv <- pxrf_pls_model$pred %>%
               filter(ncomp == 3) %>% # best model with comps = 3
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_pls_valid <- predict(pxrf_pls_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$OC) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_pls_plots <- validation_plot(pxrf_pls_cv, pxrf_pls_valid,
                                  variable = "OC", dataset = "PXRF", model = "PLS")
ggarrange(plotlist = pxrf_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pxrf_pls_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PLS importance
pxrf_pls_importance <- varImp(pxrf_pls_model)$importance %>%
                       rownames_to_column("variables") %>%
                       slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                       arrange(Overall) %>%
                       mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_pls_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - PLS Variable Importance")
ggsave("figures/OC/US/pxrf_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pxrf_rf_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
              filter(mtry == 2) %>% # best model
              arrange(rowIndex) %>%
              add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
                 as_tibble() %>%
                 add_column(obs = pxrf_valid_data$OC) %>%
                 add_column(country = pxrf_valid_data$country) %>%
                 rename(pred = value)
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "OC", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
                      rownames_to_column("variables") %>%
                      slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                      arrange(Overall) %>%
                      mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - RF Variable Importance")
ggsave("figures/OC/US/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
pxrf_cubist_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "cubist",
                           preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_cubist_cv <- pxrf_cubist_model$pred %>%
                    filter(committees == 20, neighbors == 5) %>% # filtering best model
                    arrange(rowIndex) %>%
                    add_column(country = pxrf_train_data$country)
pxrf_cubist_valid <- predict(pxrf_cubist_model, newdata = pxrf_valid_data) %>%
                       as_tibble() %>%
                       add_column(obs = pxrf_valid_data$OC) %>%
                       add_column(country = pxrf_valid_data$country) %>%
                       rename(pred = value)
pxrf_cubist_plots <- validation_plot(pxrf_cubist_cv, pxrf_cubist_valid,
                                     variable = "OC", dataset = "PXRF", model = "Cubist")
ggarrange(plotlist = pxrf_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pxrf_cubist_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## Cubist importance
pxrf_cubist_importance <- varImp(pxrf_cubist_model)$importance %>%
                          rownames_to_column("variables") %>%
                          slice_max(n = 10, order_by = Overall) %>% # 10 biggest values
                          arrange(Overall) %>%
                          mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_cubist_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - Cubist Variable Importance")
ggsave("figures/OC/US/pxrf_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
pxrf_xgb_model <- train(OC ~ ., data = pxrf_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_xgb_cv <- pxrf_xgb_model$pred %>%
               filter(nrounds == 50, max_depth == 3,
                      eta == 0.3,gamma == 0,colsample_bytree == 0.6,
                      min_child_weight == 1, subsample == 1) %>% # filtering best model
               arrange(rowIndex) %>%
               add_column(country = pxrf_train_data$country)
pxrf_xgb_valid <- predict(pxrf_xgb_model, newdata = pxrf_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pxrf_valid_data$OC) %>%
                  add_column(country = pxrf_valid_data$country) %>%
                  rename(pred = value)
pxrf_xgb_plots <- validation_plot(pxrf_xgb_cv, pxrf_xgb_valid,
                                    variable = "OC", dataset = "PXRF", model = "XGB")
ggarrange(plotlist = pxrf_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pxrf_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
pxrf_xgb_importance <- varImp(pxrf_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
                        mutate(variables = factor(variables,levels = variables))
ggplot(pxrf_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF - XGB Variable Importance")
ggsave("figures/OC/US/pxrf_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(pxrf_cubist_valid$pred, pxrf_cubist_valid$obs)) %>%
                add_row(Dataset = "PXRF", Model = "XGB",
                        n_cv = nrow(pxrf_xgb_cv), n_valid = nrow(pxrf_xgb_valid),
                        RMSE_cv = RMSE(pxrf_xgb_cv$pred, pxrf_xgb_cv$obs),
                        R2_cv = caret::R2(pxrf_xgb_cv$pred, pxrf_xgb_cv$obs),
                        RMSE_valid = RMSE(pxrf_xgb_valid$pred, pxrf_xgb_valid$obs),
                        R2_valid = caret::R2(pxrf_xgb_valid$pred, pxrf_xgb_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
pv_data_us <- pv_data %>%
              filter(country == "US")
set.seed(100)
partition_index <- createDataPartition(pv_data_us$OC, p = 0.8, list = F)
pv_train_data <- pv_data_us %>%
                 select(country, OC, all_of(pxrf_selected_predictors),
                      "band_350_nm":"band_2500_nm") %>%
                 slice(partition_index)
pv_valid_data <- pv_data_us %>%
                 select(country, OC, all_of(pxrf_selected_predictors),
                      "band_350_nm":"band_2500_nm") %>%
                 slice(-partition_index)

## PLS
set.seed(100)
pv_pls_model <- train(OC ~ ., data = pv_train_data[-1], method = "pls",
                      preProcess = preprocess, trControl = control)

## PLS results - kfold cross validation (80%) and hold-out validation (20%)
pv_pls_cv <- pv_pls_model$pred %>%
             filter(ncomp == 3) %>% # best model with comps = 3
             arrange(rowIndex) %>%
             add_column(country = pv_train_data$country)
pv_pls_valid <- predict(pv_pls_model, newdata = pv_valid_data) %>%
                as_tibble() %>%
                add_column(obs = pv_valid_data$OC) %>%
                add_column(country = pv_valid_data$country) %>%
                rename(pred = value)
pv_pls_plots <- validation_plot(pv_pls_cv, pv_pls_valid,
                                variable = "OC", dataset = "PXRF + Vis-NIR", model = "PLS")
ggarrange(plotlist = pv_pls_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pv_pls_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/pv_pls_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## RF
set.seed(100)
pv_rf_model <- train(OC ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
            filter(mtry == 114) %>% # best model
            arrange(rowIndex) %>%
            add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
               as_tibble() %>%
               add_column(obs = pv_valid_data$OC) %>%
               add_column(country = pv_valid_data$country) %>%
               rename(pred = value)
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "OC", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pv_rf_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Cubist
set.seed(100)
pv_cubist_model <- train(OC ~ ., data = pv_train_data[-1], method = "cubist",
                         preProcess = preprocess, trControl = control)

## Cubist results - kfold cross validation (80%) and hold-out validation (20%)
pv_cubist_cv <- pv_cubist_model$pred %>%
                filter(committees == 20, neighbors == 0) %>% # filtering best model
                arrange(rowIndex) %>%
                add_column(country = pv_train_data$country)
pv_cubist_valid <- predict(pv_cubist_model, newdata = pv_valid_data) %>%
                   as_tibble() %>%
                   add_column(obs = pv_valid_data$OC) %>%
                   add_column(country = pv_valid_data$country) %>%
                   rename(pred = value)
pv_cubist_plots <- validation_plot(pv_cubist_cv, pv_cubist_valid,
                                   variable = "OC", dataset = "PXRF + Vis-NIR", model = "Cubist")
ggarrange(plotlist = pv_cubist_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pv_cubist_pred_obs.png", dpi = 300, units = "mm",
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
ggsave("figures/OC/US/pv_cubist_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## XGB
set.seed(100)
pv_xgb_model <- train(OC ~ ., data = pv_train_data[-1], method = "xgbTree",
                             preProcess = preprocess, trControl = control)

## XGB results - kfold cross validation (80%) and hold-out validation (20%)
pv_xgb_cv <- pv_xgb_model$pred %>%
               filter(nrounds == 50, max_depth == 2,
                      eta == 0.3,gamma == 0,colsample_bytree == 0.8,
                      min_child_weight == 1, subsample == 0.5) %>% # filtering best model
               arrange(rowIndex) %>%
               add_column(country = pv_train_data$country)
pv_xgb_valid <- predict(pv_xgb_model, newdata = pv_valid_data) %>%
                  as_tibble() %>%
                  add_column(obs = pv_valid_data$OC) %>%
                  add_column(country = pv_valid_data$country) %>%
                  rename(pred = value)
pv_xgb_plots <- validation_plot(pv_xgb_cv, pv_xgb_valid,
                                    variable = "OC", dataset = "PXRF + Vis-NIR", model = "XGB")
ggarrange(plotlist = pv_xgb_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/OC/US/pv_xgb_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## XGB importance
pv_xgb_importance <- varImp(pv_xgb_model)$importance %>%
                        rownames_to_column("variables") %>%
                        slice_max(n = 30, order_by = Overall) %>% # 30 biggest values
                        arrange(Overall) %>%
                        mutate(variables = str_replace(variables, "^band_", "")) %>%
                        mutate(variables = str_replace(variables, "_nm$", "")) %>%
                        mutate(variables = factor(variables,levels = variables))
ggplot(pv_xgb_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
       ggtitle("PXRF + Vis-NIR - XGB Variable Importance")
ggsave("figures/OC/US/pv_xgb_importance.png", dpi = 300, units = "mm",
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
                        R2_valid = caret::R2(pv_cubist_valid$pred, pv_cubist_valid$obs)) %>%
                add_row(Dataset = "PXRF + Vis-NIR", Model = "XGB",
                        n_cv = nrow(pv_xgb_cv), n_valid = nrow(pv_xgb_valid),
                        RMSE_cv = RMSE(pv_xgb_cv$pred, pv_xgb_cv$obs),
                        R2_cv = caret::R2(pv_xgb_cv$pred, pv_xgb_cv$obs),
                        RMSE_valid = RMSE(pv_xgb_valid$pred, pv_xgb_valid$obs),
                        R2_valid = caret::R2(pv_xgb_valid$pred, pv_xgb_valid$obs))
write_excel_csv(model_scores, "tables/OC/US/model_scores.csv")
