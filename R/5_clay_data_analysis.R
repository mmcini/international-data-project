# Libs and functions ###############################################################################

source("R/_functions.R")

# Modeling　########################################################################################

raw_data <- read_excel("data/oc_texture_data.xlsx", na = "NA") %>%
   mutate(country = factor(country))

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
   select(country, clay, "K":"Pb", "350":"2500") %>%
   mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
   drop_na() %>%
   select(country, clay, "K":"Pb")

## PXRF + Vis-NIR
pv_data <- raw_data %>%
   select(country, clay, "K":"Pb", "350":"2500") %>%
   # substituting numbers by legal variable names
   rename_with(~ str_replace(., "^[0-9]+$", paste0("band_", ., "_nm"))) %>%
   mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
   drop_na()

# All countries ####################################################################################

## Vis-NIR models ##################################################################################
## Partitioning
set.seed(100)
partition_index <- createDataPartition(visnir_data$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data %>%
   slice(partition_index)
visnir_valid_data <- visnir_data %>%
   slice(-partition_index)
mtry_param <- round((ncol(visnir_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## Models
## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
   as_tibble() %>%
   add_column(obs = visnir_valid_data$clay) %>%
   add_column(country = visnir_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(visnir_rf_valid, "tables/clay/allcountries/visnir_rf_valid.csv")
write_excel_csv(visnir_rf_cv, "tables/clay/allcountries/visnir_rf_cv.csv")
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR",
                                   model = "RF", group_by = "country")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
   mutate(variables = factor(variables,levels = variables))
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("Vis-NIR - RF Variable Importance")
write_excel_csv(visnir_rf_importance, "tables/clay/allcountries/visnir_rf_importance.csv")
ggsave("figures/clay/allcountries/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "Vis-NIR", Model = "RF",
           n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
           RMSE_cv = mean_from_folds(visnir_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(visnir_rf_cv, type = "R2"),
           RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
           R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs))

## PXRF models #####################################################################################
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data[-c(1, 2)], y = pxrf_data[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_allcountries_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   as_tibble() %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/allcountries/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/allcountries/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF",
                                 model = "RF", group_by = "country")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/allcountries/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/allcountries/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))

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
mtry_param <- round((ncol(pv_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pv_valid_data$clay) %>%
   add_column(country = pv_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pv_rf_valid, "tables/clay/allcountries/pv_rf_valid.csv")
write_excel_csv(pv_rf_cv, "tables/clay/allcountries/pv_rf_cv.csv")
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR",
                               model = "RF", group_by = "country")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_replace(variables, "^band_", "")) %>%
   mutate(variables = str_replace(variables, "_nm$", "")) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pv_rf_importance, "tables/clay/allcountries/pv_rf_importance.csv")
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/clay/allcountries/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
   add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
           n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
           RMSE_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
           R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/allcountries/model_scores.csv")

# Brazil ###########################################################################################
## Vis-NIR models ##################################################################################
## Partitioning
visnir_data_br <- visnir_data %>%
   filter(country == "Brazil")
set.seed(100)
partition_index <- createDataPartition(visnir_data_br$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data_br %>%
   slice(partition_index)
visnir_valid_data <- visnir_data_br %>%
   slice(-partition_index)
mtry_param <- round((ncol(visnir_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## Models
## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
   as_tibble() %>%
   add_column(obs = visnir_valid_data$clay) %>%
   add_column(country = visnir_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(visnir_rf_valid, "tables/clay/Brazil/visnir_rf_valid.csv")
write_excel_csv(visnir_rf_cv, "tables/clay/Brazil/visnir_rf_cv.csv")
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/Brazil/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(visnir_rf_importance, "tables/clay/Brazil/visnir_rf_importance.csv")
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/clay/Brazil/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "Vis-NIR", Model = "RF",
           n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
           RMSE_cv = mean_from_folds(visnir_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(visnir_rf_cv, type = "R2"),
           RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
           R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs))

## PXRF models #####################################################################################
pxrf_data_br <- pxrf_data %>%
   filter(country == "Brazil")
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_br[-c(1, 2)], y = pxrf_data_br[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_br_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_br$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_br %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data_br %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/Brazil/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/Brazil/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/Brazil/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/Brazil/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/Brazil/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
pv_data_br <- pv_data %>%
   filter(country == "Brazil")
set.seed(100)
partition_index <- createDataPartition(pv_data_br$clay, p = 0.8, list = F)
pv_train_data <- pv_data_br %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(partition_index)
pv_valid_data <- pv_data_br %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(-partition_index)
mtry_param <- round((ncol(pv_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pv_valid_data$clay) %>%
   add_column(country = pv_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pv_rf_valid, "tables/clay/Brazil/pv_rf_valid.csv")
write_excel_csv(pv_rf_cv, "tables/clay/Brazil/pv_rf_cv.csv")
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/Brazil/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_replace(variables, "^band_", "")) %>%
   mutate(variables = str_replace(variables, "_nm$", "")) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pv_rf_importance, "tables/clay/Brazil/pv_rf_importance.csv")
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/clay/Brazil/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
   add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
           n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
           RMSE_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pv_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
           R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/Brazil/model_scores.csv")

# US ###############################################################################################
## Vis-NIR models ##################################################################################
## Partitioning
visnir_data_us <- visnir_data %>%
   filter(country == "US")
set.seed(100)
partition_index <- createDataPartition(visnir_data_us$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data_us %>%
   slice(partition_index)
visnir_valid_data <- visnir_data_us %>%
   slice(-partition_index)
mtry_param <- round((ncol(visnir_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## Models
## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
   as_tibble() %>%
   add_column(obs = visnir_valid_data$clay) %>%
   add_column(country = visnir_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(visnir_rf_valid, "tables/clay/US/visnir_rf_valid.csv")
write_excel_csv(visnir_rf_cv, "tables/clay/US/visnir_rf_cv.csv")
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/US/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(visnir_rf_importance, "tables/clay/US/visnir_rf_importance.csv")
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/clay/US/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "Vis-NIR", Model = "RF",
           n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
           RMSE_cv = mean_from_folds(visnir_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(visnir_rf_cv, type = "R2"),
           RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
           R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs))

## PXRF models #####################################################################################
pxrf_data_us <- pxrf_data %>%
   filter(country == "US")
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_us[-c(1, 2)], y = pxrf_data_us[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_us_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_us$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_us %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data_us %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/US/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/US/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/US/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/US/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/US/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
pv_data_us <- pv_data %>%
   filter(country == "US")
set.seed(100)
partition_index <- createDataPartition(pv_data_us$clay, p = 0.8, list = F)
pv_train_data <- pv_data_us %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(partition_index)
pv_valid_data <- pv_data_us %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(-partition_index)
mtry_param <- round((ncol(pv_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pv_valid_data$clay) %>%
   add_column(country = pv_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pv_rf_valid, "tables/clay/US/pv_rf_valid.csv")
write_excel_csv(pv_rf_cv, "tables/clay/US/pv_rf_cv.csv")
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/US/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_replace(variables, "^band_", "")) %>%
   mutate(variables = str_replace(variables, "_nm$", "")) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pv_rf_importance, "tables/clay/US/pv_rf_importance.csv")
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/clay/US/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
   add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
           n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
           RMSE_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pv_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
           R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/US/model_scores.csv")

# France ##############################################################################################
## Vis-NIR models ##################################################################################
## Partitioning
visnir_data_fr <- visnir_data %>%
   filter(country == "France")
set.seed(100)
partition_index <- createDataPartition(visnir_data_fr$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data_fr %>%
   slice(partition_index)
visnir_valid_data <- visnir_data_fr %>%
   slice(-partition_index)
mtry_param <- round((ncol(visnir_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## Models
## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
   as_tibble() %>%
   add_column(obs = visnir_valid_data$clay) %>%
   add_column(country = visnir_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(visnir_rf_valid, "tables/clay/France/visnir_rf_valid.csv")
write_excel_csv(visnir_rf_cv, "tables/clay/France/visnir_rf_cv.csv")
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/France/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(visnir_rf_importance, "tables/clay/France/visnir_rf_importance.csv")
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/clay/France/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "Vis-NIR", Model = "RF",
           n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
           RMSE_cv = mean_from_folds(visnir_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(visnir_rf_cv, type = "R2"),
           RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
           R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs))

## PXRF models #####################################################################################
pxrf_data_fr <- pxrf_data %>%
   filter(country == "France")
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_fr[-c(1, 2)], y = pxrf_data_fr[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_fr_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_fr$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_fr %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data_fr %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/France/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/France/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/France/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/France/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/France/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
pv_data_fr <- pv_data %>%
   filter(country == "France")
set.seed(100)
partition_index <- createDataPartition(pv_data_fr$clay, p = 0.8, list = F)
pv_train_data <- pv_data_fr %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(partition_index)
pv_valid_data <- pv_data_fr %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(-partition_index)
mtry_param <- round((ncol(pv_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pv_valid_data$clay) %>%
   add_column(country = pv_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pv_rf_valid, "tables/clay/France/pv_rf_valid.csv")
write_excel_csv(pv_rf_cv, "tables/clay/France/pv_rf_cv.csv")
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/France/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_replace(variables, "^band_", "")) %>%
   mutate(variables = str_replace(variables, "_nm$", "")) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pv_rf_importance, "tables/clay/France/pv_rf_importance.csv")
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/clay/France/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
   add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
           n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
           RMSE_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pv_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
           R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/France/model_scores.csv")

# India ############################################################################################
## Vis-NIR models ##################################################################################
## Partitioning
visnir_data_in <- visnir_data %>%
   filter(country == "India")
set.seed(100)
partition_index <- createDataPartition(visnir_data_in$clay, p = 0.8, list = F)
visnir_train_data <- visnir_data_in %>%
   slice(partition_index)
visnir_valid_data <- visnir_data_in %>%
   slice(-partition_index)
mtry_param <- round((ncol(visnir_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## Models
## RF
set.seed(100)
visnir_rf_model <- train(clay ~ ., data = visnir_train_data[-1], method = "rf",
                         preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
visnir_rf_cv <- visnir_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = visnir_train_data$country)
visnir_rf_valid <- predict(visnir_rf_model, newdata = visnir_valid_data) %>%
   as_tibble() %>%
   add_column(obs = visnir_valid_data$clay) %>%
   add_column(country = visnir_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(visnir_rf_valid, "tables/clay/India/visnir_rf_valid.csv")
write_excel_csv(visnir_rf_cv, "tables/clay/India/visnir_rf_cv.csv")
visnir_rf_plots <- validation_plot(visnir_rf_cv, visnir_rf_valid,
                                   variable = "clay", dataset = "Vis-NIR", model = "RF")
ggarrange(plotlist = visnir_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/India/visnir_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
visnir_rf_importance <- varImp(visnir_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_extract(variables, "\\d+")) %>% # extracting numbers
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(visnir_rf_importance, "tables/clay/India/visnir_rf_importance.csv")
ggplot(visnir_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("Vis-NIR - RF Variable Importance")
ggsave("figures/clay/India/visnir_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "Vis-NIR", Model = "RF",
           n_cv = nrow(visnir_rf_cv), n_valid = nrow(visnir_rf_valid),
           RMSE_cv = mean_from_folds(visnir_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(visnir_rf_cv, type = "R2"),
           RMSE_valid = RMSE(visnir_rf_valid$pred, visnir_rf_valid$obs),
           R2_valid = caret::R2(visnir_rf_valid$pred, visnir_rf_valid$obs))

## PXRF models #####################################################################################
pxrf_data_in <- pxrf_data %>%
   filter(country == "India")
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data_in[-c(1, 2)], y = pxrf_data_in[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_in_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data_in$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data_in %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data_in %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/India/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/India/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF", model = "RF")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/India/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/India/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/India/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- model_scores %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))

## PXRF + Vis-NIR (pv) models ######################################################################
## Partitioning
pv_data_in <- pv_data %>%
   filter(country == "India")
set.seed(100)
partition_index <- createDataPartition(pv_data_in$clay, p = 0.8, list = F)
pv_train_data <- pv_data_in %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(partition_index)
pv_valid_data <- pv_data_in %>%
   select(country, clay, all_of(pxrf_selected_predictors),
          "band_350_nm":"band_2500_nm") %>%
   slice(-partition_index)
mtry_param <- round((ncol(pv_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pv_rf_model <- train(clay ~ ., data = pv_train_data[-1], method = "rf",
                     preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pv_rf_cv <- pv_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pv_train_data$country)
pv_rf_valid <- predict(pv_rf_model, newdata = pv_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pv_valid_data$clay) %>%
   add_column(country = pv_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pv_rf_valid, "tables/clay/India/pv_rf_valid.csv")
write_excel_csv(pv_rf_cv, "tables/clay/India/pv_rf_cv.csv")
pv_rf_plots <- validation_plot(pv_rf_cv, pv_rf_valid,
                               variable = "clay", dataset = "PXRF + Vis-NIR", model = "RF")
ggarrange(plotlist = pv_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/India/pv_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pv_rf_importance <- varImp(pv_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 30, order_by = Overall, with_ties = F) %>% # 30 biggest values
   arrange(Overall) %>%
   mutate(variables = str_replace(variables, "^band_", "")) %>%
   mutate(variables = str_replace(variables, "_nm$", "")) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pv_rf_importance, "tables/clay/India/pv_rf_importance.csv")
ggplot(pv_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF + Vis-NIR - RF Variable Importance")
ggsave("figures/clay/India/pv_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

model_scores <- model_scores %>%
   add_row(Dataset = "PXRF + Vis-NIR", Model = "RF",
           n_cv = nrow(pv_rf_cv), n_valid = nrow(pv_rf_valid),
           RMSE_cv = mean_from_folds(pv_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pv_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pv_rf_valid$pred, pv_rf_valid$obs),
           R2_valid = caret::R2(pv_rf_valid$pred, pv_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/India/model_scores.csv")

# PXRF for all countries with Africa data ##########################################################
pxrf_data <- raw_data %>%
   select(country, clay, "K":"Pb") %>%
   mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat NAs as 0
   drop_na() %>%
   select(country, clay, "K":"Pb")

## PXRF models #####################################################################################
## Feature selection using RF
control_rfe <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
set.seed(100)
pxrf_feat_select <- rfe(x = pxrf_data[-c(1, 2)], y = pxrf_data[["clay"]],
                        sizes = c(1:16), rfeControl = control_rfe)
pxrf_selected_predictors <- predictors(pxrf_feat_select)
clay_allcountriesAfrica_pxrf_vars <- pxrf_selected_predictors

## Partitioning
set.seed(100)
partition_index <- createDataPartition(pxrf_data$clay, p = 0.8, list = F)
pxrf_train_data <- pxrf_data %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(partition_index)
pxrf_valid_data <- pxrf_data %>%
   select(country, clay, all_of(pxrf_selected_predictors)) %>%
   slice(-partition_index)
mtry_param <- round((ncol(pxrf_train_data) - 2) / 3) # var number / 3
grid <- expand.grid(mtry = mtry_param)

## RF
set.seed(100)
pxrf_rf_model <- train(clay ~ ., data = pxrf_train_data[-1], method = "rf",
                       preProcess = preprocess, trControl = control, tuneGrid = grid)

## RF results - kfold cross validation (80%) and hold-out validation (20%)
pxrf_rf_cv <- pxrf_rf_model$pred %>%
   filter(mtry == mtry_param) %>% # best model
   arrange(rowIndex) %>%
   add_column(country = pxrf_train_data$country)
pxrf_rf_valid <- predict(pxrf_rf_model, newdata = pxrf_valid_data) %>%
   as_tibble() %>%
   add_column(obs = pxrf_valid_data$clay) %>%
   add_column(country = pxrf_valid_data$country) %>%
   rename(pred = value)
write_excel_csv(pxrf_rf_valid, "tables/clay/allcountries_Africa/pxrf_rf_valid.csv")
write_excel_csv(pxrf_rf_cv, "tables/clay/allcountries_Africa/pxrf_rf_cv.csv")
pxrf_rf_plots <- validation_plot(pxrf_rf_cv, pxrf_rf_valid,
                                 variable = "clay", dataset = "PXRF",
                                 model = "RF", group_by = "country")
ggarrange(plotlist = pxrf_rf_plots, ncol = 2, common.legend = T, legend = "bottom")
ggsave("figures/clay/allcountries_Africa/pxrf_rf_pred_obs.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## RF importance
pxrf_rf_importance <- varImp(pxrf_rf_model)$importance %>%
   rownames_to_column("variables") %>%
   slice_max(n = 10, order_by = Overall, with_ties = F) %>% # 10 biggest values
   arrange(Overall) %>%
   mutate(variables = factor(variables,levels = variables))
write_excel_csv(pxrf_rf_importance, "tables/clay/allcountries_Africa/pxrf_rf_importance.csv")
ggplot(pxrf_rf_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
   ggtitle("PXRF - RF Variable Importance")
ggsave("figures/clay/allcountries_Africa/pxrf_rf_importance.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Combining model scores
model_scores <- tibble(Dataset = character(), Model = character(),
                       n_cv = numeric(), n_valid = numeric(),
                       RMSE_cv = numeric(), RMSE_valid = numeric(),
                       R2_cv = numeric(), R2_valid = numeric()) %>%
   add_row(Dataset = "PXRF", Model = "RF",
           n_cv = nrow(pxrf_rf_cv), n_valid = nrow(pxrf_rf_valid),
           RMSE_cv = mean_from_folds(pxrf_rf_cv, type = "RMSE"),
           R2_cv = mean_from_folds(pxrf_rf_cv, type = "R2"),
           RMSE_valid = RMSE(pxrf_rf_valid$pred, pxrf_rf_valid$obs),
           R2_valid = caret::R2(pxrf_rf_valid$pred, pxrf_rf_valid$obs))
write_excel_csv(model_scores, "tables/clay/allcountries_Africa/model_scores.csv")

## PXRF vars
best_pxrf_vars <- tibble(Dataset = c("All countries", "Brazil", "US",
                                     "France", "India", "All countries with Africa"),
                         Variables = c(paste(clay_allcountries_pxrf_vars, collapse = ", "),
                                       paste(clay_br_pxrf_vars, collapse = ", "),
                                       paste(clay_us_pxrf_vars, collapse = ", "),
                                       paste(clay_fr_pxrf_vars, collapse = ", "),
                                       paste(clay_in_pxrf_vars, collapse = ", "),
                                       paste(clay_allcountriesAfrica_pxrf_vars, collapse = ", ")),
                         n = c(length(clay_allcountries_pxrf_vars),
                               length(clay_br_pxrf_vars),
                               length(clay_us_pxrf_vars),
                               length(clay_fr_pxrf_vars),
                               length(clay_in_pxrf_vars),
                               length(clay_allcountriesAfrica_pxrf_vars)))
write_excel_csv(best_pxrf_vars, "tables/clay/best_pxrf_vars.csv")
