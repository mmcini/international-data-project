# Libs and functions ###############################################################################

source("R/_functions.R")

raw_data <- read_excel("data/oc_texture_data.xlsx", na = "NA")

# Target vars and land use #########################################################################
boxplot_data <- raw_data %>%
                rename(Clay = clay, Silt = silt, Sand = sand)
variable_names <- boxplot_data %>%
                  select(country, c("Clay":"Silt"), OC) %>%
                  colnames()
target_vars_boxplots <- create_boxplots(boxplot_data, variable_names, by_country = T, nrow = 1,
                        ordered = T, levels = c("OC", "Sand", "Silt", "Clay")) +
                        ggtitle("Boxplots of target variables (%)") +
                        theme(panel.grid = element_blank())

landuse_data <- boxplot_data %>%
                group_by(land_use, country) %>%
                drop_na(OC) %>%
                summarize(across("OC", mean, na.rm = T), n = n())

landuse_plot <- ggplot(landuse_data, aes(x = country, y = OC, fill = land_use)) +
                geom_bar(stat = "identity", position = "dodge", color = "black") +
                         geom_text(aes(label = paste0("n = ", n)),
                                   hjust = -0.2, vjust = 0, angle = 90,
                                   position = position_dodge(width = 0.9),
                                   family = "Times New Roman", size = 3) +
                         ggtitle("Organic carbon per landuse") +
                         ylab("Organic carbon (%)") + xlab("") +
                         scale_y_continuous(expand = expansion(mult = c(0.1, 0.4))) +
                         scale_fill_manual(values = brewer.pal(n = 5, "Spectral")) +
                         theme_bw() +
                         theme(text = element_text(family = "Times New Roman"),
                               legend.title = element_blank(),
                               panel.grid = element_blank())

target_vars_boxplots + landuse_plot + plot_layout(nrow = 2, ncol = 1)
ggsave("figures/publication/target_vars_boxplots_bycountries.png", dpi = 300, units = "mm",
       width = 230, height = 150, bg = "white")

# Obs pred plots ###################################################################################
## models that had the best results
oc_allcountries_valid <- read_csv("tables/OC/allcountries/pv_rf_valid.csv") %>%
                         mutate(country = factor(country))
oc_allcountries_cv <- read_csv("tables/OC/allcountries/pv_rf_cv.csv") %>%
                      mutate(country = factor(country))
best_oc_plots <- validation_plot(oc_allcountries_cv, oc_allcountries_valid,
                                 variable = "OC", dataset = "PXRF + Vis-NIR",
                                 model = "RF", group_by = "country")

sand_allcountries_valid <- read_csv("tables/sand/allcountries/pxrf_rf_valid.csv") %>%
                         mutate(country = factor(country))
sand_allcountries_cv <- read_csv("tables/sand/allcountries/pxrf_rf_cv.csv") %>%
                      mutate(country = factor(country))
best_sand_plots <- validation_plot(sand_allcountries_cv, sand_allcountries_valid,
                                   variable = "sand", dataset = "PXRF",
                                   model = "RF", group_by = "country")

silt_allcountries_valid <- read_csv("tables/silt/allcountries/pxrf_rf_valid.csv") %>%
                         mutate(country = factor(country))
silt_allcountries_cv <- read_csv("tables/silt/allcountries/pxrf_rf_cv.csv") %>%
                      mutate(country = factor(country))
best_silt_plots <- validation_plot(silt_allcountries_cv, silt_allcountries_valid,
                                   variable = "silt", dataset = "PXRF",
                                   model = "RF", group_by = "country")

clay_allcountries_valid <- read_csv("tables/clay/allcountries/pxrf_rf_valid.csv") %>%
                         mutate(country = factor(country))
clay_allcountries_cv <- read_csv("tables/clay/allcountries/pxrf_rf_cv.csv") %>%
                      mutate(country = factor(country))
best_clay_plots <- validation_plot(clay_allcountries_cv, clay_allcountries_valid,
                                   variable = "clay", dataset = "PXRF",
                                   model = "RF", group_by = "country")

ggarrange(plotlist = c(best_oc_plots, best_sand_plots, best_silt_plots, best_clay_plots),
          ncol = 2, nrow = 4, common.legend = T, legend = "bottom")
ggsave("figures/publication/best_pred_allcountries.png", dpi = 300, units = "mm",
       width = 210, height = 297, bg = "white")

# Importance plots #################################################################################
## models that had the best results
imp_plots <- list()

oc_importance <- read_csv("tables/OC/allcountries/pv_rf_importance.csv") %>%
                 slice_max(n = 15, order_by = Overall, with_ties = F) %>%
                 arrange(Overall) %>%
                 mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[1]] <- ggplot(oc_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("PXRF + Vis-NIR - OC")

sand_importance <- read_csv("tables/sand/allcountries/pxrf_rf_importance.csv") %>%
                   mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[2]] <- ggplot(sand_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("PXRF - Sand")

silt_importance <- read_csv("tables/silt/allcountries/pxrf_rf_importance.csv") %>%
                 mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[3]] <- ggplot(silt_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("PXRF - Silt")

clay_importance <- read_csv("tables/clay/allcountries/pxrf_rf_importance.csv") %>%
                 mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[4]] <- ggplot(clay_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("PXRF - Clay")

ggarrange(plotlist = imp_plots, ncol = 2, nrow = 2)
ggsave("figures/publication/var_importances.png", dpi = 300, units = "mm",
       width = 150, height = 150, bg = "white")
