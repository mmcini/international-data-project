# Libs and functions ###############################################################################

source("R/_functions.R")

raw_data <- read_excel("data/oc_texture_data.xlsx", na = "NA")

# Target vars and land use #########################################################################
boxplot_data <- raw_data %>%
                mutate(country = factor(country)) %>%
                rename(Clay = clay, Silt = silt, Sand = sand)
boxplot_data$Sand[which(boxplot_data$country == "Mozambique")] <- NA
boxplot_data$Silt[which(boxplot_data$country == "Mozambique")] <- NA
boxplot_data$Clay[which(boxplot_data$country == "Mozambique")] <- NA
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
                         ggtitle("Organic carbon per land use") +
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
oc_allcountries_valid <- read_csv("tables/OC/allcountries/pxrf_rf_valid.csv") %>%
                         mutate(country = factor(country))
oc_allcountries_cv <- read_csv("tables/OC/allcountries/pxrf_rf_cv.csv") %>%
                      mutate(country = factor(country))
best_oc_plots <- validation_plot(oc_allcountries_cv, oc_allcountries_valid,
                                 variable = "OC", dataset = "PXRF",
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
# Important pXRF variables
## shows the number of times each variable appears as most important for each dataset
OC_imp_vars <- read_csv("tables/OC/best_pxrf_vars.csv")$Variables %>%
               extract_pxrf_vars() %>%
               table(dnn = c("pxrf_variables")) %>%
               as_tibble() %>%
               add_column(target_variable = "OC") %>%
               arrange(desc(n)) %>%
               slice_head(n = 10)
clay_imp_vars <- read_csv("tables/clay/best_pxrf_vars.csv")$Variables %>%
                 extract_pxrf_vars() %>%
                 table(dnn = c("pxrf_variables")) %>%
                 as_tibble() %>%
                 add_column(target_variable = "Clay") %>%
                 arrange(desc(n)) %>%
                 slice_head(n = 10)
silt_imp_vars <- read_csv("tables/silt/best_pxrf_vars.csv")$Variables %>%
                 extract_pxrf_vars() %>%
                 table(dnn = c("pxrf_variables")) %>%
                 as_tibble() %>%
                 add_column(target_variable = "Silt") %>%
                 arrange(desc(n)) %>%
                 slice_head(n = 10)
sand_imp_vars <- read_csv("tables/sand/best_pxrf_vars.csv")$Variables %>%
                 extract_pxrf_vars() %>%
                 table(dnn = c("pxrf_variables")) %>%
                 as_tibble() %>%
                 add_column(target_variable = "Sand") %>%
                 arrange(desc(n)) %>%
                 slice_head(n = 10)

imp_table <- rbind(OC_imp_vars, clay_imp_vars, silt_imp_vars, sand_imp_vars)

pxrf_imp_plots <- ggplot(data = imp_table) +
                  geom_bar(stat = "identity", width = 0.5, fill = "steelblue",
                           aes(x = reorder_within(pxrf_variables, -n, target_variable), y = n)) +
                  xlab("") + ylab("Number of occurrences") +
                  facet_wrap(.~target_variable, scales = "free_x") +
                  scale_x_reordered() +
                  theme_bw() +
                  theme(text = element_text(family = "Times New Roman"),
                        panel.grid = element_blank())

## models that had the best results
imp_plots <- list()

oc_importance <- read_csv("tables/OC/allcountries/visnir_rf_importance.csv") %>%
                 slice_max(n = 15, order_by = Overall, with_ties = F) %>%
                 arrange(Overall) %>%
                 mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[1]] <- ggplot(oc_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("Vis-NIR - OC") + theme(panel.grid = element_blank())

sand_importance <- read_csv("tables/sand/allcountries/visnir_rf_importance.csv") %>%
                   slice_max(n = 15, order_by = Overall, with_ties = F) %>%
                   arrange(Overall) %>%
                   mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[2]] <- ggplot(sand_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("Vis-NIR - Sand") + theme(panel.grid = element_blank())

silt_importance <- read_csv("tables/silt/allcountries/visnir_rf_importance.csv") %>%
                   slice_max(n = 15, order_by = Overall, with_ties = F) %>%
                   arrange(Overall) %>%
                   mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[3]] <- ggplot(silt_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("Vis-NIR - Silt") + theme(panel.grid = element_blank())

clay_importance <- read_csv("tables/clay/allcountries/visnir_rf_importance.csv") %>%
                   slice_max(n = 15, order_by = Overall, with_ties = F) %>%
                   arrange(Overall) %>%
                   mutate(variables = factor(variables,levels = variables, ordered = T))
imp_plots[[4]] <- ggplot(clay_importance, aes(x = variables, y = Overall)) + importance_plot_layout +
                  ggtitle("Vis-NIR - Clay") + theme(panel.grid = element_blank())

visnir_imp <- ggarrange(plotlist = imp_plots, ncol = 2, nrow = 2)

pxrf_imp_plots + visnir_imp + plot_layout(ncol = 2)
ggsave("figures/publication/var_importances.png", dpi = 300, units = "mm",
       width = 200, height = 130, bg = "white")

cor_visnir_allvariables <- raw_data %>%
                           select("country", "OC", "sand", "silt", "clay", c("350":"2500")) %>%
                           filter(country != "Africa") %>% # Africa has no Vis-NIR data
                           drop_na()

## Vis-NIR plots
## No continuum removal
visnir_data <- raw_data %>%
               select(country, c("350":"2500")) %>%
               filter(country != "Africa") %>% # Africa has no Vis-NIR data
               group_by(country) %>%
               summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
               pivot_longer(c("350":"2500"), names_to = "wavelength", values_to = "reflectance") %>%
               mutate(wavelength = as.numeric(wavelength))
visnir_plot <- ggplot(visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
               geom_line() + visnir_layout + ggtitle("Vis-NIR Original Spectra") + xlab("")
visnir_plot <- visnir_labels(visnir_plot)

## With continuum removal
cr_visnir_data <- raw_data %>%
                  select(country, c("350":"2500")) %>%
                  filter(country != "Mozambique") %>% # Africa has no Vis-NIR data
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
                  geom_line() + visnir_layout + ggtitle("Vis-NIR Continuum Removal") +
                  xlab("Wavelength (nm)")
cr_visnir_plot <- visnir_labels(cr_visnir_plot)

visnir_plot + cr_visnir_plot + plot_layout(nrow = 2, guides = "collect")
ggsave("figures/publication/visnir_cr_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 130, bg = "white")

x_axis <- seq(350, 2500, 10)
cor_matrix_allcountries <- cor_visnir_allvariables %>%
                           select(-country) %>%
                           cor() %>%
                           as_tibble() %>%
                           select("OC":"clay") %>%
                           slice(5:n()) %>%
                           add_column(x = x_axis) %>%
                           rename(clay = clay, silt = silt, sand = sand) %>%
                           pivot_longer(c("OC":"clay"),
                                        names_to = "variables", values_to = "values") %>%
                           mutate(variables = factor(variables,
                                                     labels = c("Clay", "Silt", "Sand", "OC"),
                                                     ordered = T))

cor_matrix_br <- cor_visnir_allvariables %>%
                 filter(country == "Brazil") %>%
                 select(-country) %>%
                 cor() %>%
                 as_tibble() %>%
                 select("OC":"clay") %>%
                 slice(5:n()) %>%
                 add_column(x = x_axis) %>%
                 rename(clay = clay, silt = silt, sand = sand) %>%
                 pivot_longer(c("OC":"clay"),
                              names_to = "variables", values_to = "values") %>%
                 mutate(variables = factor(variables,
                                           labels = c("Clay", "Silt", "Sand", "OC"),
                                           ordered = T))

cor_matrix_fr <- cor_visnir_allvariables %>%
                 filter(country == "France") %>%
                 select(-country) %>%
                 cor() %>%
                 as_tibble() %>%
                 select("OC":"clay") %>%
                 slice(5:n()) %>%
                 add_column(x = x_axis) %>%
                 rename(clay = clay, silt = silt, sand = sand) %>%
                 pivot_longer(c("OC":"clay"),
                              names_to = "variables", values_to = "values") %>%
                 mutate(variables = factor(variables,
                                           labels = c("Clay", "Silt", "Sand", "OC"),
                                           ordered = T))

cor_matrix_in <- cor_visnir_allvariables %>%
                 filter(country == "India") %>%
                 select(-country) %>%
                 cor() %>%
                 as_tibble() %>%
                 select("OC":"clay") %>%
                 slice(5:n()) %>%
                 add_column(x = x_axis) %>%
                 rename(clay = clay, silt = silt, sand = sand) %>%
                 pivot_longer(c("OC":"clay"),
                              names_to = "variables", values_to = "values") %>%
                 mutate(variables = factor(variables,
                                           labels = c("Clay", "Silt", "Sand", "OC"),
                                           ordered = T))

cor_matrix_us <- cor_visnir_allvariables %>%
                 filter(country == "US") %>%
                 select(-country) %>%
                 cor() %>%
                 as_tibble() %>%
                 select("OC":"clay") %>%
                 slice(5:n()) %>%
                 add_column(x = x_axis) %>%
                 rename(clay = clay, silt = silt, sand = sand) %>%
                 pivot_longer(c("OC":"clay"),
                              names_to = "variables", values_to = "values") %>%
                 mutate(variables = factor(variables,
                                           labels = c("Clay", "Silt", "Sand", "OC"),
                                           ordered = T))

visnir_corplot_allcountries <- ggplot(cor_matrix_allcountries,
                                      aes(x = x, y = variables, fill = values)) +
                  geom_tile() + cor_visnir_layout +
                  ggtitle("Pearson's coefficient - All countries (except Mozambique)")
visnir_corplot_allcountries <- visnir_cor_labels(visnir_corplot_allcountries)

visnir_corplot_br <- ggplot(cor_matrix_br,
                            aes(x = x, y = variables, fill = values)) +
                  geom_tile() + cor_visnir_layout +
                  ggtitle("Pearson's coefficient - Brazil")
visnir_corplot_br <- visnir_cor_labels(visnir_corplot_br)

visnir_corplot_fr <- ggplot(cor_matrix_fr,
                            aes(x = x, y = variables, fill = values)) +
                  geom_tile() + cor_visnir_layout +
                  ggtitle("Pearson's coefficient - France")
visnir_corplot_fr <- visnir_cor_labels(visnir_corplot_fr)

visnir_corplot_in <- ggplot(cor_matrix_in,
                            aes(x = x, y = variables, fill = values)) +
                  geom_tile() + cor_visnir_layout +
                  ggtitle("Pearson's coefficient - India")
visnir_corplot_in <- visnir_cor_labels(visnir_corplot_in)

visnir_corplot_us <- ggplot(cor_matrix_us,
                            aes(x = x, y = variables, fill = values)) +
                  geom_tile() +
                  xlab("Wavelength (nm)") + cor_visnir_layout +
                  ggtitle("Pearson's coefficient - US")
visnir_corplot_us <- visnir_cor_labels(visnir_corplot_us)

visnir_corplot_allcountries + visnir_corplot_br + visnir_corplot_fr +
visnir_corplot_in + visnir_corplot_us + plot_layout(nrow = 5, guides = "collect")

ggsave("figures/publication/visnir_corplots.png", dpi = 300, units = "mm",
       width = 170, height = 200, bg = "white")

tex_data <- raw_data %>%
            select(country, sand, silt, clay, "K":"Pb", "350":"2500") %>%
            mutate(across(c("K":"Pb"), ~ replace_na(., 0))) %>% # treat nas as 0
            drop_na() %>%
            mutate(country = factor(country,
                                    labels = c("Brazil", "France", "India", "US"))) %>%
            select(clay, sand, silt, country) %>%
            filter(country != "Mozambique") %>%
            rename(Sand = sand, Silt = silt, Clay = clay)
usda_tex_classes <- read_xlsx("data/tex_triangle/USDA.xlsx")

# Put tile labels at the midpoint of each tile
usda_labels <- usda_tex_classes %>%
               group_by(Label) %>%
               summarise(across(c(1:3), mean))

# Tweak
usda_labels$Angle = 0
usda_labels$Angle[which(usda_labels$Label == 'Loamy Sand')] = -35

# Construct the plot.
ggplot(data = usda_tex_classes, ggtern::aes(y = Clay, x = Sand, z = Silt)) +
       coord_tern(L = "x", T = "y", R = "z") +
       Llab("") + Tlab("") + Rlab("") + Larrowlab("Sand (%)") +
       Tarrowlab("Clay (%)") + Rarrowlab("Silt (%)") +
       geom_polygon(ggtern::aes(fill = Label), show.legend = F,
                    alpha = 0.000002, size = 0.6, colour = "black") +
       geom_point(data = tex_data, size = 4.5, alpha = 0.9,
                  ggtern::aes(colour = country)) +
       geom_text(data = usda_labels,
                 ggtern::aes(label = Label),
                 color = 'black',
                 size = 4.8,
                 family = "serif",
                 alpha = 1) +
       scale_color_manual(limits = levels(tex_data$country),
                          values = country_colors) +
       ggtern::theme_bw() + theme_clockwise() +
       theme(legend.position = "bottom",
             legend.title = element_blank(),
             axis.text = element_text(family = "serif", size = 15,
                                      color = "black"),
             text = element_text(family = "serif", size = 17,
                                 color = "black")) +
       theme_showarrows()

ggtern::ggsave("figures/publication/tex_plot.png", dpi = 300, units = "mm",
       width = 250, height = 250, bg = "white")
