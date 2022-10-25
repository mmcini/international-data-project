# Libs and functions ###############################################################################

source("R/_functions.R")

# Descriptive stats　###############################################################################

raw_data <- read_excel("data/oc_texture_data.xlsx", na = "NA")

## All countries
desc_stats_allcountries <- raw_data %>%
                           select("OC", "sand", "silt", "clay", "K":"Pb") %>%
                           descriptive_stats()
write_excel_csv(desc_stats_allcountries, "tables/stats/descriptive_stats_all_countries.csv")

## By country
desc_stats_bycountries <- raw_data %>%
                          select("country", "OC", "sand", "silt", "clay", "K":"Pb") %>%
                          descriptive_stats(group_by = "country")
write_excel_csv(desc_stats_bycountries, "tables/stats/descriptive_stats_by_country.csv")

## Variable plots
histogram_data <- raw_data %>%
                  select("country", "OC", "sand", "silt", "clay")
target_vars_hists <- create_histograms(histogram_data, c(2:5), "country")
ggarrange(plotlist = target_vars_hists, ncol = 4, common.legend = T, legend = "bottom")
ggsave("figures/stats/target_variables_histograms.png", dpi = 300, units = "mm",
       width = 150, height = 50, bg = "white")

## Boxplots
variable_names <- raw_data %>%
                  select(c("K":"Pb")) %>%
                  colnames()
create_boxplots(raw_data, variable_names, by_country = F)
ggsave("figures/stats/pxrf_boxplots_allcountries.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")
variable_names <- raw_data %>%
                  select(country, c("K":"Pb")) %>%
                  colnames()
create_boxplots(raw_data, variable_names, by_country = T)
ggsave("figures/stats/pxrf_boxplots_bycountries.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

## PXRF plots
pxrf_data <- raw_data %>%
             select(country, c(K:Pb))
pxrf_hists <- create_histograms(pxrf_data, c(2:17), "country")
ggarrange(plotlist = pxrf_hists, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")
ggsave("figures/stats/pxrf_variables_histograms.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

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
               geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra")
visnir_labels(visnir_plot)
ggsave("figures/stats/visnir_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 100, bg = "white")

## With continuum removal
cr_visnir_data <- raw_data %>%
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
ggsave("figures/stats/visnir_cr_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 100, bg = "white")

## Savitzky-Golay + first derivative
deriv_visnir_data <- raw_data %>%
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
ggsave("figures/stats/visnir_deriv_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

## Correlations
## PXRF correlations
cor_pxrf <- raw_data %>%
            select("OC", "sand", "silt", "clay", c("K":"Pb")) %>%
            drop_na() %>%
            cor()
cor_pxrf_pvalue <- raw_data %>%
                   select("OC", "sand", "silt", "clay", c("K":"Pb")) %>%
                   drop_na() %>%
                   cor.mtest(conf.level = 0.95) # calculates p-value
png("figures/stats/pxrf_corrplot_all_countries.png", units = "mm", res = 300,
    width = 200, height = 200, bg = "white")
par(family = "Times New Roman")
corrplot::corrplot(cor_pxrf, method = "color", order = "hclust",
                   p.mat = cor_pxrf_pvalue$p, addrect = 2,
                   tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05),
                   insig = "label_sig", title = "PXRF vs. OC and Texture", mar = c(0, 0, 1.5, 0))
dev.off()

## Vis-NIR correlations
## OC
cor_visnir_plot_allcountries <- create_cor_visnir(raw_data,
                                                  bands = c(27:242), property = "OC")
cor_visnir_bycountry <- raw_data %>%
                        select("country", "OC", c("350":"2500")) %>%
                        filter(country != "Africa") %>% # Africa has no Vis-NIR data
                        drop_na()
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(3:218),
                                                property = "OC", group_by = "country")
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")
ggsave("figures/stats/visnir_corrplot_oc_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

# Sand
cor_visnir_plot_allcountries <- create_cor_visnir(raw_data,
                                                  bands = c(27:242), property = "sand")
cor_visnir_bycountry <- raw_data %>%
                        select("country", "sand", c("350":"2500")) %>%
                        filter(country != "Africa") %>% # Africa has no Vis-NIR data
                        drop_na()
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(3:218),
                                                property = "sand", group_by = "country")
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")
ggsave("figures/stats/visnir_corrplot_sand_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

# Silt
cor_visnir_plot_allcountries <- create_cor_visnir(raw_data,
                                                  bands = c(27:242), property = "silt")
cor_visnir_bycountry <- raw_data %>%
                        select("country", "silt", c("350":"2500")) %>%
                        filter(country != "Africa") %>% # Africa has no Vis-NIR data
                        drop_na()
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(3:218),
                                                property = "silt", group_by = "country")
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")
ggsave("figures/stats/visnir_corrplot_silt_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

# Clay
cor_visnir_plot_allcountries <- create_cor_visnir(raw_data,
                                                  bands = c(27:242), property = "clay")
cor_visnir_bycountry <- raw_data %>%
                        select("country", "clay", c("350":"2500")) %>%
                        filter(country != "Africa") %>% # Africa has no Vis-NIR data
                        drop_na()
cor_visnir_plot_bycountry <- create_cor_visnir(cor_visnir_bycountry, bands = c(3:218),
                                                property = "clay", group_by = "country")
cor_visnir_plot_list <- append(cor_visnir_plot_allcountries, cor_visnir_plot_bycountry)
ggarrange(plotlist = cor_visnir_plot_list,ncol = 1, nrow = 5, common.legend = T, legend = "bottom")
ggsave("figures/stats/visnir_corrplot_clay_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 200, bg = "white")

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

ggplot(data = imp_table) +
       geom_bar(stat = "identity", width = 0.5, fill = "steelblue",
                aes(x = reorder_within(pxrf_variables, -n, target_variable), y = n)) + xlab("") +
       facet_wrap(.~target_variable, scales = "free_x") +
       scale_x_reordered() +
       theme_bw() +
       theme(text = element_text(family = "Times New Roman"),
             panel.grid = element_blank())

ggsave("figures/stats/var_importance_count.png", dpi = 300, units = "mm",
       width = 110, height = 110, bg = "white")
