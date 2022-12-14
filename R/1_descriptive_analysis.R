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
create_boxplots(raw_data, variable_names, by_country = T) + theme(panel.grid = element_blank())
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
               geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra") + xlab("")
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
                  geom_line() + visnir_layout + ggtitle("Vis-NIR Spectra - Continuum Removal")
cr_visnir_plot <- visnir_labels(cr_visnir_plot)
ggarrange(visnir_plot, cr_visnir_plot, nrow = 2, common.legend = T, legend = "bottom")
ggsave("figures/stats/visnir_cr_spectra_all_countries.png", dpi = 300, units = "mm",
       width = 200, height = 150, bg = "white")

x_axis <- seq(350, 2500, 10)
cor_matrix <- cor_visnir_allvariables %>%
              select(-country) %>%
              cor() %>%
              as_tibble() %>%
              select("oc":"clay") %>%
              slice(5:n()) %>%
              add_column(x = x_axis) %>%
              rename(clay = clay, silt = silt, sand = sand) %>%
              pivot_longer(c("oc":"clay"), names_to = "variables", values_to = "values") %>%
              mutate(variables = factor(variables, c("clay", "silt", "sand", "oc"), ordered = t))

ggplot(cor_matrix, aes(x = x, y = variables, fill = values)) +
       geom_tile() + ylab("") + cor_visnir_layout + xlab("wavelength (nm)") +
       ggtitle("vis-nir - pearson's correlation coefficients")

## savitzky-golay + first derivative
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

# All variables using all data
cor_visnir_allvariables <- raw_data %>%
                           select("country", "OC", "sand", "silt", "clay", c("350":"2500")) %>%
                           filter(country != "Africa") %>% # Africa has no Vis-NIR data
                           drop_na()

properties <- c("OC", "sand", "silt", "clay")
cor_visnir_plot_list_allvariables <- list()
for (property in properties) {
  plot <- create_cor_visnir(cor_visnir_allvariables, bands = c(6:221), property = property)
  plot <- plot[[1]] + ylab(str_to_title(property)) + ggtitle("")
  cor_visnir_plot_list_allvariables[[property]] <- plot
}
cor_visnir_plot_list_allvariables[["OC"]] <- cor_visnir_plot_list_allvariables[["OC"]] +
                                             ylab("OC") + ggtitle("Vis-NIR correlations") +
                                             theme(axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   plot.margin = margin(r = 5,
                                                                        unit = "mm"))
cor_visnir_plot_list_allvariables[["sand"]] <- cor_visnir_plot_list_allvariables[["sand"]] +
                                             theme(axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   plot.margin = margin(r = 5,
                                                                        unit = "mm"))
cor_visnir_plot_list_allvariables[["silt"]] <- cor_visnir_plot_list_allvariables[["silt"]] +
                                             theme(axis.text.x = element_blank(),
                                                   axis.ticks.x = element_blank(),
                                                   axis.title.x = element_blank(),
                                                   plot.margin = margin(r = 5,
                                                                        unit = "mm"))
cor_visnir_plot_list_allvariables[["clay"]] <- cor_visnir_plot_list_allvariables[["clay"]] +
                                               xlab("Wavelength (nm)") +
                                               theme(plot.margin = margin(r = 5, b = 2,
                                                                          unit = "mm"))
ggarrange(plotlist = cor_visnir_plot_list_allvariables,ncol = 1,
          nrow = 4, common.legend = T, legend = "bottom")

ggsave("figures/stats/visnir_corrplot_clay_all_variables.png", dpi = 300, units = "mm",
       width = 130, height = 130, bg = "white")
