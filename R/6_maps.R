source("R/_functions.R")

# Data #############################################################################################
data <- read_xlsx("data/oc_texture_data.xlsx", na = "NA") %>%
        filter(!is.na(X))
points <- st_as_sf(data, coords = c("X", "Y"), crs = 4326)

world_map <- ne_countries(scale = "medium", returnclass = "sf")
studied_countries <- world_map %>%
                     filter(admin %in% c("Brazil", "India", "United States of America",
                                         "France", "Mozambique"))
studied_countries <- cbind(studied_countries, st_coordinates(st_centroid(studied_countries)))

### Warning: 10 Gb download!
clim_vars <- worldclim_global(var =  "bio", res = 0.5, path = "data/raster/bio/")
precip <- clim_vars$wc2.1_30s_bio_12
temp <- clim_vars$wc2.1_30s_bio_1
points <- cbind(points, extract(precip, vect(points))) %>%
          rename(precipitation = wc2.1_30s_bio_12)
points <- cbind(points, extract(temp, vect(points))) %>%
          rename(temperature = wc2.1_30s_bio_1)

# Maps #############################################################################################
## Samples
ggplot() +
       ylab("Latitude") + xlab("Longitude") +
       geom_sf(data = world_map) +
       geom_sf(data = studied_countries, fill = "steelblue") +
       geom_sf(data = points, aes(size = precipitation, color = temperature)) +
       geom_label_repel(data = studied_countries, aes(X, Y, label = admin), size = 3, force = 110,
                        family = "Times New Roman") +
       coord_sf(expand = F) +
       scale_size_continuous(range = c(1, 5)) +
       scale_color_distiller(palette = "RdBu") +
       guides(size = guide_legend(title = "Precipitation (mm)", title.position = "top",
                                  title.hjust = 0.5, label.hjust = 0.5),
              color = guide_colorbar(title = "Temperature (°C)", title.position = "top",
                                     title.hjust = 0.5, label.hjust = 0.5)) +
       theme_bw() +
       theme(text = element_text(family = "Times New Roman"),
            legend.position = "bottom",
            legend.key.size = unit(4, "mm"))

ggsave("figures/maps/samples.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")

# Climate vars vs. studied vars ####################################################################
clim_data_summary <- points %>%
                     dplyr::select(country, clay, silt, sand, OC, precipitation, temperature) %>%
                     mutate(across(c(clay:temperature), range_preprocess, na.rm = T)) %>%
                     group_by(country, levels = c("")) %>%
                     summarise(across(everything(), mean, na.rm = T)) %>%
                     pivot_longer(cols = clay:temperature,
                                  names_to = "variables",
                                  values_to = "values") %>%
                     mutate(variables = factor(variables, ordered = T,
                                               levels = c("clay", "silt", "sand",
                                                          "OC", "temperature",
                                                          "precipitation")))

color_palette <- c("#DC493A", "#467599", "#F5F749", "#000000", "#C14953", "#1D3354")
ggplot(data = clim_data_summary) +
       xlab("") + ylab("") +
       geom_bar(stat = "identity", position = "dodge",
                aes(x = country, y = values, fill = variables), color = "black", width = 0.7) +
      scale_fill_manual(name = "", values = color_palette,
                        labels = c("Clay", "Silt", "Sand", "OC", "Temperature", "Precipitation")) +
      theme_bw() +
      theme(text = element_text(family = "Times New Roman"),
            panel.grid = element_blank(),
            legend.key.size = unit(4, "mm"))

ggsave("figures/maps/climate_vars.png", dpi = 300, units = "mm",
       width = 150, height = 100, bg = "white")
