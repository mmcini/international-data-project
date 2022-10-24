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

clim_vars <- worldclim_global(var =  "bio", res = 0.5, path = "data/raster/bio/")
precip <- clim_vars$wc2.1_30s_bio_12
temp <- clim_vars$wc2.1_30s_bio_1
points <- merge(points, extract(precip, vect(points)))
points <- merge(points, extract(temp, vect(points)))

# Maps #############################################################################################
## Samples
ggplot() +
      ylab("Latitude") + xlab("Longitude") +
      geom_sf(data = world_map) +
      geom_sf(data = studied_countries, fill = "steelblue") +
      geom_sf(data = points, shape = 23, fill = "chocolate1") +
      geom_label_repel(data = studied_countries, aes(X, Y, label = admin), size = 3, force = 120,
                       family = "Times New Roman") +
      coord_sf(expand = F) +
      theme_bw() +
      theme(text = element_text(family = "Times New Roman"))

ggsave("figures/maps/samples.png", dpi = 300, units = "mm",
       width = 250, height = 150, bg = "white")
