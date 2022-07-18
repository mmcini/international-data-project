# Libs #############################################################################################

library(RColorBrewer)
library(tidyverse)
library(prospectr)
library(extrafont)
library(corrplot)
library(ggpubr)
library(readxl)

# Merging Vis-NIR data/fixing format　##############################################################

## Some data from Brazil has 2151 rows
## but most data was binned in intervals of 10 (350, 360 ... 2500)
unbinned_oc_data <- read_excel("../data/_unbinned_oc_data.xlsx", na = "NA")
unbinned_oc_spectra <- unbinned_oc_data %>%
                       select(c("355":"2500"))

## Creating bins between 360 and 2498 and adding columns 350 and 2500
## This is the format of the rest of the data
binned_oc_spectra <- binning(unbinned_oc_spectra, bin.size = 10) %>%
                     as_tibble()
binned_oc_spectra <- binned_oc_spectra[-215] %>% # removing band 2498 and replacing for 2500
                     add_column("350" = unbinned_oc_data$`350`, .before = "360") %>%
                     add_column("2500" = unbinned_oc_data$`2500`)

write_csv(binned_oc_spectra, "../data/_binned_oc_data.csv")

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
                select(variables) %>%
                colnames()
  }
  plots <- list()
  for (i in seq_len(length(variables))) {
    plots[[i]] <- ggplot(data, aes_string(x = variables[i], fill = group_by)) +
                  geom_histogram() + hist_layout
  }
  return(plots)
}

## Vis-NIR
visnir_layout <- list(theme_bw(), ylab("Reflectance factor"), xlab("Wavelength (nm)"),
                      scale_x_continuous(breaks = seq(0, 2500, 250)),
                      theme(text = element_text(family = "Times New Roman"),
                            legend.title = element_blank()))
visnir_cor_layout <- list(ylab(""), xlab("Wavelength (nm)"),
                          scale_x_continuous(breaks = seq(350, 2500, 150)),
                          scale_fill_distiller(limits = c(-1, 1), type = "div",
                                               palette = "RdBu", aesthetics = "fill"),
                          theme_bw(),
                          theme(text = element_text(family = "Times New Roman"),
                                panel.grid = element_blank(),
                                axis.text.x = element_text(angle = 90),
                                legend.title = element_blank()))

# Descriptive stats　###############################################################################

raw_oc_data <- read_excel("../data/oc_data.xlsx", na = "NA")

## All countries ###################################################################################
## OC plots
ggplot(raw_oc_data, aes(x = OC, fill = country)) +
  geom_histogram() + hist_layout

## PXRF plots
pxrf_data <- raw_oc_data %>%
             select(country, c(K:Pb))
pxrf_hists <- pxrf_histograms(pxrf_data, c(2:17), "country")
ggarrange(plotlist = pxrf_hists, ncol = 4, nrow = 4, common.legend = T, legend = "bottom")

## Vis-NIR plots
## No continuum removal
visnir_data <- raw_oc_data %>%
               select(country, c("350":"2500")) %>%
               filter(country != "Africa") %>% # Africa has no Vis-NIR data
               group_by(country) %>%
               summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
               pivot_longer(c("350":"2500"), names_to = "wavelength", values_to = "reflectance") %>%
               mutate(wavelength = as.numeric(wavelength))
ggplot(visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
       geom_line() + visnir_layout

## With continuum removal
cr_visnir_data <- raw_oc_data %>%
                  select(country, c("350":"2500")) %>%
                  filter(country != "Africa") %>% # Africa has no Vis-NIR data
                  group_by(country) %>%
                  summarise(across(c("350":"2500"), mean, na.rm = T)) %>%
                  select(c("350":"2500")) %>%
                  continuumRemoval() %>%
                  as.tibble()
colnames(cr_visnir_data) <- seq(350, 2500, 10)
cr_visnir_data <- add_column(cr_visnir_data, country = c("Brazil", "France", "India", "US"),
                             .before = "350") %>%
                             pivot_longer(c("350":"2500"), names_to = "wavelength",
                                          values_to = "reflectance") %>%
                             mutate(wavelength = as.numeric(wavelength))
ggplot(cr_visnir_data, aes(x = wavelength, y = reflectance, color = country)) +
       geom_line() + visnir_layout

## Correlations
## PXRF correlations
pxrf_cor <- raw_oc_data %>%
            select(OC, c(K:Pb)) %>%
            drop_na() %>%
            cor()
pxrf_cor_pvalue <- raw_oc_data %>%
                   select(OC, c(K:Pb)) %>%
                   drop_na() %>%
                   cor.mtest(conf.level = 0.95) # calculates p-value
par(family = "Times New Roman")
corrplot(pxrf_cor, method = "color", order = "hclust", p.mat = pxrf_cor_pvalue$p, addrect = 2,
         tl.col = "black", pch.cex = 1, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig",
         title = "PXRF vs. OC", mar = c(0, 0, 1.5, 0))

## Vis-NIR correlations
visnir_cor <- raw_oc_data %>%
              select(OC, c("350":"2500")) %>%
              drop_na() %>%
              cor() %>%
              as.tibble()
x_axis <- seq(350, 2500, 10)
ggplot(visnir_cor[c(2:217), 1], aes(x = x_axis, y = "OC", fill = OC)) +
      geom_tile() + visnir_cor_layout
