# Libs #############################################################################################

library(tidyverse)
library(readxl)

# Merging Vis-NIR data/fixing format　##############################################################

## Some data from Brazil has 2151 rows
## but most data was binned in intervals of 10 (350, 360 ... 2500)
unbinned_oc_data <- read_excel("data/_unbinned_data.xlsx", na = "NA")
unbinned_oc_spectra <- unbinned_oc_data %>%
                       select(c("355":"2500"))

## Creating bins between 360 and 2498 and adding columns 350 and 2500
## This is the format of the rest of the data
binned_oc_spectra <- binning(unbinned_oc_spectra, bin.size = 10) %>%
                     as_tibble()
binned_oc_spectra <- binned_oc_spectra[-215] %>% # removing band 2498 and replacing for 2500
                     add_column("350" = unbinned_oc_data$`350`, .before = "360") %>%
                     add_column("2500" = unbinned_oc_data$`2500`)

write_csv(binned_oc_spectra, "data/_binned_data.csv")
