

library(sf)


library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)

data_raw = fread("Glider data/cabot_20200717_114_delayed_corrected_v4.csv")

data_renamed = data_raw %>% rename(
    temp = sea_water_temperature,
    salt = sea_water_practical_salinity,
    density = sea_water_density,
    oxygen = micromoles_of_oxygen_per_unit_mass_in_sea_water
)


test = data_renamed[1:10000,]



sf_data_small = st_as_sf(test, coords = c("longitude", "latitude"), crs = 4326)

sf_data = st_as_sf(data_renamed, coords = c("longitude", "latitude"), crs = 4326)

plot(sf_data["time"], axes=TRUE)
