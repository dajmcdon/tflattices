# code to prepare `world-temperatures` dataset goes here
#
# original data downloaded from ECMWF as here
# https://github.com/akhodadadi/VolatilityTrend/blob/master/code/volatilitytrend/data_utils/ecmwf_dataUtils.py
# this dataset is only the Northern hemisphere for the year 2010 sampled at weekly frequency

# load("ecmf-weekly.Rdata")

library(tidyverse)
library(lubridate)

y2010 <- year(d) > 2009
temps_2010 <- raw[,y2010]
colnames(temps_2010) <- as.character(d[y2010])

world_temperatures <- as_tibble(temps_2010)
world_temperatures$longitude <- lons
world_temperatures$latitude <- lats

world_temperatures <- world_temperatures %>%
  relocate(longitude, latitude) %>%
  mutate(longitude = longitude - 360 * (longitude > 180))

usethis::use_data(world_temperatures, overwrite = TRUE)
