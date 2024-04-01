#' Temperatures for the northern hemisphere in 2010
#'
#' This data represents 4-degree x 2-degree gridded temperature for the
#' northern hemisphere in 2010. Each week averages the daily data from the
#' original source.
#'
#' @format A 2880 row by 53 column tibble. Columns are `longitude`, `latitude`
#' (duplicated as necessary) and the 51 complete weeks in 2010. This represents
#' 80 longitudes (0 - 360 by 4) and 36 latitudes (0 - 80 by 2).
#'
#' @source This is a subset of subsampled ERA-20C data available at <https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-20th-century>
"world_temperatures"


#' Spatial polygon map of Canadian Provinces
#' 
#' Canadian provincial borders are rather complex (many islands and inlets).
#' The official map from StatCan is around 20 Mb and challenging to use to make
#' quick and dirty maps. This dataset fills the gap. It's a lower resolution 
#' version.
#' 
#' @format map of class SpatialPolygons
#' 
#' @examples
#' library(ggplot2)
#' ggplot() + 
#'   geom_path(data=canada_map,aes(x=long,y=lat,group=group),
#'   color="darkgreen", size = .5) +
#'   coord_map(projection = 'lambert', parameters=c(49,77)) +
#'   theme_bw()
"canada_map"