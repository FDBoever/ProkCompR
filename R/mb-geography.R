metadata = read.delim('~/DATA/MarbGenomics/testCAT.txt')
metadata = metadata[grepl('Marinobacter',metadata$Organism),]

ggplot(metadata,aes(lon,lat))+geom_point()



library(ggplot2)
library(dplyr)
require(maps)
require(viridis)

world <- map_data("world")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1
  ) +
  geom_point(
    data = metadata,
    aes(lon,lat,color= lifestyle),
    alpha = 0.7
  )
