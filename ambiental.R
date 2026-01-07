#FJ 2025

rm(list=ls()) # borra el ambiente
graphics.off() # Limpiar la lista de graficos

library(raster)
library(ggplot2)
library(ggrepel)

#Set files folder
setwd("D:/Trabajo/Paper/uruguay/ambiental/rasters/")
setwd("/media/oxus/LS/Trabajo/Paper/uruguay/ambiental")

#Import ambiental data - average data for the period 2000-2018/19/20. Data from: https://bio-oracle.org/
temp <- raster("temp_baseline_2000_2019_depthsurf_74ff_39fa_9adc_U1743720863576.nc")
salinity <- raster("salinity_baseline_2000_2019_depthsurf_49ed_5fc4_602c_U1743721486412.nc")
no3 <- raster("no3_baseline_2000_2018_depthsurf_8486_b388_df7c_U1743721538633.nc")
po4 <- raster("po4_baseline_2000_2018_depthsurf_6006_d51b_00e9_U1743721574986.nc")
si <- raster("si_baseline_2000_2018_depthsurf_395f_f84b_becc_U1743721600925.nc")
sithick <- raster("sithick_baseline_2000_2020_depthsurf_4bd1_0c3d_abde_U1743721639802.nc")
siconc <- raster("siconc_baseline_2000_2020_depthsurf_8197_c5be_1634_U1743721642133.nc")
o2 <- raster("o2_baseline_2000_2018_depthsurf_c6cd_716c_f3fc_U1743721626862.nc")
ph <- raster("ph_baseline_2000_2018_depthsurf_f606_6dc8_6180_U1743721635178.nc")
dfe <- raster("dfe_baseline_2000_2018_depthsurf_f8b5_f794_fbd9_U1743721632985.nc")
chl <- raster("chl_baseline_2000_2018_depthsurf_91d8_fb73_4955_U1743721638750.nc")
rad <- raster("rad_mean_baseline_2000_2020_depthsurf_08fa_60bc_f7a7_U1743996718946.nc")
cloud <- raster("cloud_baseline_2000_2020_depthsurf_4269_0d23_17d7_U1743996678983.nc")
air <- raster("air_baseline_2000_2020_depthsurf_0103_d100_33a3_U1743996692937.nc")

print(extent(chl))
print(res(chl))
print(crs(chl))

#Joint ambiental data in one file
rasStack = stack(temp,salinity,no3,po4,si,sithick,siconc,o2,ph,dfe,chl,rad,cloud,air)

#coordinates
pointCoordinates=read.csv("pointfile.txt",sep = ",", dec= ".")

#check coordinates mapping them

world <- rworldmap::getMap(resolution = "high")
uru <- as(extent(-77, -51, -58, -33), "SpatialPolygons")
proj4string(uru) <- CRS(proj4string(world))
world_clip <- raster::intersect(world,uru)
world_clip_f <- fortify(world_clip)

ggplot() + 
  geom_polygon(data = world_clip_f,aes(x = long, y = lat, group = group),fill = "#E5E5E5", colour = "black") +
  geom_label_repel(data = pointCoordinates,aes(label=place, x = longitude, y = latitude),max.overlaps = 20)

#convert coordinates into a spatial points dataframe
coordinates(pointCoordinates)= ~ latitude+longitude

#extract raster value by points
rasValue=extract(rasStack, pointCoordinates)

#combine rasters values with locations
combinePointValue=cbind(pointCoordinates,rasValue)

#save values as txt file
write.table(combinePointValue,file="data_ambiental.txt", append=FALSE, sep= ",",dec=".", row.names = FALSE, col.names=TRUE)