# 
# Description:
#
# Script to import environmental data form nc files and extract the specific values for each variable to each sampling coordinates
#

# Clean Rstudio environment
rm(list=ls()) 
graphics.off()

#Import libraries
library(ggplot2)  # ggplot() fortify()
library(raster)  # intersect()
library(rworldmap)  # getMap()
library(rworldxtra)
library(ggrepel)
library(patchwork)
library(rnaturalearth)

#Set files folder
setwd("D:/Trabajo/Paper/uruguay/mapas/")

# Get world map
world <- getMap(resolution = "high")

##### Zonas

#### Sudamerica

sudamerica <- ne_countries(scale = "medium", continent = "south america", returnclass = "sf")

lim_g1 <- c(xmin = -62, xmax = -53, ymin = -41, ymax = -33)
lim_g2 <- c(xmin = -78, xmax = -71, ymin = -53, ymax = -41)
lim_g3 <- c(xmin = -73, xmax = -58, ymin = -57, ymax = -51)

mapa_zonas_limpio <- ggplot() +
  # 1. Map of SudAmerica
  geom_sf(data = sudamerica, fill = "#E5E5E5", color = "grey80") +
  
  # 2. Squares for sampling zones
  geom_rect(aes(xmin = lim_g1["xmin"], xmax = lim_g1["xmax"], 
                ymin = lim_g1["ymin"], ymax = lim_g1["ymax"]), 
            color = "black", fill = NA, size = 0.6) +
  
  geom_rect(aes(xmin = lim_g2["xmin"], xmax = lim_g2["xmax"], 
                ymin = lim_g2["ymin"], ymax = lim_g2["ymax"]), 
            color = "black", fill = NA, size = 0.6) +
  
  geom_rect(aes(xmin = lim_g3["xmin"], xmax = lim_g3["xmax"], 
                ymin = lim_g3["ymin"], ymax = lim_g3["ymax"]), 
            color = "black", fill = NA, size = 0.6) +
  
  # 3. Lines to connect squeares with the specific maps
  annotate("segment", x = -53, xend = -10, y = -37, yend = -37, color = "black", size = 0.5) + # Rio de la Plata
  annotate("segment", x = -78, xend = -115, y = -47, yend = -47, color = "black", size = 0.5) + # South of Chile
  annotate("segment", x = -65, xend = -65, y = -57, yend = -85, color = "black", size = 0.5) + # Magallanes
  
  # 4. Map only with SudAmerica, no axis, no coordinates
  theme_void() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  
  # 5. Framing
  coord_sf(xlim = c(-90, -30), ylim = c(-60, 15), expand = FALSE, clip = "off")

# Show map
print(mapa_zonas_limpio)

ggsave(path = getwd(), filename = "sudamerica.tiff", width = 16, height = 8.3, device='tiff', dpi="print")

############# rio de la plata

myti_pop <-data.frame(pop=c("LPUY","PEUY","PBUY","PRUY","BSAR","BBAR","MPAR",
                            "Rio de la Plata", "Atlantic Ocean"),
                      Lat=c(-34.6687,-34.9639,-34.5950,-34.8817,-36.5400,-39.3750,-37.9427,-35.1141,-38.0000),
                      Lon=c(-54.1557,-54.9513,-54.1235,-55.2798,-56.1958,-60.5520,-57.5900,-56.9810,-54.3728))
# Preliminary plot
ggplot(myti_pop, mapping = aes(x = Lon, y = Lat)) + 
  geom_point(alpha = 0.5)

# Get map data
#extent
#vector (length=4; order= xmin, xmax, ymin, ymax)

clipper_chile <- as(extent(-63, -51, -42, -33), "SpatialPolygons")
proj4string(clipper_chile) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_chile)
world_clip_f <- fortify(world_clip)

#### Mapa con eventos geograficos

ggplot() + 
  geom_polygon(data = world_clip_f,aes(x = long, y = lat, group = group),fill = "#E5E5E5", colour = "black") +
  geom_label_repel(data = myti_pop[1:7,],aes(label=pop, x = Lon, y = Lat),
                   fill = "white",
                   alpha = 1, size = 4, fontface=2, #alpha manipula la transparencia, 1 es sin 
                   nudge_x = c(0.5,0.5,0.8,0,0,0,0,0), 
                   nudge_y = c(-0.2,-0.2,0.2,0,0,0,0,0,0),
  ) +
  geom_text_repel(data = myti_pop[8:9,],aes(label=pop, x = Lon, y = Lat),
                  alpha = 1, size = 4, fontface=2,segment.color = c(rep("black",2)), #alpha manipula la transparencia, 1 es sin 
                  nudge_x = c(1,0), 
                  nudge_y = c(-0.5,0),
  ) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x="longitud",y="Latitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("b) Sampling locations Rio de la Plata")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  coord_quickmap() +
  geom_point(data = myti_pop,aes(x = Lon, y = Lat),size=c(rep(6,7),rep(0,2)),
             color=c(colorRampPalette(c("#38F4E5","#64D5F7"))(4),"#208A8E","#238289","#267A85",
                                     "white","#E5E5E5")) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size=3))) +
  theme(legend.position = "none")

ggsave(path = getwd(), filename = "sampling_riodelaplata.tiff", width = 16, height = 8.3, device='tiff', dpi="print")

############ canales del sur

myti_pop <-data.frame(pop=c("CTCL","RMCL","PCCL","IPCL","Pacific Ocean"),
                      Lat=c(-42.9120,-43.7691,-44.7247,-50.8416,-45.1922),
                      Lon=c(-72.7187,-72.9512,-72.6886,-74.0115,-77.9859))

# Preliminary plot
ggplot(myti_pop, mapping = aes(x = Lon, y = Lat)) + 
  geom_point(alpha = 0.5)

# Get map data
#extent
#vector (length=4; order= xmin, xmax, ymin, ymax)

clipper_chile <- as(extent(-85, -72, -52, -41), "SpatialPolygons")
proj4string(clipper_chile) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_chile)
world_clip_f <- fortify(world_clip)

#### Map

ggplot() + 
  geom_polygon(data = world_clip_f,aes(x = long, y = lat, group = group),fill = "#E5E5E5", colour = "black") +
  geom_label_repel(data = myti_pop[1:4,],aes(label=pop, x = Lon, y = Lat),
                   fill = "white",
                   alpha = 1, size = 4, fontface=2, #alpha manipula la transparencia, 1 es sin 
                   nudge_x = c(0,0,0,0,0,0,0,0), 
                   nudge_y = c(-0,-0,0,0,0,0,0,0,0),
  ) +
  geom_text_repel(data = myti_pop[5,],aes(label=pop, x = Lon, y = Lat),
                  alpha = 1, size = 4, fontface=2,segment.color = "black", #alpha manipula la transparencia, 1 es sin 
                  nudge_x = c(0), 
                  nudge_y = c(0),
  ) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x="longitud",y="Latitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("a) Sampling locations \n South  of Chile")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  coord_quickmap() +
  geom_point(data = myti_pop,aes(x = Lon, y = Lat),size=c(rep(6,4),rep(0,1)),
             color=c(colorRampPalette(c("#EC7B74","#BF635A"))(4),"#E5E5E5")) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size=3))) +
  theme(legend.position = "none")

ggsave(path = getwd(), filename = "sampling_southofchile.tiff", width = 16, height = 8.3, device='tiff', dpi="print")

############ magallanes

myti_pop <-data.frame(pop=c("FKUK","PMAR","BACL","SGCL","BICL","BLCL","PACL","CMCL","ALAR",
                            "PPCL","CAAR","PNCL","PWCL",
                            "Falkland Islands","Isla Grande de Tierra del Fuego","Bahía Inutil",
                            "Strait of Magellan","Almirantazgo Sound","Beagle Channel","Atlantic Ocean","Pacific Ocean","Choiseul Sound","Navarino Island"),
                      Lat=c(-51.8266,-53.9468,-52.4941,-52.5608,-53.4727,-52.9694,-53.1651,-54.4812,-54.8689,
                            -52.8898,-54.8539,-54.9249,-54.9340,
                            -51.8002,-54.5000,-53.5306,-53.4881,-54.2949,-54.8745,-52.4224,-55.8089,-51.9185,-55.0285),
                      Lon=c(-58.9721,-67.4819,-69.5220,-70.0428,-69.3400,-70.8255,-70.9008,-68.9934,-67.5522,
                            -70.1226,-67.5033,-68.3248,-67.6146,
                            -60.0000,-66.9768,-69.8365,-70.7594,-69.5677,-68.3266,-64.0418,-72.0080,-58.7523,-67.6287))

# Preliminary plot
ggplot(myti_pop, mapping = aes(x = Lon, y = Lat)) + 
  geom_point(alpha = 0.5)

# Get map data
#extent
clipper_chile <- as(extent(-74, -57, -56, -51), "SpatialPolygons")
proj4string(clipper_chile) <- CRS(proj4string(world))
world_clip <- raster::intersect(world, clipper_chile)
world_clip_f <- fortify(world_clip)

#### Map

ggplot() + 
  geom_polygon(data = world_clip_f,aes(x = long, y = lat, group = group),fill = "#E5E5E5", colour = "black") +
  geom_label_repel(data = myti_pop[1:13,],aes(label=pop, x = Lon, y = Lat),
                   fill = "white",
                   alpha = 1, size = 4, fontface=2, #alpha manipula la transparencia, 1 es sin 
                   nudge_x = c(-1,1,1.5,-1,0.4,-0.2,-0.6,0.4,-0.5,0.3,0.2,0.4,1.5), 
                   nudge_y = c(-0.7,0.5,0,0.3,0.3,0.3,0.2,0.2,0.2,-0.15,0.2,-0.2,-0.3),
  ) +
  geom_text_repel(data = myti_pop[14:23,],aes(label=pop, x = Lon, y = Lat),
                   alpha = 1, size = 4, fontface=2,segment.color = c(NA,rep("black",9)), #alpha manipula la transparencia, 1 es sin 
                   nudge_x = c(0.5,4,0.4,-2.5,-3,-1,rep(0,2),0.5,1.5), 
                   nudge_y = c(0.7,0.6,-0.3,-1.2,-1,-1,rep(0,2),-0.8,-0.5),
  ) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  labs(x="longitud",y="Latitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("c) Sampling locations \n Strait of Magellan, Isla Grande de Tierra del Fuego and Falkland Islands")+
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) +
  coord_quickmap() +
  geom_point(aes(x = Lon, y = Lat),size=c(rep(6,13),rep(0,2),rep(0.5,4),rep(0,2),0.5,0.5),
             data = myti_pop,color=c("#297380","#2C6B7C","#306478",colorRampPalette(c("#5E8B48","#DAC753"))(10),"white","#E5E5E5",rep("black",4),rep("white",2),"black","black")) + 
  guides(color = guide_legend(override.aes = list(linetype = 0, size=3))) +
  theme(legend.position = "none")

ggsave(path = getwd(), filename = "sampling_patagonia.tiff", width = 16, height = 8.3, device='tiff', dpi="print")
