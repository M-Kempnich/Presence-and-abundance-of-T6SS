library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

#Import world map info from rnaturalearth (https://cran.r-project.org/web/packages/rnaturalearth/README.html)
world <- ne_coastline(scale = "medium", returnclass = "sf")
class(world)

#Read coordinates into data frame from csv file
lat_long<-read.csv("lat-long.csv",row.names=1,check.names=FALSE)

#create map using world data
p<-ggplot(data = world) +
  geom_sf()+
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("TARA Oceans T6SS")+
  #add a layer of points based on the lat_long data frame
  geom_point(data = lat_long, aes(x=Mean_Long,y=Mean_Lat,colour=RPM))+
  scale_color_gradient(low="blue",high="red")
  

png("Tara_map.png", width=6, height=4, units="in", res=600)
print(p)
dev.off()

pdf("Tara_map.pdf")
print(p)
dev.off()