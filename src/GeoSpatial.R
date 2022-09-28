library(sf)
library(tmap)
library(tidyverse)
library(highcharter)
library(plotly)
library(ggmap)
library(maptools)
library(rgdal)
require(GISTools)

#Loading Required Dataset
resale <- read.csv("/Users/cavin/Desktop/GeoSpatial 2/resale.csv", header = T)
hdb_coord <- read.csv("/Users/cavin/Desktop/GeoSpatial 2/GeoCoded.csv", header = T)
mpsz <- st_read(dsn = "/Users/cavin/Desktop/GeoSpatial 2/master-plan-2014-subzone-boundary-web-shp", layer = "MP14_SUBZONE_WEB_PL")

#Filter 2021 data to be analyzed
resale2021 <- resale %>%
  filter(year == 2021) %>%
  select('town','flat_type','storey_range', 'remaining_lease', 'resale_price','block','street_name')

plotmean <- resale2021 %>%
  group_by(town) %>% 
  summarise(Mean_Resale_Price = mean(resale_price))

mpsz_resale2021 <- left_join(mpsz, plotmean, by = c("PLN_AREA_N" = "town"))

### GeoCoding
resale2021$add <- paste("Block ",resale2021$block," ",resale2021$street_name," Singapore")
# Get unique addresses
all_adds <- unique(resale2021$add)
# Scrape coordinates
register_google(key = "AIzaSyDmGZpQq_WOL7XdE0nMllwZosd-mREiJWw")
hdb_geocoded <- data.frame()
for (i in all_adds){
  append_data <- data.frame(
    add = i,
    geocode(i, source = "google")
  )
  hdb_geocoded <- rbind(hdb_geocoded, append_data)
}
# write to file
write.csv(hdb_geocoded,"GeoCoded.csv")
###


### Align Spatial Point CRS with Shapefile CRS
# convert SpatialPointsDataFrame to SpatialPoints Object
hdb_sp = SpatialPoints(cbind(hdb_coord$lon, hdb_coord$lat), proj4string = CRS("+proj=longlat"))
# align CRS of SpatialPoints with shapefile by converting coordinate to UTM using CRS SVY21
hdb_sp_UTM <- spTransform(hdb_sp, CRS("SVY21"))
# plot SpatialPoints with polygons of subzones
tm_shape(mpsz) + tm_borders(alpha=0.5) + tm_shape(hdb_sp_UTM) + tm_dots()
###

#### Point Pattern Analysis using Quadrat Density
# convert SpatialPointsDataFrame to sf Object
hdb_point_pattern <- st_as_sf(hdb_coord, coords = c("lon", "lat"), crs = st_crs(mpsz))
# convert sf object as ppp object
pp <-as.ppp(st_geometry(hdb_point_pattern))
pp
plot(pp, main=NULL, cols=rgb(0,0,0,.2), pch=20)
Q <- quadratcount(pp, nx= 10, ny=8)
plot(pp,main = "", col=rgb(0,0,0,.2), pch=20)
plot(Q, add=TRUE)
Q.d <- intensity(Q)
plot(intensity(Q, image=TRUE), main=NULL, las=1)  
plot(pp, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)  

# rescale unit from kilometer to meter 
pp.m <- rescale(pp, 0.001, "m")
pp.m
Q <- quadratcount(pp.m, nx= 10, ny=8)
plot(pp.m,main = "", col=rgb(0,0,0,.2), pch=20)
plot(Q, add=TRUE)
Q.d <- intensity(Q)
plot(intensity(Q, image=TRUE), main=NULL, las=1)  
plot(pp.m, pch=20, cex=0.6, col=rgb(0,0,0,.5), add=TRUE)  

#### Point Pattern Analysis using KDE  (Non-parametric)
K1 <- density(pp.m) # Using the default bandwidth
plot(K1, main=NULL, las=1)
contour(K1, add=TRUE)

# with a 10m bandwidth
K2 <- density(pp.m, sigma=100) # Using a 100m bandwidth
plot(K2, main=NULL, las=1)
contour(K2, add=TRUE)
###

### Autocorelation
library(spdep)
library(tidyr)
mpsz_resale2021_new <- drop_na(mpsz_resale2021)
mpsz_resale2021_new <- mpsz_resale2021_new[,c('SUBZONE_N','Mean_Resale_Price')]
# define neighbors
nb <- poly2nb(mpsz_resale2021_new, queen=TRUE)
# test neighbors
nb[[1]]
mpsz_resale2021_new$SUBZONE_N[1]
mpsz_resale2021_new$SUBZONE_N[c(2,4,8,12,17,19)]
# 
lw <- nb2listw(nb, style="W", zero.policy=TRUE)
# calculate lagged resale price
Prce.lag <- lag.listw(lw, mpsz_resale2021_new$Mean_Resale_Price)
# create a regression model
M <- lm(Prce.lag ~ mpsz_resale2021_new$Mean_Resale_Price)
# plot the data
plot(Prce.lag ~ mpsz_resale2021_new$Mean_Resale_Price, pch=20, asp=1, las=1)
abline(M)
coef(M)[2]
moran.test(mpsz_resale2021_new$Mean_Resale_Price,lw)

# Monte Carlo test
MC<- moran.mc(mpsz_resale2021_new$Mean_Resale_Price, lw, nsim=999)
MC
plot(MC, main="Monte-Carlo test", las=1)

