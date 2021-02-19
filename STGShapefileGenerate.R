library(raster)
library(sf)
setwd("H:/Research/MW_ClimateChange")


PML_STGA_N.ext <- as(extent(c(-135, -55, 70, 89)), 'SpatialPolygons')
crs(PML_STGA_N.ext) <- "+init=epsg:4326"
shapefile(PML_STGA_N.ext, 'PML_STGA_N.shp')


PML_STGA_S.ext <- as(extent(c(-135, -55, 30, 50)), 'SpatialPolygons')
crs(PML_STGA_S.ext) <- "+init=epsg:4326"
shapefile(PML_STGA_S.ext, 'PML_STGA_S.shp')


WA_LSTCA_W.ext <- as(extent(c( -90, -80, 30, 40)), 'SpatialPolygons')
crs(WA_LSTCA_W.ext) <- "+init=epsg:4326"
shapefile(WA_LSTCA_W.ext, 'WA_LSTCA_W.shp')


WA_LSTCA_E.ext <- as(extent(c( -75, -65, 30, 40)), 'SpatialPolygons')
crs(WA_LSTCA_E.ext) <- "+init=epsg:4326"
shapefile(WA_LSTCA_E.ext, 'WA_LSTCA_E.shp')

EA_LSTCA_W.ext <- as(extent(c(-15, 0, 7.5, 32)), 'SpatialPolygons')
crs(EA_LSTCA_W.ext) <- "+init=epsg:4326"
shapefile(EA_LSTCA_W.ext, 'EA_LSTCA_W.shp')

EA_LSTCA_E.ext <- as(extent(c(-30, -15, 7.5, 32)), 'SpatialPolygons')
crs(EA_LSTCA_E.ext) <- "+init=epsg:4326"
shapefile(EA_LSTCA_E.ext, 'EA_LSTCA_E.shp')

EP_LSTCA_E.ext<- as(extent(-110,-105,25,32), 'SpatialPolygons')
crs(EP_LSTCA_E.ext) <- "+init=epsg:4326"
shapefile(EP_LSTCA_E.ext, 'EP_LSTCA_E.shp')

EP_LSTCA_W.ext<- as(extent(-120,-115,25,32), 'SpatialPolygons')
crs(EP_LSTCA_W.ext) <- "+init=epsg:4326"
shapefile(EP_LSTCA_W.ext, 'EP_LSTCA_W.shp')

