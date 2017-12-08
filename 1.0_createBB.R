# To create bounding boxes for extracting OA fluxes
# Six BB protracting from the station coords seaward

# From the station create an box going throught the coords of the station such that the
#   station is in the centre of the line
# Determine the number of land pixels within the box

library(raster)
library(rasterVis)
library(rgeos)
library(maptools)

maxLat <- -15
minLat <- -50
maxLon <- 60
minLon <- -10

setwd("~/R/myFunctions/")

# setup bathymetry map -----------------------------------------------------
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(minLon,minLat,maxLon,maxLat),nrow=2,ncol=2)
x <- extent(locbb)
bathy.land <- crop(bathy,x)
bathy.nan <- bathy.land
bathy.nan[bathy.nan>=0] <- NA
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

mapThemeBathy <- colorRampPalette(c("blue","white","brown"),alpha=F)
at <- c(10000,0,-50,-100,-200,-500,-1000,-2000,-5000,-10000)
rm(bathy,locbb,x)

# get station locations
station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata")
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

for (i in 1:nrow(station.loc)){
  
  stat.name <- station.loc$station[i]
  cat(paste0("Starting bbox for ",stat.name,'\n'))
  
  st <- station.pts@data[i,1:2]
  coordinates(st) <- ~lon+lat
  
  lon.st <- st@coords[1]
  lat.st <- st@coords[2]

  # add degrees to the north and south of the lat and save points
  NE.15st <- c(lon.st,lat.st+7.5)   # 15 deg block
  NE.7st <- c(lon.st,lat.st+3.5)    # 7 deg block
  NE.2st <- c(lon.st,lat.st+1)      # 2 deg block
  NE.1st <- c(lon.st,lat.st+0.5)    # 1 deg block
  NE.05st <- c(lon.st,lat.st+0.25)   # 0.5 deg block
  NE.025st <- c(lon.st,lat.st+0.125)   # 0.25 deg block
  
  SE.15st <- c(lon.st,lat.st-7.55)
  SE.7st <- c(lon.st,lat.st-3.5)
  SE.2st <- c(lon.st,lat.st-1)
  SE.1st <- c(lon.st,lat.st-0.5)    # 1 deg block
  SE.05st <- c(lon.st,lat.st-0.25)   # 0.5 deg block
  SE.025st <- c(lon.st,lat.st-0.125)   # 0.25 deg block
  
  # subtract twice the degrees from each of the lon coords to create remaining two corners
  NW.15st <- c(NE.15st[1]-15,NE.15st[2])
  NW.7st <- c(NE.7st[1]-7,NE.7st[2])
  NW.2st <- c(NE.2st[1]-2,NE.2st[2])
  NW.1st <- c(NE.1st[1]-1,NE.1st[2])
  NW.05st <- c(NE.05st[1]-0.5,NE.05st[2])
  NW.025st <- c(NE.025st[1]-0.25,NE.025st[2])
  
  SW.15st <- c(NE.15st[1]-15,SE.15st[2])
  SW.7st <- c(NE.7st[1]-7,SE.7st[2])
  SW.2st <- c(NE.2st[1]-2,SE.2st[2])
  SW.1st <- c(NE.1st[1]-1,SE.1st[2])
  SW.05st <- c(NE.05st[1]-0.5,SE.05st[2])
  SW.025st <- c(NE.025st[1]-0.25,SE.025st[2])
  
  # repeat the first value
  bb.15 <- as.matrix(rbind(NE.15st,SE.15st,SW.15st,NW.15st,NE.15st))
  bb.7 <- as.matrix(rbind(NE.7st,SE.7st,SW.7st,NW.7st,NE.7st))
  bb.2 <- as.matrix(rbind(NE.2st,SE.2st,SW.2st,NW.2st,NE.2st))
  bb.1 <- as.matrix(rbind(NE.1st,SE.1st,SW.1st,NW.1st,NE.1st))
  bb.05 <- as.matrix(rbind(NE.05st,SE.05st,SW.05st,NW.05st,NE.05st))
  bb.025 <- as.matrix(rbind(NE.025st,SE.025st,SW.025st,NW.025st,NE.025st))
  
  bb.15poly <- Polygon(bb.15)
  bb.7poly <- Polygon(bb.7)
  bb.2poly <- Polygon(bb.2)
  bb.1poly <- Polygon(bb.1)
  bb.05poly <- Polygon(bb.05)
  bb.025poly <- Polygon(bb.025)
  
  bb.15sp <-  SpatialPolygons(list(Polygons(list(bb.15poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bb.7sp <-  SpatialPolygons(list(Polygons(list(bb.7poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bb.2sp <-  SpatialPolygons(list(Polygons(list(bb.2poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bb.1sp <-  SpatialPolygons(list(Polygons(list(bb.1poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bb.05sp <-  SpatialPolygons(list(Polygons(list(bb.05poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bb.025sp <-  SpatialPolygons(list(Polygons(list(bb.025poly), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  # rotate the polygons according to longitude
  if (lon.st > 20 & lon.st <= 28){
    bb.15sp <- elide(bb.15sp,rotate=-90,center=c(lon.st,lat.st))
    bb.7sp <- elide(bb.7sp,rotate=-90,center=c(lon.st,lat.st))
    bb.2sp <- elide(bb.2sp,rotate=-90,center=c(lon.st,lat.st))
    bb.1sp <- elide(bb.1sp,rotate=-90,center=c(lon.st,lat.st))
    bb.05sp <- elide(bb.05sp,rotate=-90,center=c(lon.st,lat.st))
    bb.025sp <- elide(bb.025sp,rotate=-90,center=c(lon.st,lat.st))
    
  } else if (lon.st > 28) {
    bb.15sp <- elide(bb.15sp,rotate=-135,center=c(lon.st,lat.st))
    bb.7sp <- elide(bb.7sp,rotate=-135,center=c(lon.st,lat.st))
    bb.2sp <- elide(bb.2sp,rotate=-135,center=c(lon.st,lat.st))
    bb.1sp <- elide(bb.1sp,rotate=-135,center=c(lon.st,lat.st))
    bb.05sp <- elide(bb.05sp,rotate=-135,center=c(lon.st,lat.st))
    bb.025sp <- elide(bb.025sp,rotate=-135,center=c(lon.st,lat.st))
  }
  
  # print flux boundaring box
  filename1 <- paste0("~/R_projects-CS/ame-temporalbabe/ExtractOA/Output/",stat.name,"_",
                      "OAfluxBBox.png")
  
  png(filename1, width = 7, height = 5, units = 'in', res = 300)
  
  p <- contourplot(bathy.land, main=list(paste0(stat.name," Extracted Sub-regions"), x=.6, just="center"),
              region=T,col.regions=mapThemeBathy, alpha.regions=.3,
              at=at, margin = F, labels = F, ylab="Latitude", xlab="Longitude") +
    layer(sp.polygons(bb.15sp, col="red",lwd=2)) +
    layer(sp.polygons(bb.7sp, col="red",lwd=2)) +
    layer(sp.polygons(bb.2sp, col="red",lwd=2)) +
    layer(sp.polygons(bb.1sp, col="red2",lwd=1.8)) +
    layer(sp.polygons(bb.05sp, col="red4",lwd=1.5)) +
    layer(sp.polygons(bb.025sp, col="black",lwd=1)) +
    #layer(sp.points(st, col="blue", pch=19, cex=1.2, fill=F, lwd=2))   # solid dot
    layer(sp.points(st, col="blue", pch=1, cex=1.2, fill=F, lwd=2))
  p <- update(p, aspect="iso")
  print(p)
  dev.off()
  
  # save boundaring box spatial polygons
  filename2 <- paste0("~/R_projects-CS/ame-temporalbabe/ExtractOA/Output/",stat.name,"_",
                      "OAfluxPolygons.Rdata")

  save(bb.15sp, bb.7sp, bb.2sp, bb.1sp, bb.05sp, bb.025sp, file=filename2)
  rm(filename2,filename1)
  
  
  # extract data and count number of land values ---------------------------------------------
  # land.ext <- unlist(extract(bathy.nan,bb.10sp))
  # land.cnt <- sum(is.na(land.ext))
  # rm(land.ext)
  # 
  # # rotate polygon -45 deg and recount
  # new.cnt <- 0
  # rot.poly <- bb.10sp
  # while (new.cnt<land.cnt){
  #   rot.poly <- elide(rot.poly,rotate=-15,center=c(lon.st,lat.st))
  #   land.ext <- unlist(extract(bathy.nan,bb.10sp))
  #   new.cnt <- sum(is.na(land.ext))
  #     if (new.cnt < land.cnt){
  #       land.cnt <- new.cnt
  #       new.cnt <- 0
  #     }
  # }
  # bb.10sp <- rot.poly
  # ------------------------------------------------------------------------------------------
  
}


