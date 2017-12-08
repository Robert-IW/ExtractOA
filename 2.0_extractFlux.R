# To extract the median value of latent, sensible heat, humidity, wind strength,
#   lower atmospher temp & SST

library(raster)
library(rasterVis)
library(rgeos)
library(maptools)
library(grid)
library(data.table)
library(doParallel)

doGraphic <- F

setwd("~/R/myFunctions/")

maxLat <- -15
minLat <- -50
maxLon <- 60
minLon <- -10

# setup bathymetry map -----------------------------------------------------
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.4/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(minLon,minLat,maxLon,maxLat),nrow=2,ncol=2)
x <- extent(locbb)
bathy.land <- crop(bathy,x)
bathy.nan <- bathy.land
bathy.nan[bathy.nan>=0] <- NA
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

rm(bathy,locbb,x)

# get SA borders
borders <- readShapePoly("SAfrica-borders.shp", proj4string = crs)
borders <- SpatialPolygons(borders@polygons)
minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40
mapThemeBathy <- colorRampPalette(c("blue","white"),alpha=T)

#mapTheme <- rasterTheme(region = rev(brewer.pal(8, "YlOrRd")))
myColor0 <- as.list(colorRampPalette(c("red","darkorange1","yellow","white","palegreen","lightseagreen","royalblue"))(1e3)) # neg-> 0-> pos
myColor1 <- as.list(colorRampPalette(c("red","yellow","springgreen","royalblue"))(1e3)) # 0-> pos
myColor2 <- as.list(colorRampPalette(c("royalblue","springgreen","yellow","red"))(1e3)) # temp 0 <- pos
myColor3 <- as.list(colorRampPalette(c("white","springgreen","royalblue"))(1e3))  # wind 0 <- pos
leg.title <- gpar(fontsize=10,cex=1,lineheight=0.4,fontface="bold")    # legend title settings

# get station locations
station.loc <- read.table("~/R/myFunctions/station_locations.Rdata",stringsAsFactors = F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

fluxes <- list(flux.latent=list("latentHt","Latent Heat",
                                quote(W.m^{-2}),-200,600,myColor0,T), # wild card .* and escape the point \\ and $ indicates end of the string
               flux.humid=list("specificHum","Specific Humidity at 2 m",
                               quote(g.Kg^{-1}),0,25,myColor1,F),
               flux.specific=list("sensibleHt","Sensible Heat",
                                  quote(W.m^{-2}),-150,300,myColor0,T),
               flux.tatmos=list("tempAtmos","Atmospheric Temp. at 2 m",
                                quote(~degree~C),-5,35,myColor2,F),
               flux.tsea=list("tempSea","Ocean Temp.",
                              quote(~degree~C),0,35,myColor2,F),
               flux.wind=list("wind10m","Wind Spped at 10 m",
                              quote(m.s^{-1}),0,25,myColor3,F))

for (h in 1:nrow(station.loc)){
  
  stat.name <- station.loc$station[h]
  cat(paste0("Starting bbox for ",stat.name,'\n'))
  
  st <- station.pts@data[h,1:2]
  coordinates(st) <- ~lon+lat
  
  # load the spatial polygons for each station
  load(paste0("~/R_projects-CS/ame-temporalbabe/ExtractOA/Output/",stat.name,"_OAfluxPolygons.Rdata"))
  
  # for each heat flux variable
  for (i in 1:6){
    # list all files in the directory
    product <- fluxes[[i]][[1]]           # for search pattern
    product.label <- fluxes[[i]][[2]]     # for graph label
    product.unit <- fluxes[[i]][[3]]      # units for graph plot
    product.min <- fluxes[[i]][[4]]
    product.max <- fluxes[[i]][[5]]
    myColor <- unlist(fluxes[[i]][[6]])
    negpos <- fluxes[[i]][[7]]
    
    cat(paste0("Extracting ",product.label," data\n"))
    
    listFiles <- list.files(paste0("/media/robert/KELP-HDD-Portable//OAFlux/csv/",product,"/"), pattern=product,full.names = TRUE)
    
    #registerDoParallel(cl, cores = detectCores()-2)
    cl <- makeCluster(6)
    registerDoParallel(cl, cores = 6)
    
    date.names <- lapply(listFiles, function(x) strsplit(basename(x),"_")[[1]][[2]])
    #date.names <- date.names[1:5113]                                   # for test run
    
    cat("Starting parallel processing\n")
    data.list <- foreach(j = 1:length(listFiles),.inorder=TRUE,
    #data.list <- foreach(j = 1:5113,.inorder=TRUE,                     # to test on 10 samples
                         .final = function(x) setNames(x, date.names),
                         .packages=c("raster","rasterVis","data.table","grid","rgeos","sp","zoo","chron","lattice")) %dopar% {
                           
                           saveData <- "OAFlux/extract/latent/"
                           saveImage <- "OAFlux/extrImage/latent"
                           
                           # convert csv to raster
                           inter1= fread(listFiles[j], header=TRUE)
                           #convert to spatial points
                           coordinates(inter1) = ~lon + lat
                           #gridify your set of points
                           gridded(inter1) <- TRUE 
                           #convert to raster
                           r <- raster(inter1)
                           
                           median15 <- extract(r, bb.15sp, fun=median, na.rm=T)
                           median7 <- extract(r, bb.7sp, fun=median, na.rm=T)
                           median2 <- extract(r, bb.2sp, fun=median, na.rm=T)
                           
                           stdDev15 <- extract(r, bb.15sp, fun=sd, na.rm=T)
                           stdDev7 <- extract(r, bb.7sp, fun=sd, na.rm=T)
                           stdDev2 <- extract(r, bb.2sp, fun=sd, na.rm=T)
                           
                           list(c("med15deg"=median15,"std15deg"=stdDev15),
                                c("med7deg"=median7,"std7deg"=stdDev7),
                                c("med2deg"=median2,"std2deg"=stdDev2),"raster"=r)
                         }
    
    stopCluster(cl)
    cat("Processing complete\n")
    
    df.pre <- lapply(data.list, function(x) matrix(unlist(x[1:3]),ncol=6))
    
    dt <- data.table(row.names = names(df.pre),station.loc[h,1],station.loc[h,2],matrix(unlist(df.pre), ncol=6, byrow=T),
                     stringsAsFactors=FALSE)
    names(dt) <- c("date","lon","lat","median15deg","std15deg","median7deg","std7deg","median2deg","std2deg")
    
    # save list of heat flux for all days at station
    filename1 <- paste0("/media/robert/KELP-HDD-Portable/OAFlux/extractCSV/",product,"/",
                        stat.name,"_",product,".csv")
    #save(data.list, file=filename1)
    fwrite(dt, file.path = filename1, col.names=T)
    rm(df.pre, dt,filename1)
    gc()
  
    
    if (doGraphic){
      cat("Starting printing of images\n")
      
      # extropolate colours from legend that match the median values
      if (negpos){
        limits1 <- round(unique(c(seq(product.min, 0, length=50), seq(0, product.max, length=50))))
        limits2 <- round(unique(c(seq(product.min, 0, length=500), seq(0, product.max, length=500))))
      } else {
        limits1 <- round(unique(c(seq(product.min, product.max, length=100))))
        limits2 <- round(unique(c(seq(product.min, product.max, length=1000))))
      }
      
      # rounding means that some integers may be omitted so nearest match is used
      temp <- sapply(limits1, function(x) which(abs(limits2-x)==min(abs(limits2-x)))[1])
      myColor.sub <- myColor[temp]
      for (k in 1:5){
        
        data.plot <- data.list[[k]]
        
        median15 <- data.plot[[1]][[1]]
        median7 <- data.plot[[2]][[1]]
        median2 <- data.plot[[3]][[1]]
        r <- data.plot[[4]]
        
        col.median15 <- myColor.sub[which(abs(limits1-median15)==min(abs(limits1-median15)))]
        col.median7 <- myColor.sub[which(abs(limits1-median7)==min(abs(limits1-median7)))]
        col.median2 <- myColor.sub[which(abs(limits1-median2)==min(abs(limits1-median2)))]
        
        #rm(median15,median7,median2)
        
        date.name <- names(data.list[k])
        
        filename2 <- paste0("/media/robert/KELP-HDD-Portable/OAFlux/extractPNG/",product,"/",
                            stat.name,"_",product,"_",date.name,".png")
        
        date.name <- as.Date(sprintf("%08s", date.name), format="%Y%m%d")
        
        png(filename2, width = 7, height = 5, units = 'in', res = 300)
        
        if (negpos){
          obj <- levelplot(r, col.regions = myColor, margin = F,
                           xlab="Latitude",ylab="Longitude",
                           at= unique(c(seq(product.min, 0, length=50), seq(0, product.max, length=50))), 
                           main = paste(stat.name,"\nOAflux ",product.label,"\n",date.name, sep = '')) +
            layer(sp.lines(borders, col = "black", lwd = 1)) +
            layer(sp.polygons(bb.15sp, col="red",lwd=3,fill=col.median15)) +
            layer(sp.polygons(bb.7sp, col="red",lwd=3,fill=col.median7)) +
            layer(sp.polygons(bb.2sp, col="red",lwd=3,fill=col.median2)) +
            layer(sp.points(st, col="blue", pch=19, cex=1.5, fill=F, lwd=2))
          
          obj <- update(obj, aspect="iso")
          print(obj)
          trellis.focus("legend",side ="right",clip.off = TRUE,highlight = FALSE)
          grid.text(product.unit, 0.3, 1.08, rot = 0, gp=leg.title)
          trellis.unfocus()
          
          dev.off()
          
        } else {
          obj <- levelplot(r, col.regions = myColor, margin = F,
                           xlab="Latitude",ylab="Longitude",
                           at= unique(c(seq(product.min, product.max, length=100))), 
                           main = paste(stat.name,"\nOAflux ",product.label,"\n",date.name, sep = '')) +
            layer(sp.lines(borders, col = "black", lwd = 1)) +
            layer(sp.polygons(bb.15sp, col="red",lwd=3,fill=col.median15)) +
            layer(sp.polygons(bb.7sp, col="red",lwd=3,fill=col.median7)) +
            layer(sp.polygons(bb.2sp, col="red",lwd=3,fill=col.median2)) +
            layer(sp.points(st, col="blue", pch=19, cex=1.5, fill=F, lwd=2))
          
          obj <- update(obj, aspect="iso")
          print(obj)
          trellis.focus("legend",side ="right",clip.off = TRUE,highlight = FALSE)
          grid.text(product.unit, 0.3, 1.08, rot = 0, gp=leg.title)
          trellis.unfocus()
          
          dev.off()
        } # if negpos
      } # for k in 
    } # if doGraphic
    rm(data.list)
    gc()
  } # for i in 
} # for h in