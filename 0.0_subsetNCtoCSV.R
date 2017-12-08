## To extract a subset of data from the netcdf files and save as csv
## Individual heat flux tyoes are stored in annual netcdf files

library(maptools)
library(ncdf4)
library(raster)
library(rasterVis)
library(rgeos)
library(sp)
library(lattice)
library(grid)
library(stringr)
library(doParallel)
library(data.table)

maxLat <- -15
minLat <- -50
maxLon <- 60
minLon <- -10

setwd("~/R/myFunctions/")

# setup bathymetry map -----------------------------------------------------
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(10,-40,40,-20),nrow=2,ncol=2)
x <- extent(locbb)
bathy.cont <- crop(bathy,x)
bathy.loc <- crop(bathy,x)      # crop bathy map
bathy.land <- crop(bathy,x)
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
rm(bathy,locbb,x)

#mapTheme <- rasterTheme(region = rev(brewer.pal(8, "YlOrRd")))
myColor0 <- as.list(colorRampPalette(c("red","yellow","white","springgreen","royalblue"))(1e3)) # neg-> 0-> pos
myColor1 <- as.list(colorRampPalette(c("red","yellow","springgreen","royalblue"))(1e3)) # 0-> pos
myColor2 <- as.list(colorRampPalette(c("royalblue","springgreen","yellow","red"))(1e3)) # temp 0 <- pos
myColor3 <- as.list(colorRampPalette(c("white","springgreen","royalblue"))(1e3))  # wind 0 <- pos
leg.title <- gpar(fontsize=10,cex=1,lineheight=0.4,fontface="bold")    # legend title settings

# get SA borders
borders <- readShapePoly("SAfrica-borders.shp", proj4string = crs)
borders <- SpatialPolygons(borders@polygons)
minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40
mapThemeBathy <- colorRampPalette(c("blue","white"),alpha=T)

# get station locations
station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata")
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

# get the file directories --------------------------------------------------
baseURL <- "/media/robert/KELP-HDD-Portable"
type <- "/OAFlux/data"
fluxes <- list(flux.latent=list("lh.*\\.nc$","lhtfl","latentHt","Latent Heat",
                                quote(W.m^{-2}),-200,600,myColor0,T), # wild card .* and escape the point \\ and $ indicates end of the string
               flux.humid=list("qa.*\\.nc$","hum2m","specificHum","Specific Humidity at 2 m",
                               quote(g.Kg^{-1}),0,25,myColor1,F),
               flux.specific=list("sh.*\\.nc$","shtfl","sensibleHt","Sensible Heat",
                                  quote(W.m^{-2}),-150,300,myColor0,T),
               flux.tatmos=list("ta.*\\.nc$","tmp2m","tempAtmos","Atmospheric Temp. at 2 m",
                                quote(~degree~C),-5,35,myColor2,F),
               flux.tsea=list("ts.*\\.nc$","tmpsf","tempSea","Ocean Temp.",
                              quote(~degree~C),0,35,myColor2,F),
               flux.wind=list("ws.*\\.nc$","wnd10","wind10m","Wind Spped at 10 m",
                              quote(m.s^{-1}),0,25,myColor3,F))

# for each heat flux type ------------------------------------------------------
# # get the min and max for each flux
# getLim <- function(filename){
#   nc_stor <- nc_open(filename)
#   data.flux <- ncvar_get(nc_stor, product.short) # 1 (360 x 180) deg by 365 day
#   ret <- list(c(min(data.flux,na.rm=T),max(data.flux,na.rm=T)))
#   rm(data.flux)
#   nc_close(nc_stor)
#   return(ret)
# }
# 
# minmax <- list()
# for (i in 1:6){
#   product <- fluxes[[i]][[1]]           # for search pattern
#   product.short <- fluxes[[i]][[2]]     # for nc variable
#   listFiles <- list.files(paste0(baseURL,type),pattern=product, full.names = TRUE)
#   output <- lapply(listFiles, function(r) getLim(r))
#  minmax[[i]] <- output
#  rm(output)
# }

for (i in 1:6){
  # list all files in the directory
  product <- fluxes[[i]][[1]]           # for search pattern
  product.short <- fluxes[[i]][[2]]     # for nc variable
  product.read <- fluxes[[i]][[3]]      # for file name
  product.label <- fluxes[[i]][[4]]     # for graph label
  product.unit <- fluxes[[i]][[5]]      # units for graph plot
  product.min <- fluxes[[i]][[6]]
  product.max <- fluxes[[i]][[7]]
  myColor <- unlist(fluxes[[i]][[8]])
  negpos <- fluxes[[i]][[9]]
  listFiles <- list.files(paste0(baseURL,type),pattern=product, full.names = TRUE)
  
  cat("Working with ", product.label, " data\n")
  for (j in 1:length(listFiles)){                  # for each year
  #for (j in 1:1){
    
    # for each annual file in the list
    nc_stor <- nc_open(listFiles[[j]])
    
    data.flux <- ncvar_get(nc_stor, product.short) # 1 (360 x 180) deg by 365 days
    
    lon <- nc_stor$dim$lon$vals
    ind <- which(lon>180 & lon<360)
    lon[ind] <- lon[ind]-360
    
    lat <- nc_stor$dim$lat$vals
    time <- nc_stor$dim$time$vals
    time_length <- nc_stor$dim$time$len
    
    lonlat <- expand.grid(lon, lat)
    
    yr <- str_sub(listFiles[[j]],-7,-4)       # get the year from the file name
    time_labels <- as.Date(strptime(paste(yr, time), format="%Y %j"))
    
    nc_close(nc_stor)
    rm(nc_stor)
    gc(verbose = FALSE)
    
    ########################################### Extract subset and save csv and png in parallel
    sourceURL <- ("/media/robert/KELP-HDD-Portable/OAFlux/")
    saveImage <- "images/"
    saveData <- "csv/"
    sourceTitle <- "OAFlux"
    
    # for time test
    strt<-Sys.time()
    
    #cl <- makeCluster(2)
    registerDoParallel(cores = 4)
    cat("Starting parallel processing\n")
    
    data.list <- foreach(k = 1:time_length,.inorder=TRUE,
    #data.list <- foreach(k = 1:10,.inorder=TRUE,                    
                         .packages=c("raster","rasterVis","grid","rgeos","sp","zoo","chron","lattice")) %dopar% {
                           
                           date_label <- time_labels[k]
                           date_file <- gsub("-", "", strftime(date_label, format = "%Y%m%d"))
                           
                           # create csv filename
                           filename1 <- paste0(sourceURL,saveImage,product.read,"/",sourceTitle,"_",date_file,"_",product.read,".png")
                           filename2 <-paste0(sourceURL,saveData,product.read,"/",sourceTitle,"_",date_file,"_",product.read,".csv")
                           
                           flux <- data.flux[,,k]
                           df <- data.frame(cbind(lonlat, as.vector(flux)))
                           rm(flux)
                           names(df) <- c("lon", "lat", "flux")
                           sub.df <- subset(df, lat <= maxLat & lat >= minLat & lon <= maxLon & lon >= minLon)
                           rm(df)
                           
                           # replace NA values with 0
                           #sub.df[is.na(sub.df$flux), 3] <- 99
                           
                           r  <- rasterFromXYZ(sub.df, crs = proj, digits = 3)
                           
                           if (negpos){
                             obj <- levelplot(r, col.regions = myColor, margin = F,
                                              xlab="Latitude",ylab="Longitude",
                                              at= unique(c(seq(product.min, 0, length=50), seq(0, product.max, length=50))), 
                                              main = paste("OAflux ",product.label,"\n",date_label, sep = '')) +
                               layer(sp.lines(borders, col = "black", lwd = 1))
                           } else {
                             obj <- levelplot(r, col.regions = myColor, margin = F,
                                              xlab="Latitude",ylab="Longitude",
                                              at= unique(c(seq(product.min, product.max, length=100))), 
                                              main = paste("OAflux ",product.label,"\n",date_label, sep = '')) +
                               layer(sp.lines(borders, col = "black", lwd = 1))
                           }
                           
                           gc(verbose = FALSE)
                           
                           # return data frame and image with filenames
                           list(sub.df, filename2, obj, filename1)
                         } # end foreach
    #stopCluster(cl)
    rm(data.flux)
    gc(verbose = FALSE)
    
    cat("Beginning writing files\n")
    for (l in 1:length(data.list)){
      
      # save data frame as CSV file
      fwrite(data.list[[l]][[1]], file.path = data.list[[l]][[2]],col.names = T)
      gc(verbose = FALSE)
      
      # save image as png
      png(data.list[[l]][[4]],width = 9,height = 5,units = 'in',res = 300)
      print(data.list[[l]][[3]])
      trellis.focus("legend",side ="right",clip.off = TRUE,highlight = FALSE)
      grid.text(product.unit, 0.3, 1.08, rot = 0, gp=leg.title)
      trellis.unfocus()
      dev.off()
      
    } # for l
    rm(data.list)
    gc(verbose = FALSE)
  }
}

