srtm.slr <- function(srtm,slr,pos=FALSE,type="",out="output",cca.s=2,mode=3,plot=FALSE){

srtm[ is.na(srtm) ] <- -9999

if(pos[1]==FALSE){
  nas <- which(srtm@data@values==-9999)
  nas <- nas[round(length(nas)/2)] 
  y <- ceiling(nas/srtm@ncols)
  x <- nas-((y-1)*srtm@ncols)}else{
    #no resolution??
    x <- (pos[1]-srtm@extent@xmin)/((srtm@extent@xmax-srtm@extent@xmin)/srtm@ncols)
    y <- (srtm@extent@ymax-pos[2])/((srtm@extent@ymax-srtm@extent@ymin)/srtm@nrows)
    if(srtm[y*srtm@ncols+x]!=-9999){print("!!! YOUR POSITION CONTAINS NO OCEAN and will be replaced!!!")
                              nas <- which(srtm@data@values==-9999)
                              nas <- nas[round(length(nas)/2)] 
                              y <- ceiling(nas/srtm@ncols)
                              x <- nas-((y-1)*srtm@ncols)
                              }
  }

# reclassify the values into two groups
m <- c(-Inf, slr, 1, slr, Inf, 0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(srtm, rclmat)

# slower as cca !!!
#clump <- clump(rc, filename="", directions=8, gaps=TRUE)

#library(cca)
matrix <- matrix(rc@data@values, ncol = rc@ncols, byrow = TRUE)
matrix_cluster <- raster(cca.single(matrix, s=cca.s, x=y,y=x, mode = mode))
matrix_cluster@data@values[ which(matrix_cluster@data@values==0)] <- NA

rc@data@values <- matrix_cluster@data@values
srtm@data@values[ which(srtm@data@values==-9999)] <- NA

final <- overlay(srtm, rc, fun=function(x,y){(x*y)})

if(plot==TRUE){
  X11()
  plot(final, col=rev(terrain.colors(25)), maxpixels=100000, axes = TRUE, xlab="",ylab="", alpha=1, useRaster=TRUE)
}
if(type=="tif"){
  out <- writeRaster(final, out, "GTiff", options=c("COMPRESS=NONE", "TFw=YES"), overwrite=TRUE)
  }

if(type=="shp"){
  m <- c(-Inf, Inf, 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  final_s <- reclassify(final, rclmat)
  polygon <- rasterToPolygons(final_s, dissolve=TRUE)
  #plot(polygon, add=T, col="red")
  writeOGR(polygon, dsn=".", layer=out, driver = "ESRI Shapefile")
  }
return(final)
#simple <- gSimplify(polygon, tol=10)
#plot(simple)
#drv <- "ESRI Shapefile"
#writeOGR(simple, dsn=".", layer="simple", driver =drv)
}