srtm.region <- function(w,e,s,n,path=getwd(),del.tmp=TRUE,plot=FALSE,out="output",type=""){

if(sum(dir(path)=="tmp")==0){dir.create(paste(path,"/tmp",sep=""))}
  
srtm_x_min <- ceiling( (w + 180.00001) / 5 )
srtm_x_max <- ceiling( (e + 179.99999) / 5 )

srtm_y_min <- ceiling( (59.99999 - s) / 5 )
srtm_y_max <- ceiling( (60.00001 - n) / 5 )

for(x in srtm_x_min:srtm_x_max){
  for(y in srtm_y_min:srtm_y_max){
    if(x<10 && y<10){ name <- paste("srtm_0",x,"_0",y,sep="") }
    if(x<10 && y>9) { name <- paste("srtm_0",x,"_",y,sep="") }
    if(x>9 && y<10) { name <- paste("srtm_",x,"_0",y,sep="") }
    if(x>9 && y>9)  { name <- paste("srtm_",x,"_",y,sep="") }
    print(name)
    
    if(sum(name==srtm.list())==1){
      if(sum(dir(paste(path,"/tmp/",sep=""))==paste(name,".tif",sep=""))==0){
	unzip(files=paste(name,".tif",sep=""), zipfile = paste(path,"/",name,".zip",sep=""), exdir = "tmp")
      }
    
      #extract meta information
#      meta <- GDALinfo(paste(path,"/tmp/",name,".tif",sep=""),silent=TRUE)
#   
#     rows <- as.numeric(meta[1])
#     cols <- as.numeric(meta[2])
#     x_ori <- as.numeric(meta[4])
#     y_ori <- as.numeric(meta[5])
#     x_res <- as.numeric(meta[6])
#     y_res <- as.numeric(meta[7])
  
############ use celing and floor ???? Maybe helps.    

#      x_off <- (w-x_ori)/x_res
#      if(x_off<1){x_off<-0}
#      y_off <- rows-((n-y_ori)/y_res)
#      if(y_off<1){y_off<-0}
#    
#      x_reg <- ((e-x_ori)/x_res)-x_off+1
#      if(x_reg>(cols-x_off)){x_reg<-cols-x_off}
#      y_reg <- rows-((s-y_ori)/y_res)-y_off+1
#      if(y_reg>(rows-y_off)){y_reg<-rows-y_off}

    
#      srtm_read <- raster(readGDAL(paste(path,"/tmp/",name,".tif",sep=""),offset=c(y_off,x_off),region.dim=c(y_reg,x_reg)))

      srtm_read <- raster(paste(path,"/tmp/",name,".tif",sep=""))
      
      if(sum(ls()=="srtm")==0){
	srtm <- crop(srtm_read,extent(w,e,s,n))}else{
	  srtm_crop <- crop(srtm_read,extent(w,e,s,n));
	  srtm <- merge(srtm,srtm_crop)
	  #srtm <- expand(srtm,srtm_crop,value=srtm_crop@data@values)
	  #srtm <- expand(srtm,srtm_crop,value=c(srtm_crop@data@values,srtm@data@values),filename="temp",overwrite=TRUE)
	  }
      }else{print("At least a part of your requested area is outside the valid SRTM. Only ocean?")}
  }
}

if(del.tmp==TRUE){unlink(paste(path,"/tmp",sep=""),recursive=TRUE)}

if(sum(ls()=="srtm")==1){
  if(plot==TRUE){
    X11()
    plot(srtm, col=rev(terrain.colors(25)), maxpixels=100000, axes = TRUE, xlab="",ylab="", alpha=1, useRaster=TRUE)
    }
  if(type=="tif"){
    out <- writeRaster(srtm, out, "GTiff", options=c("COMPRESS=NONE", "TFw=YES"), overwrite=TRUE)}
  if(type=="shp"){
    polygon <- rasterToPolygons(srtm, dissolve=TRUE)
    writeOGR(polygon, dsn=".", layer="polygon", driver="ESRI Shapefile")}
  return(srtm)
  }
}