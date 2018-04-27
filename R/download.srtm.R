download.srtm <- function(url=3, all=TRUE, w, e, s, n, dest=getwd()){
  
  if(url==1){url <- "ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/"}
  if(url==2){url <- "http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/"}
  if(url==3){url <- "http://droppr.org/srtm/v4.1/6_5x5_TIFs/"}
  if(url==4){url <- "ftp://xftp.jrc.it/pub/srtmV4/tiff/"}
  
  if(all==TRUE){
    srtm_x_min <- 1
    srtm_x_max <- 72
    srtm_y_min <- 1
    srtm_y_max <- 24
  }
  else{
    srtm_x_min <- ceiling( (w + 180.00001) / 5 )
    srtm_x_max <- ceiling( (e + 179.99999) / 5 )

    srtm_y_min <- ceiling( (59.99999 - s) / 5 )
    srtm_y_max <- ceiling( (60.00001 - n) / 5 )
  }

  for(x in srtm_x_min:srtm_x_max){
    for(y in srtm_y_min:srtm_y_max){
      if(x<10 && y<10){ name <- paste("srtm_0",x,"_0",y,sep="") }
      if(x<10 && y>9) { name <- paste("srtm_0",x,"_",y,sep="") }
      if(x>9 && y<10) { name <- paste("srtm_",x,"_0",y,sep="") }
      if(x>9 && y>9)  { name <- paste("srtm_",x,"_",y,sep="") }
      print(name)
      url_c <- paste(url,name,".zip",sep="")
      try(download.file(url_c, destfile=paste(dest,"/",name,".zip",sep="")))
    }
  }
}