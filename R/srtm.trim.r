
library(raster)

srtm.trim <- function(raster)
{	
  na <- -999
  raster[is.na(raster[])] <- na
  coord <- rep(0,4)
  nrc <- c(raster@nrows,raster@ncols)
  print("NA replacement done")
	fromC <- .C("srtmtrim",
           raster=as.double(raster[]),
	   #test=as.integer(list(raster@nrows,raster@ncols)),
           #nrow=as.integer(raster@nrows),
           #ncol=as.integer(raster@ncols),
           #n=as.integer(0),s=as.integer(0),
           #w=as.integer(0),e=as.integer(0)
	   nrc=as.integer(nrc),
	   coord=as.integer(coord)
	   )

  #print(fromC)

  extent <- extent(xFromCol(raster,fromC$coord[3]+1)-0.5*xres(raster),
		   xFromCol(raster,fromC$coord[4]+1)+0.5*xres(raster),
		   yFromRow(raster,fromC$coord[2]+1)-0.5*yres(raster),
		   yFromRow(raster,fromC$coord[1]+1)+0.5*yres(raster)
		  )

  #print(extent)
  
  final <- raster(ncols=fromC$nrc[2], nrows=fromC$nrc[1], extent)
  final[] <- fromC$raster[1:(fromC$nrc[1] * fromC$nrc[2])]
  final[which(final[]==na)] <- NA
  return(final);
}

