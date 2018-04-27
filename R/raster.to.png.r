raster.to.png <- function(raster, name="test", mypal=c("wheat1", "red3","green"),n=10, style="equal", fixed=NULL,trans=1,legend=FALSE){
  png.path <- getOption("png.path",default=getwd())
  
  matrix <- matrix(raster[],nrow=raster@nrows,ncol=raster@ncols,byrow=TRUE)
  h5 <- classIntervals(as.numeric(matrix[]), n=n, style=style, fixedBreaks=fixed)
  h5Colours <- findColours(h5, mypal)
  
  if(legend==TRUE){
    pdf(paste(png.path,"/",name,"_legend.pdf",sep=""))
    plot(1,1,col="white",axes=FALSE, xlab="",ylab="")
    legend("top", fill=attr(h5Colours, "palette"), legend=names(attr(h5Colours, "table")), bg="white")
    dev.off()
  }
  
  rgba <- array(dim=c(dim(matrix),4))
  # transperancy
  rgba[,,4] <- trans
  rgba[,,4][is.na(matrix)] <- 0
  # colours
  rgb.data <- col2rgb(h5Colours)
  rgba[,,1] <- rgb.data[1,]/255
  rgba[,,2] <- rgb.data[2,]/255
  rgba[,,3] <- rgb.data[3,]/255
  
  writePNG(rgba, target=paste(png.path,"/",name,".png",sep=""))   
}
