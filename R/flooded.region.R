      
flooded.region <-function(x,y,w,e,s,n,write_path){
      print("... load flooded region ...")
      meta <- GDALinfo(paste(write_path,x,"_",y,".tif",sep=""),silent=TRUE)
    
      rows <- as.numeric(meta[1])
      cols <- as.numeric(meta[2])
      x_ori <- as.numeric(meta[4])
      y_ori <- as.numeric(meta[5])
      x_res <- as.numeric(meta[6])
      y_res <- as.numeric(meta[7])
  
    
      x_off <- (w-x_ori)/x_res
      if(x_off<1){x_off<-0}
      y_off <- rows-((n-y_ori)/y_res)
      if(y_off<1){y_off<-0}
    
      x_reg <- ((e-x_ori)/x_res)-x_off
      if(x_reg>(cols-x_off)){x_reg<-cols-x_off}
      y_reg <- rows-((s-y_ori)/y_res)-y_off
      if(y_reg>(rows-y_off)){y_reg<-rows-y_off}
      
      #print(w)
      #print(e)
      #print(s)
      #print(n)
      
      #print(x_off)
      #print(y_off)
      #print(x_reg)
      #print(y_reg)
      neighbour_border <- raster(readGDAL(paste(write_path,x,"_",y,".tif",sep=""),offset=c(y_off,x_off),region.dim=c(y_reg,x_reg)))
      
      #print("... readed...")
      return(neighbour_border)
}