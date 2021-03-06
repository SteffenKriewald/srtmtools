compare.srtm.borders <- function(region,x,y,slr,write_path){
 
### north border
   if(sum(dir(write_path)==paste(x,"_",y+1,".tif",sep=""))==1){
      print(" ... load flooded north border ...")
      
      n_index <- NULL

      neighbour <- raster(paste(writepath,"/tmp/",name,".tif",sep=""))
      n_neighbour <- getValues(neighbour, row=6000)
      rm(neighbour)
      n_border <- getValues(region, row=1)
      
      for(i in 1:6000){
	if(is.na(n_neighbour[i])==FALSE){
	  for(j in -1:1){
	    if(is.na(n_border[i+j])==FALSE){
	      if(n_border[i+j]<=slr){n_index <- c(n_index,i+j)}
	    }
	  }
	}
      }
    if(is.null(n_index)==FALSE){n_values <- n_border[n_index]; n_border[n_index] <- NA}
    }
### west border    
    if(sum(dir(write_path)==paste(x-1,"_",y,".tif",sep=""))==1){
      print(" ... load flooded west border ...")

      w_index <- NULL

      neighbour <- raster(paste(writepath,"/tmp/",name,".tif",sep=""))
      w_neighbour <- getValuesBlock(neighbour, row=1, nrow=6000, col=6000, ncol=1)
      rm(neighbour)
      w_border <- getValuesBlock(region, row=1, nrow=6000, col=1, ncol=1)
     
      for(i in 1:6000){
	if(is.na(w_neighbour[i])==FALSE){
	  for(j in -1:1){
	    if(is.na(w_border[i+j])==FALSE){
	      if(w_border[i+j]<=slr){w_index <- c(w_index,i+j)}
	    }
	  }
	}
      }
      if(is.null(w_index)==FALSE){w_values <- w_border[w_index]; w_border[w_index] <- NA}
   }
    
    


    return(list(region,neighbour_border))
}