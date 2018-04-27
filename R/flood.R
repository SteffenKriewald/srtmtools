flood <- function(srtm,x,y,slr,write_path,write=TRUE){
    ocean_srtm <- srtm
    # reclass
    nas <- which(is.na(ocean_srtm[])==TRUE)
    ocean_srtm[nas] <- 1
    ocean_srtm@data@values[-nas] <- 0

    rm(nas)

    # clump
    matrix <- matrix(ocean_srtm[], ncol=ocean_srtm@ncols, byrow=TRUE)
    matrix <- raster(cca(matrix))
    ocean_srtm[] <- matrix[]
    #ocean_srtm <- clump(ocean_srtm)
    max <- max(na.omit(ocean_srtm[]))

    rm(matrix)

    ocean_pos <- 1:max

    for(i in 1:max){
      position <- which(ocean_srtm[]==i)
      ocean_pos[i]<-position[1]
    }

    rm(ocean_srtm)
    
    ############

    srtm_slr <- srtm
    srtm_slr[ is.na(srtm_slr[]) ] <- -9999

    # reclass
    srtm_slr[which(srtm_slr[]<slr)] <- 1
    srtm_slr[which(srtm_slr[]>=slr)] <- 0

    #m <- c(-Inf, slr, 1, slr, Inf, 0)
    #rclmat <- matrix(m, ncol=3, byrow=TRUE)
    #srtm_slr <- reclass(srtm_slr, rclmat)

    # clump
    matrix <- matrix(srtm_slr[], ncol=srtm_slr@ncols,byrow=TRUE)
    matrix <- raster(cca(matrix))
    srtm_slr[] <- matrix[]
    #srtm_slr <- clump(srtm_slr)

    rm(matrix)

    clump <- 1:max

    flooded_land <- NULL

    for(i in 1:max){
      clump[i] <- srtm_slr[ocean_pos[i]]
      #if(i==1){flooded_land <- which(srtm_slr[]==clump[i])}else{flooded_land <- c(flooded_land,which(srtm_slr[]==clump[i]))}
      flooded_land <- c(flooded_land,which(srtm_slr[]==clump[i]))
      flooded_land <- unique(flooded_land)
    }


    #srtm@data@values[-flooded_land] <- NA
    values <- srtm[flooded_land]
    srtm[] <- NA
    srtm[flooded_land] <- values

    if(write==TRUE){writeRaster(srtm, paste(write_path,x,"_",y,sep=""), "GTiff", overwrite=TRUE)}
    #writeRaster(srtm, paste(write_path,x,"_",y,sep=""), "GTiff", options=c("COMPRESS=NONE", "TFw=YES"), overwrite=TRUE)
    return(srtm)
}