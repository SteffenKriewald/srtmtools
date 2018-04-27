compare.border <- function(srtm,x,y,n,s,w,e,write_path){
    # compare border
    neighbour_border <- NULL
    if(sum(dir(write_path)==paste(x,"_",y,".tif",sep=""))==1){
      print(" ... load flooded border ...")
      #print(x)
      #print(y)
      neighbour_border <- flooded.region(x=x,y=y,n=n,s=s,w=w,e=e,write_path=write_path)
      flooded_land <- which(is.na(neighbour_border[])==FALSE)
      border <- crop(srtm,neighbour_border)
      border[flooded_land] <- NA
      srtm <- merge(border,srtm,overlap=FALSE)
    }
    return(list(srtm,neighbour_border))
}