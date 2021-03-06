cca <- function(pop, s=1, mode=2, count.cells=FALSE, count.max=ncol(pop)*3, scale=FALSE, res=0.00277778, lat=60){
	#do checks
	stopifnot(is.numeric(pop))
	stopifnot(is.matrix(pop))
	stopifnot(is.numeric(s))
	stopifnot(mode==1 | mode==2 | mode==3)
	the.pop <- as.integer(t(pop))
	clu <- as.integer(rep(0, ncol(pop)*nrow(pop)))
        count.max <- as.integer(count.max)
	count <- as.double(rep(0, count.max))
	if(count.cells==TRUE && scale==FALSE){
	  count <- as.integer(count)
	  out <- .C("callburn_count",  s=as.integer(s), xmax=nrow(pop), ymax=ncol(pop), mode=as.integer(mode)[1], pop=the.pop, clu=clu, count=count, count.max=count.max,CLASSES=c("integer", "integer", "integer", "integer", "integer","integer"))
	  print("max ist")
	  print(max(out$clu))
	  return(list(clusters=matrix(out$clu, ncol=ncol(pop), byrow=TRUE), cluster.count=out$count[1:max(out$clu)]))
	  }
	if(count.cells==FALSE && scale==FALSE){
	  out <- .C("callburn",  s=as.integer(s), xmax=nrow(pop), ymax=ncol(pop), mode=as.integer(mode)[1], pop=the.pop, clu=clu, count=count, count.max=count.max,CLASSES=c("integer", "integer", "integer", "integer", "integer","integer"))
	  return(clusters=matrix(out$clu, ncol=ncol(pop), byrow=TRUE))
	  }
	if(count.cells==TRUE && scale==TRUE){
	  out <- .C("callburn_scale",  s=as.double(s), xmax=nrow(pop), ymax=ncol(pop), mode=as.integer(mode)[1], pop=the.pop, clu=clu, res=res, lat=lat, count=count, count.max=count.max, CLASSES=c("double", "integer", "integer", "integer", "integer","integer","double","double","double","integer"))
	  return(list(clusters=matrix(out$clu, ncol=ncol(pop), byrow=TRUE), cluster.count=out$count[1:max(out$clu)]))
	  }
}

