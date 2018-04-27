
buffer <- function(m, sz)
{	
	m1 <- .C("ccaBuffED", m=as.integer(m), nr=as.integer(dim(m)[1]), nc=as.integer(dim(m)[2]), sz=as.integer(sz));
	m1 <- m1$m
	m1[which(m1<0)] <- -1;
	m1 <- matrix(m1, nrow=dim(m)[1]);
	return(m1);
}



coord <- function(x,y,x_res,y_res, vector, rows, cols)
{	
	length <- length(vector)
	x_min <- x_max <- y_min <- y_max <-  rep(0,length)
	m1 <- .C("coord", x=as.double(x), y=as.double(y), x_res=as.double(x_res), y_res=as.double(y_res), vect=as.double(vector), rows=as.double(rows), cols=as.double(cols), x_min=as.double(x_min), x_max=as.double(x_max), y_min=as.double(y_min), y_max=as.double(y_max), length=as.integer(length) );
	#m1 <- cbind(vector=m1$vector,x_min=m1$x_min,x_max=m1$x_max,y_min=m1$y_min,y_max=m1$y_max)
	return(m1);
}


#sample <- function(x_min, x_max, y_min, y_max, x, y, x_res, y_res, rows, cols, vector_x_max, vector_y_min, vector_x_res, vector_y_res)
#{
	
#	length <- length(x_min)
#	vector <- rep(0,length(vector_x_max))
#	m1 <- .C("sample", x_min=as.double(x_min), x_max=as.double(x_max), y_min=as.double(y_min), y_max=as.double(y_max), x=as.double(x), y=as.double(y), x_res=as.double(x_res), y_res=as.double(y_res), rows=as.double(rows), cols=as.double(cols), length=as.integer(length), vector=as.double(vector), vector_x_max=as.double(vector_x_max), vector_y_min=as.double(vector_y_min), vector_x_res=as.double(vector_x_res), vector_y_res=as.double(vector_y_res) );
#	#m1 <- cbind(vector=m1$vector,x_min=m1$x_min,x_max=m1$x_max,y_min=m1$y_min,y_max=m1$y_max)
#	return(m1);
#}


sample <- function(r1,r2,cell.constraint=0,area=TRUE,crop=TRUE)
{	
	# r1 = low resolution #### r2 =high resolution #

	# coordinates from all cells of r1 (low resolution)
	vr1 <- 1:(r1@nrows*r1@ncols)
	cr1 <- coord(x=r1@extent@xmin,y=r1@extent@ymax,x_res=res(r1)[1],y_res=res(r1)[2],vector=vr1,rows=r1@nrows,cols=r1@ncols)
	length_stop <- length(cr1$x_min)	

	# crop r2 to the extent of r1 if extent(r2)>extent(r1)
	if(crop==TRUE){
	  if(r1@extent@xmin  > r2@extent@xmin | r1@extent@xmax < r2@extent@xmax | r1@extent@ymin  > r2@extent@ymin | r1@extent@ymax < r2@extent@ymax ){	
		  print("crop raster")
		  r2 <- crop(r2,r1,snap="in")
		  print("start sampling")
		  }
	}
		  
	# coordinates from all cells > cell.constraint of r2 (high resolution)
	vr2 <- which(r2[]>cell.constraint)
	cr2 <- coord(x=r2@extent@xmin,y=r2@extent@ymax,x_res=res(r2)[1],y_res=res(r2)[2],vector=vr2,rows=r2@nrows,cols=r2@ncols)

	if(area==TRUE){area <- area(r2)[vr2]}else{
	  area <- r2[vr2]
	}
  
  # internal check variable
  sum_org <- sum(area)
	
  length <- length(cr2$x_min)
	vector <- rep(0,length(vr1))
	m1 <- .C("sample", x_min=as.double(cr2$x_min), x_max=as.double(cr2$x_max), y_min=as.double(cr2$y_min), y_max=as.double(cr2$y_max), x=as.double(r1@extent@xmin), y=as.double(r1@extent@ymax), x_res=as.double(res(r2)[1]), y_res=as.double(res(r2)[2]), rows=as.integer(r1@nrows), cols=as.integer(r1@ncols), length=as.integer(length), vector=as.double(vector), vector_x_max=as.double(cr1$x_max), vector_y_min=as.double(cr1$y_min), vector_x_res=as.double(res(r1)[1]), vector_y_res=as.double(res(r1)[2]), area=as.double(area), length_stop=as.integer(length_stop));
	#m1 <- cbind(vector=m1$vector,x_min=m1$x_min,x_max=m1$x_max,y_min=m1$y_min,y_max=m1$y_max)
	
	rs <- r1
	rs[] <- matrix(m1$vector,ncol=r1@ncols,byrow=TRUE)
  
  # internal check result
  sum_new <- sum(rs[])
  if((sum_org-sum_new)>0.01*(sum_org)){print(paste("!!! WARNING !!! Internal variation is bigger than 1% !!!",sum_org,sum_new))}
  
	
	return(list(m1,rs,c(sum_org,sum_new)));
}
