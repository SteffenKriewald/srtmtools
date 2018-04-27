#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void coord(double *x, double *y, double *x_res, double *y_res, double *vector, double *rows, double *cols, double *x_min, double *x_max, double *y_min, double *y_max, int *length)
{
	int i;
	double j;
	for(i=0; i<*length; i++){
	    j = ceil( (vector[i] / *cols) -1);
	    //printf("j: %f\n",j);
	    y_max[i] = *y - (j * *y_res);
	    y_min[i] = y_max[i] - *y_res;
	    x_min[i] = *x + ( ( vector[i] - (j * *cols) - 1 ) * *x_res );
	    x_max[i] = x_min[i] + *x_res;
	}
}


void sample( double *x_min, double *x_max, double *y_min, double *y_max, double *x, double *y, double *x_res, double *y_res, int *rows, int *cols, int *length, double *vector, double *vector_x_max, double *vector_y_min, double *vector_x_res, double *vector_y_res, double *area, int *length_stop)
{
	int i, min,max,col_min,col_max,row_min,row_max;
	double per_x, per_y,tmp;
	for(i=0; i<*length; i++){
//	for(i=3617; i<3680; i++){
	  
	col_min = (int) floor( floor( ((x_min[i] - *x) / *vector_x_res )*10 )/10);
	row_min = (int) floor( floor( ((*y - y_max[i]) / *vector_y_res )*10 )/10); 

	col_max = (int) floor( floor( ((x_max[i] - *x) / *vector_x_res )*10 )/10);
	row_max = (int) floor( floor( ((*y - y_min[i]) / *vector_y_res )*10 )/10); 

	//if(col_min >= *cols)printf("i: %d col_min: %d\n",i,col_min);
	//if(col_max >= *cols)printf("i: %d col_max: %d\n",i,col_max);
	//if(row_min >= *rows)printf("i: %d row_min: %d\n",i,row_min);
	//if(row_max >= *rows)printf("i: %d row_max: %d\n",i,row_max);

	//printf("x_min:%f y_max:%f x_max:%f y_min:%f\n",x_min[i],y_max[i],x_max[i],y_min[i]);
	//printf("y=%f\n",*y);
	//printf("v_y_res=%f\n", *vector_y_res);
	//printf("y-ymax=%f\n",*y - y_max[i]);
	//printf("round=%f\n",round( ((*y - y_max[i]) / *vector_y_res )*10 ));
	
	//if(col_min < 0)printf("i: %d col_min: %d\n",i,col_min);
	//if(col_max < 0)printf("i: %d col_max: %d\n",i,col_max);
	//if(row_min < 0)printf("i: %d row_min: %d\n",i,row_min);
	//if(row_max < 0)printf("i: %d row_max: %d\n",i,row_max);

	
	if( ((0 <= col_min && col_min < *cols) && (0 <= row_min && row_min < *rows)) || ((0 <= col_max && col_max < *cols) && (0 <= row_max && row_max < *rows)) )
	  {
	 
	   min = col_min + row_min* *cols;
	   max = col_max + row_max* *cols;
	   
	   if( (col_min == col_max) && (row_min == row_max) )
	    { 
	      vector[min]=vector[min] + area[i];	      
	    } 
	  else
	    { 
	      
	      // If min is part of the coarse raster
	      if((0 <= col_min && col_min < *cols) && (0 <= row_min && row_min < *rows))
		{
		per_x = (vector_x_max[min] - x_min[i]) / *x_res;
		if(per_x > 1)per_x=1;
		//printf("per_x: v_x_max_min=%f x_min_i=%f xres=%f\n",vector_x_max[min],x_min[i],*x_res);
		  
		per_y = (y_max[i] - vector_y_min[min]) / *y_res;
		if(per_y > 1)per_y=1;
		//printf("per_y: y_max_i=%f vector_y_min=%f yres=%f\n",y_max[i],vector_y_min[min],*y_res);
		//printf("min i:%d per_x=%f per_y=%f\n",i,per_x,per_y);
		
		// the first one
		vector[min] = vector[min] + ( per_x * per_y * area[i] );
		// one to the right
		if( col_min + 1 < *cols) vector[min+1] = vector[min+1] + ( (1.0 - per_x) * per_y * area[i] );
		// one down
		if( row_min + 1 < *rows )
		  {
		  vector[min + *cols] = vector[min + *cols] + ( per_x * (1.0 - per_y) * area[i] );
		  // still down and one to the right
		  if( col_min + 1 < *cols ) vector[min + *cols + 1] = vector[min + *cols + 1] + ( (1.0 - per_x) * (1.0 - per_y) * area[i] );
		  }
		}
	      // otherwise is only max part of the coarse raster
	      else
		{
		per_x = (x_max[i] - (vector_x_max[max] - *vector_x_res) ) / *x_res;
		if(per_x > 1)per_x=1;
		//printf("per_x: v_x_max_max=%f vxres=%f x_max_i=%f xres=%f\n",vector_x_max[max],*vector_x_res,x_max[i],*x_res);
		  
		per_y = ( (vector_y_min[max] + *vector_y_res) - y_min[i]) / *y_res;
		if(per_y > 1)per_y=1;
		//printf("per_y: y_max_i=%f vector_y_min=%f yres=%f\n",y_max[i],vector_y_min[min],*y_res);
		//printf("max i:%d per_x=%f per_y=%f\n",i,per_x,per_y);
	      
	        // the last one
		vector[max] = vector[max] + ( per_x * per_y * area[i] );
		// one to the left
		if( col_max - 1 >= 0) vector[max-1] = vector[max-1] + ( (1.0 - per_x) * per_y * area[i] );
		// one up
		if( row_max - 1 >= 0)
		  {
		  vector[max - *cols] = vector[max - *cols] + ( per_x * (1.0 - per_y) * area[i] );
		  // still up and one to the lrftt
		  if( col_max - 1 >= 0) vector[max - *cols - 1] = vector[max - *cols - 1] + ( (1.0 - per_x) * (1.0 - per_y) * area[i] );
		  }
		
		}
		
		
	      
		
	      }
	  //printf("vector nachher: %f\n", vector[min]);
	  	  
	  }
	  //if(abs(tmp+area[i]-vector[min])>1e-06){ printf("Area loss at: %d\n",i);};
	 }
}

//[1] "!!! WARNING !!! Internal check failed !!! 18.133082456059 17.8961274236594"
//[1] "!!! WARNING !!! Internal check failed !!! 46.0528818176111 45.0932215762422



void sample_till_22_04_2013( double *x_min, double *x_max, double *y_min, double *y_max, double *x, double *y, double *x_res, double *y_res, double *rows, double *cols, int *length, double *vector, double *vector_x_max, double *vector_y_min, double *vector_x_res, double *vector_y_res, double *area, int *length_stop)
{
	int i, min,max;
	double per_x, per_y;
	for(i=0; i<*length; i++){
	//for(i=3990; i<3991; i++){//6744
	min = (int) ( floor( round( ((x_min[i] - *x) / *vector_x_res )*10 )/10) + ( ( floor( round( ((*y - y_max[i])  / *vector_y_res)*10)/10) ) * *cols) );
	max = (int) ( floor( round( ((x_max[i] - *x) / *vector_x_res )*10 )/10) + ( ( floor( round( ((*y - y_min[i])  / *vector_y_res)*10)/10) ) * *cols) );
	
	/* printf("i: %d\n", i);

	printf("y_max[i]: %.15f\n", y_max[i]);	
	
	printf("x: %d\n", *x);

	printf("min: %d\n", min);
	printf("x_min[i]: %f\n", x_min[i]);
	printf("v_x_res]: %f\n", *vector_x_res);
    printf("floor1: %f\n", floor(x_min[i] / *vector_x_res));
	printf("floor2: %f\n", floor(*x / *vector_x_res ));
	printf("dritte: %f\n",  ( ( floor(*y / *vector_y_res) - floor(y_max[i]  / *vector_y_res) ) * *cols) );
	
	printf("max: %d\n", max);
	printf("x_max[i]: %f\n", x_max[i]);
	printf("v_x_res]: %f\n", *vector_x_res);
	printf("floor1: %f\n", floor(x_max[i] / *vector_x_res));
	printf("floor2: %f\n", floor(*x / *vector_x_res ));
	printf("dritte: %f\n", ( ( floor(*y/ *vector_y_res ) - floor(y_min[i]  / *vector_y_res) ) * *cols) );
	
	*/	

	if(min == max)vector[min]=vector[min] + area[i] ;
	else{ per_x = (vector_x_max[min] - x_min[i]) / *x_res;

		//printf("i: %d\n", i);
	    //printf("vector vorher: %f\n", vector[min+1]);
	    //printf("area: %f\n", area[i]);
	    //printf("min: %d\n", min);
	    //printf("max: %d\n", max);
	    //printf("v_x_max[min]: %f\n", vector_x_max[min]);
	    //printf("v_x_max[min+1]: %f\n", vector_x_max[min+1]);
	    //printf("v_x_max[min-1]: %f\n", vector_x_max[min-1]);
	    //printf("x_min[i]: %f\n", x_min[i]);
	    //printf("x_max[i]: %f\n", x_max[i]);
	    //printf("x_res: %f\n", *x_res);
	    //printf("per_x: %f\n", per_x);
	    
	      if(per_x > 1)per_x=1;
	      per_y = (y_max[i] - vector_y_min[min]) / *y_res;

		//printf("per_x: %f\n", per_x);

	      if(per_y > 1)per_y=1;

	    /*printf("v_y_min[min]: %f\n", vector_y_min[min]);
	    printf("y_min[i]: %f\n", y_min[i]);
	    printf("y_max[i]: %f\n", y_max[i]);
	    printf("y_res: %f\n", *y_res);
	    printf("per_y: %f\n", per_y);
		*/
	      //if((min <= (int)*cols))vector[min] = vector[min] + ( per_x * per_y * area[i] );
		  if(min < *length_stop)vector[min] = vector[min] + ( per_x * per_y * area[i] );
	      //if((min+1 <= (int)*cols))vector[min+1] = vector[min+1] + ( (1.0 - per_x) * per_y * area[i] );
		  if(min+1 < *length_stop)vector[min+1] = vector[min+1] + ( (1.0 - per_x) * per_y * area[i] );
	      if((min + (int)*cols) < *length_stop)vector[min + (int)*cols] = vector[min + (int)*cols] + ( per_x * (1.0 - per_y) * area[i] ); // here it is
	      if((min + (int)*cols + 1) < *length_stop)vector[min + (int)*cols + 1] = vector[min + (int)*cols + 1] + ( (1.0 - per_x) * (1.0 - per_y) * area[i] );

	   printf("vector nachher: %f\n", vector[min+1]);
	  

	     }
	}
}






void sample_org( double *x_min, double *x_max, double *y_min, double *y_max, double *x, double *y, double *x_res, double *y_res, double *rows, double *cols, int *length, double *vector, double *vector_x_max, double *vector_y_min, double *vector_x_res, double *vector_y_res, double *area)
{
	int i, min,max;
	double per_x, per_y;
	for(i=0; i<*length; i++){
	//printf("y_max[i]: %.15f\n", y_max[i]);
	//min = (int) ( floor(x_min[i] / *vector_x_res) - floor(*x / *vector_x_res ) + ( ( floor(*y / *vector_y_res) - floor(y_max[i]  / *vector_y_res) ) * *cols) );
	min = (int) ( floor( round( ((x_min[i] - *x) / *vector_x_res )*10 )/10) + ( ( floor( round( ((*y - y_max[i])  / *vector_y_res)*10)/10) ) * *cols) );
	//printf("min: %d\n", min);
	//printf("x_min[i]: %f\n", x_min[i]);
	//printf("v_x_res]: %f\n", *vector_x_res);
	//printf("floor1: %f\n", floor(x_min[i] / *vector_x_res));
	//printf("floor2: %f\n", floor(*x / *vector_x_res ));
	//printf("dritte: %f\n",  ( ( floor(*y / *vector_y_res) - floor(y_max[i]  / *vector_y_res) ) * *cols) );
	//max = (int) ( floor(x_max[i] / *vector_x_res) - floor(*x / *vector_x_res ) + ( ( floor(*y/ *vector_y_res ) - floor(y_min[i]  / *vector_y_res) ) * *cols) );
	max = (int) ( floor( round( ((x_max[i] - *x) / *vector_x_res )*10 )/10) + ( ( floor( round( ((*y - y_min[i])  / *vector_y_res)*10)/10) ) * *cols) );
	//printf("max: %d\n", max);
	//printf("x_max[i]: %f\n", x_max[i]);
	//printf("v_x_res]: %f\n", *vector_x_res);
	//printf("floor1: %f\n", floor(x_max[i] / *vector_x_res));
	//printf("floor2: %f\n", floor(*x / *vector_x_res ));
	//printf("dritte: %f\n", ( ( floor(*y/ *vector_y_res ) - floor(y_min[i]  / *vector_y_res) ) * *cols) );
	//if(min == max)vector[min]=vector[min] + 1.0 ;
	if(min == max)vector[min]=vector[min] + area[i] ;
	else{ per_x = (vector_x_max[min] - x_min[i]) / *x_res;
	      //if(per_x < 0)per_x=1;
	      //printf("i: %d\n", i);
	      //printf("vector vorher: %f\n", vector[min+1]);
	      //printf("area: %f\n", area[i]);
	      //printf("min: %d\n", min);
	      //printf("max: %d\n", max);
	      //printf("v_x_max[min]: %f\n", vector_x_max[min]);
	      //printf("v_x_max[min+1]: %f\n", vector_x_max[min+1]);
	      //printf("v_x_max[min-1]: %f\n", vector_x_max[min-1]);
	      //printf("x_min[i]: %f\n", x_min[i]);
	      //printf("x_max[i]: %f\n", x_max[i]);
	      //printf("x_res: %f\n", *x_res);
	      //printf("per_x: %f\n", per_x);
	      if(per_x > 1)per_x=1;
	      //printf("per_x: %f\n", per_x);
	      //per_y = (vector_y_min[min] - y_min[i]) / *y_res;
	      per_y = (y_max[i] - vector_y_min[min]) / *y_res;
	      //per_y = (vector_y_min[min] - y_max[i]) / *y_res;  ### ??? ###
	      //printf("v_y_min[min]: %f\n", vector_y_min[min]);
	      //printf("y_min[i]: %f\n", y_min[i]);
	      //printf("y_max[i]: %f\n", y_max[i]);
	      //printf("y_res: %f\n", *y_res);
	      //if(per_y < 0)per_y=1;
	      //printf("per_y: %f\n", per_y);
	      if(per_y > 1)per_y=1;
	      //printf("per_y: %f\n", per_y);
	      vector[min] = vector[min] + ( per_x * per_y * area[i] );
	      vector[min+1] = vector[min+1] + ( (1.0 - per_x) * per_y * area[i] );
	      vector[min + (int)*cols] = vector[min + (int)*cols] + ( per_x * (1.0 - per_y) * area[i] );
	      vector[min + (int)*cols + 1] = vector[min + (int)*cols + 1] + ( (1.0 - per_x) * (1.0 - per_y) * area[i] );
	      //printf("vector nachher: %f\n", vector[min+1]);
	    }
	}
}

