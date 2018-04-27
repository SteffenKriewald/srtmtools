#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>   /* INT_MIN und INT_MAX */

//void srtmtrim(double *raster, int nrc[0], int nrc[1], int *n, int *s, int *w, int *e)
void srtmtrim(double *raster, int *nrc, int *coord)
{
  int r,c,tncol,tnrow,var,na,north,south,west,east;
  double sum,mult;
  
  //printf("int Gr��e : %d Byte\n",sizeof(int));
  //printf("Wertebereich von %d bis %d\n",INT_MIN,INT_MAX);
  


  //printf("na=%f\n",raster[1]);
  na=-999;

  // find north
  r=0;
  mult=(na)* nrc[1];
  //printf("mult=%f\n",mult);
  do{
     sum=0;
     for(c=0; c<nrc[1]; c++)
      {
	sum = sum + raster[r* nrc[1] + c];
	//printf("sum=%f\n",sum);
      }
      //printf("r=%i\n",r);
      r++;
      //printf("r=%i\n",r);
  }while(sum == mult);
  //printf("r=%i\n",r);
  r=r-1;
  north=r;
  printf("north=%i\n",north);
  
    
  // find south
  r=nrc[0]-1;
  //mult=(na)* nrc[1];
  do{
     sum=0;
     for(c=0; c<nrc[1]; c++)
      {
	sum = sum + raster[r* nrc[1] + c]; 
      }
      r--;
  }while(sum == mult);
  r=r+1;
  south=r;
  printf("south=%i\n",south);


  tnrow = south - north +1;
  // find west
  c=0;
  //printf("col_s=%i,\n",c);
  mult=(na)* (tnrow);
  //printf("mult=%f\n",mult);
  do{
     sum=0;
       //printf("n=%i\n",*n);
       //printf("s=%i\n",*s);
     //for(r=2; r=2; r++)
     var = north;
    //printf("north=%i,\n",north);
    //printf("row_n=%i,\n",var);
     while(var <= south)
      {
	//printf("row=%i,\n",var);
	//printf("col=%i,\n",c);
	//printf("raster=%f,\n",raster[var* nrc[1] + c]);
	sum = sum + raster[var* nrc[1] + c];
	var++;
      }
      //printf("sum=%f\n",sum);
      c++;
  }while(sum == mult);
  c=c-1;
  west=c;
  printf("west=%i\n",west);
  
  // find east
  c=nrc[1]-1;
  do{
     sum=0;
     //for(r=*n; r=*s; r++)
     var = north;
     while(var <= south)
      {
	sum = sum + raster[var* nrc[1] + c];
	var++;
      }
      //printf("sum=%f\n",sum);
      c--;
  }while(sum == mult);
  c=c+1;
  east=c;
  printf("east=%i\n",east);


  tncol = east - west + 1;
  // collect trimed data
  for(r=north;r<=south;r++){
  //while(y <= *)
    for(c=west;c<=east;c++){
      raster[(r-north) * tncol + (c-west)] = raster[r* nrc[1] + c];
      //printf("%f ",raster[(r-north) * tncol + (c-west)]);
    }
    //printf("\n");
  }
  
 /* ncol = &tncol;
  printf("ncol=%i\n",nrc[1]);
  nrow = &tnrow;
  printf("nrow=%i\n",nrc[0]);
*/
  nrc[0]=tnrow;
  nrc[1]=tncol;

  coord[0]=north;
  coord[1]=south;
  coord[2]=west;
  coord[3]=east;


}