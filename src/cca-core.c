
#include <R.h>
#include <math.h>
#include <stdio.h>

#define PI 3.1415926

//Required 
void getrow(int *x, int *xmax, int *ymax, int *pop, int *ret) {

	int y;
	for(y=0; y<= *ymax; y++){
		ret[y] = pop[*x * *ymax + y];
	}
}
void getcol(int *y, int *xmax, int *ymax, int *pop, int *ret) {
	int x;
	for(x=0; x<= *xmax; x++){
		ret[x] = pop[x * *ymax + *y];
	}
}
void getblock(int *y, int *x, int *dy, int *dx, int *xmax, int *ymax, int *pop, int *ret) {
	int e, f;
	for(e=0; e< *dx; e++){
		for(f=0; f< *dy; f++){
			ret[f * *dy + e] = pop[f * *ymax + e];
		}
	}
}


void burnn(int *x,int *y,int *c, int *xmax, int *ymax, int *pop,int *clu) {
	// burn nearest neigbors
	//Rprintf("xycrxmaxymax\t%i\t%i\t%i\t%i\t%i\t%i\n",*x,*y,*c,1,*xmax,*ymax);
	int d, e,f,g;
	d=*x;					// steps to the left
	//while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
	while(d>=0 && pop[d * *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		d--;

	}
	e=*x+1;					// steps to the right
	//while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
	while(e<*xmax&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		e++;
	}
	//Rprintf("  de pop\t%i\t%i\t", d,e);
	//for(f=d-1;f<=e+1;f++) Rprintf("%i\t",pop[f * *ymax + *y]);

	for(f=*y+1;f>=*y-1;f=f-2){// one step up or down
		if(f<*ymax && f>=0){				
			for(g=d+1; g<e; g++){ //I think we must not include d and e to avoid diagonals
				if(clu[g * *ymax + f]==0&&pop[g * *ymax + f]>0) {
					burnn(&g,&f,c,xmax,ymax,pop,clu);
				}
			}
		}
	}
}

void burns(int *pop,int *clu,int *x,int *y,int *c,int *s, int *xmax, int *ymax) {
	// burn s shells
	//Rprintf("xycrxmaxymax\t%i\t%i\t%i\t%i\t%i\t%i\n",*x,*y,*c,*s,*xmax,*ymax);
	long a,b;
	int dx,dy;
	int d, e,g;
	d=*x;					// steps to the left
	while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		d--;

	}
	e=*x+1;					// steps to the right
	while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		e++;
	}
	for(g=d+1; g<e; g++){ // go allong all marked cells and use algorithm from diego
		for(a=-*s;a<=*s;a++) {
			dx=g+a;
			//Rprintf("dx\t%i\n",dx);
			if(dx>=0&&dx<*xmax) {
				for(b=-*s;b<= *s;b++) {
					dy=*y+b;
					//Rprintf("dy\t%i\n",dy);
					if(dy>=0&&dy<*ymax) {
						//Rprintf("clu[%i,%i], pop[%i,%i]\t%i\t%i\n",dx, dy, dx,dy,clu[dx * *ymax +dy],pop[dx * *ymax +dy]);
						if(clu[dx * *ymax +dy]==0&&pop[dx * *ymax +dy]>0) {
							burns(pop,clu,&dx,&dy,c,s,xmax,ymax);
						}
					}
				}
			}
		}
	}
}

void burnr(int *pop,int *clu,int *x,int *y,int *c,int *r, int *xmax, int *ymax) {
	// burn with radius r
	//	r=1 should correspond burnn
	long a,b;
	int dx,dy;
	double rd;
	int d, e,g;
	d=*x;					// steps to the left
	while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		d--;

	}
	e=*x+1;					// steps to the right
	while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		e++;
	}
	for(g=d+1; g<e; g++){ // go allong all marked cells and use algorithm from diego

		for(a=-*r;a<=*r;a++) {
			for(b=-*r;b<=*r;b++) {
				rd=sqrt((double)a*(double)a+(double)b*(double)b);
				//d=1;
				if(rd<=*r) {
					dx=g+a;
					dy=*y+b;
					if((dx>=0&&dx<*xmax)&&(dy>=0&&dy<*ymax)) {
						if(clu[dx * *ymax + dy]==0&& pop[dx * *ymax + dy]>0) {
							burnr(pop,clu,&dx,&dy,c,r, xmax, ymax);
						}
					}
				}
			}
		}
	}
}

 /* Note: R fills matrices per column, so first index is row,
     second is column */
void callburn(int *s, int *xmax, int *ymax, int *mode, int *pop,int *clu) {
	int a,b;
	int c;
	//long s=5;
	c=0;
    	//printf("max: %d %d\n",*xmax, *ymax);
	for(a=0;a<*xmax;a++) {
		//printf("s: %ld\tcell: %ld %ld: %ld\n",*s,a, b,pop[a * *ymax+b]);
		for(b=0;b<*ymax;b++) {
			if(pop[a* *ymax + b]>0&&clu[a* *ymax + b]==0) {
				c++;
				//printf("burn:\t%d\n",c);
				if(*mode == 1){
					burnn(&a,&b,&c,xmax,ymax,pop,clu);		// burn nearest
				} else if(*mode == 2){
					burns(pop,clu,&a,&b,&c,s,xmax,ymax);		// burn shells
				} else if(*mode == 3){
					burnr(pop,clu,&a,&b,&c,s,xmax,ymax);		// burn radius
				} else {
					printf("unknown mode: %d\n", *mode);
				}
			}
		}
	}
}

			


void burnn_count(int *x,int *y,int *c, int *xmax, int *ymax, int *pop,int *clu, int *count) {
	// burn nearest neigbors
	//Rprintf("xycrxmaxymax\t%i\t%i\t%i\t%i\t%i\t%i\n",*x,*y,*c,1,*xmax,*ymax);
	int d, e,f,g;
	d=*x;					// steps to the left
	//while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
	while(d>=0 && pop[d * *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		count[*c-1]++;
		d--;
	}
	e=*x+1;					// steps to the right
	//while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
	while(e<*xmax&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		count[*c-1]++;
		e++;
	}
	//Rprintf("  de pop\t%i\t%i\t", d,e);
	//for(f=d-1;f<=e+1;f++) Rprintf("%i\t",pop[f * *ymax + *y]);

	for(f=*y+1;f>=*y-1;f=f-2){// one step up or down
		if(f<*ymax && f>=0){				
			for(g=d+1; g<e; g++){ //I think we must not include d and e to avoid diagonals
				if(clu[g * *ymax + f]==0&&pop[g * *ymax + f]>0) {
					burnn_count(&g,&f,c,xmax,ymax,pop,clu,count);
				}
			}
		}
	}
}

void burns_count(int *pop,int *clu,int *x,int *y,int *c,int *s, int *xmax, int *ymax, int *count) {
	// burn s shells
	//Rprintf("xycrxmaxymax\t%i\t%i\t%i\t%i\t%i\t%i\n",*x,*y,*c,*s,*xmax,*ymax);
	long a,b;
	int dx,dy;
	int d, e,g;
	d=*x;					// steps to the left
	while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		count[*c-1]++;
		d--;

	}
	e=*x+1;					// steps to the right
	while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		count[*c-1]++;
		e++;
	}
	for(g=d+1; g<e; g++){ // go allong all marked cells and use algorithm from diego
		for(a=-*s;a<=*s;a++) {
			dx=g+a;
			//Rprintf("dx\t%i\n",dx);
			if(dx>=0&&dx<*xmax) {
				for(b=-*s;b<= *s;b++) {
					dy=*y+b;
					//Rprintf("dy\t%i\n",dy);
					if(dy>=0&&dy<*ymax) {
						//Rprintf("clu[%i,%i], pop[%i,%i]\t%i\t%i\n",dx, dy, dx,dy,clu[dx * *ymax +dy],pop[dx * *ymax +dy]);
						if(clu[dx * *ymax +dy]==0&&pop[dx * *ymax +dy]>0) {
							burns_count(pop,clu,&dx,&dy,c,s,xmax,ymax,count);
						}
					}
				}
			}
		}
	}
}

void burnr_count(int *pop,int *clu,int *x,int *y,int *c,int *r, int *xmax, int *ymax, int *count) {
	// burn with radius r
	//	r=1 should correspond burnn
	long a,b;
	int dx,dy;
	double rd;
	int d, e,g;
	d=*x;					// steps to the left
	while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		count[*c-1]++;
		d--;

	}
	e=*x+1;					// steps to the right
	while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		count[*c-1]++;
		e++;
	}
	for(g=d+1; g<e; g++){ // go allong all marked cells and use algorithm from diego

		for(a=-*r;a<=*r;a++) {
			for(b=-*r;b<=*r;b++) {
				rd=sqrt((double)a*(double)a+(double)b*(double)b);
				//d=1;
				if(rd<=*r) {
					dx=g+a;
					dy=*y+b;
					if((dx>=0&&dx<*xmax)&&(dy>=0&&dy<*ymax)) {
						if(clu[dx * *ymax + dy]==0&& pop[dx * *ymax + dy]>0) {
							burnr_count(pop,clu,&dx,&dy,c,r, xmax, ymax, count);
						}
					}
				}
			}
		}
	}
}

 /* Note: R fills matrices per column, so first index is row,
     second is column */
void callburn_count(int *s, int *xmax, int *ymax, int *mode, int *pop,int *clu, int *count, int *countmax) {
	int a,b;
	int c;
	//long s=5;
	c=0;
    	//printf("max: %d %d\n",*xmax, *ymax);
	for(a=0;a<*xmax;a++) {
		//printf("s: %ld\tcell: %ld %ld: %ld\n",*s,a, b,pop[a * *ymax+b]);
		for(b=0;b<*ymax;b++) {
			if(pop[a* *ymax + b]>0&&clu[a* *ymax + b]==0) {
				c++;
				if(*countmax <= c){
				    Rprintf("count.max is too little\n");
				    return;
                                }
				Rprintf("burn:\t%d\n",c);
				if(*mode == 1){
					burnn_count(&a,&b,&c,xmax,ymax,pop,clu,count);	// burn nearest
				} else if(*mode == 2){
					burns_count(pop,clu,&a,&b,&c,s,xmax,ymax,count);		// burn shells
				} else if(*mode == 3){
					burnr_count(pop,clu,&a,&b,&c,s,xmax,ymax,count);		// burn radius
				} else {
					printf("unknown mode: %d\n", *mode);
				}
			}
		}
	}
}

void burnr_scale(int *pop,int *clu,int *x,int *y,int *c,double *r, int *xmax, int *ymax, double *res, double *lat, double *count) {
	// burn with radius r
	//	r=1 should correspond burnn
	long a,b;
	int dx,dy;
	double rd,lat_local,scale;
	int d, e,g;
	
	lat_local = *lat;
	
	d=*x;					// steps to the left ?? runter
	while(d>=0 && clu[d* *ymax + *y]==0 && pop[d* *ymax + *y]>0) {
		clu[d * *ymax + *y]=*c;
		scale=cos(lat_local*PI/180);
		count[*c-1]=count[*c-1]+(1*scale);
		d--;

	}
	e=*x+1;					// steps to the right ?? rauf
	while(e<*xmax&&clu[e * *ymax + *y]==0&&pop[e * *ymax + *y]>0) {
		clu[e * *ymax + *y]=*c;
		scale=cos(lat_local*PI/180);
		count[*c-1]=count[*c-1]+(1*scale);
		e++;
	}
	for(g=d+1; g<e; g++){ // go allong all marked cells and use algorithm from diego
	//printf("*r: %f, \n", *r);
		for(a=-*r;a<=*r;a++) {
			lat_local=*res*a+lat_local;
	//		printf("lat_local: %f, \n", lat_local);
			scale=cos(lat_local*PI/180);
	//		printf("scale: %f, \n", scale);
			for(b=-(int)(*r/scale);b<=(int)(*r/scale);b++) {
				rd=sqrt((double)a*(double)a+(double)b*scale*(double)b*scale);
				//printf("rd: %f, \n", rd);
				//d=1;
				if(rd<=*r) {
					dx=g+a;
					dy=*y+b;
					if((dx>=0&&dx<*xmax)&&(dy>=0&&dy<*ymax)) {
						if(clu[dx * *ymax + dy]==0&& pop[dx * *ymax + dy]>0) {
							burnr_scale(pop,clu,&dx,&dy,c,r, xmax, ymax,res,&lat_local,count);
						}
					}
				}
			}
		}
	}
	
}


void callburn_scale(double *s, int *xmax, int *ymax, int *mode, int *pop,int *clu, double *res, double *lat, double *count, int *countmax) {
	int a,b;
	int c;
	double lat_local;
	//long s=5;
	c=0;
    	//printf("max: %d %d\n",*xmax, *ymax);
	for(a=0;a<*xmax;a++) {
		//printf("s: %ld\tcell: %ld %ld: %ld\n",*s,a, b,pop[a * *ymax+b]);
		for(b=0;b<*ymax;b++) {
			if(pop[a* *ymax + b]>0&&clu[a* *ymax + b]==0) {
				lat_local = *lat - (*res/2) - (*res *(double)a);
				c++;
				if(*countmax <= c){
				    Rprintf("count.max is too little\n");
				    return;
				}
				Rprintf("burn:\t%d\n",c);
				if(*mode == 1){
					//burnn_count(&a,&b,&c,xmax,ymax,pop,clu,count);	// burn nearest
				} else if(*mode == 2){
					//burns_count(pop,clu,&a,&b,&c,s,xmax,ymax,count);		// burn shells
				} else if(*mode == 3){
					burnr_scale(pop,clu,&a,&b,&c,s,xmax,ymax,res,&lat_local,count);		// burn radius
				} else {
					printf("unknown mode: %d\n", *mode);
				}
			}
		}
	}
}
