/* C functions in the sROC library*/

#include <R.h>
#include <Rmath.h>

void NKern(double *x, int *n, double *xgrid, int *ngrid, double *bw, double *Fhat){
  int i, j;
  double d, ksum;
  for (i=0; i < *ngrid; i++){
    ksum=0.0;
    for (j=0; j <*n; j++){
      d = xgrid[i]-x[j];
      ksum += pnorm(d / *bw, 0, 1, 1, 0);
    }
    Fhat[i] = ksum/(*n);  
  }
}

double Epan(double a){
    double k;
    if(a < -1.0) k=0.0;
		else if(a > 1.0) k=1.0;
    		else k= 0.75*a - 0.25*pow(a,3) + 0.5;
    return k;
}


void EKern(double *x, int *n, double *xgrid, int *ngrid, double *bw, double *Fhat){
  int i, j;
  double d, ksum;
  for (i=0; i < *ngrid; i++){
    ksum=0.0;
    for (j=0; j <*n; j++){
      d = xgrid[i]-x[j];
      ksum += Epan(d / *bw);
    }
    Fhat[i] = ksum/(*n);  
  }
}

double pnorm2(double x){
    double y;
	y =dnorm(x, 0, 1, 0)*(pow(x,2)-1.0);
    return y;
}

void phi2(double *x, int *n, double *xgrid, int *ngrid, double *bw, double *phi){
  int i, j;
  double d, ksum;
  for (i=0; i < *ngrid; i++){
    ksum=0.0;
    for (j=0; j <*n; j++){
      d = xgrid[i]-x[j];
      ksum += pnorm2(d / *bw);
    }
    phi[i] = ksum/((*n)*pow(*bw,3));  
  }
}
