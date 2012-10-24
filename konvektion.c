#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define FALSE 0
#define TRUE 1

static int sortfunc (const void *a, const void *b) {
  
  if(*((double*)a)<*((double*)b)) {
    return 1;
  }
  if(*((double*)a)>*((double*)b)) {
    return -1;
  }  
  return 0;
  
}
  

int main() {
  
  double g=9.8065;
  double cp=1004;
  double H=100;
  double p0=1000;
  double deltaT=0;
  double deltat=10;
  double Ra=287;
  
  int timesteps = 0;
  
  int ilev=0;
  int ilyr=0;
  int nlyr=10;
  int nlev=nlyr+1;
  int instabil=FALSE;
  
  double deltap=p0/nlyr;
  double kappa=Ra/cp;
  
  double *p=calloc(nlev,sizeof(double));
  double *T=calloc(nlyr,sizeof(double));
  double *theta=calloc(nlev,sizeof(double));
  double *plyr=calloc(nlev,sizeof(double));

  
  for(ilev=0; ilev<nlev; ilev++) {
    p[ilev]=p0*ilev/(nlev-1);
    
   printf("%d %f\n", ilev, p[ilev]);
   
  }
  
  for(ilyr=0; ilyr<nlyr; ilyr+=1) {
    T[ilyr]=255;
    
  }
  
  while (T[nlyr-1]<400) {
    printf("\nNew time %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
    timesteps++;
    deltaT=H/(deltap*100.0)/ cp*g*deltat;
    //deltaT = H * g / (deltap * cp) * deltat;
    T[nlyr-1]+= deltaT;
  
  
    for(ilyr=0; ilyr<nlyr; ilyr++) {
      plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
      //printf("%d %f\n", ilyr, plyr[ilyr]);
    }
  
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {
      theta[ilyr]=pow(p0/plyr[ilyr], kappa)*T[ilyr];
      //printf("%d %f\n", ilyr, theta[ilyr]);
    }
    
    instabil = FALSE;
    
    for (ilyr=0; ilyr<nlyr;ilyr++) {
      if (theta[ilyr+1]>theta[ilyr]) {
	instabil = TRUE;
	break;
      }
    }
      
    if (instabil) {
      qsort (theta, nlyr, sizeof(double), sortfunc);
      printf ("Konvektion :-)\n");
    }
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {
      //theta[ilyr]=pow(p0/plyr[ilyr], kappa)*T[ilyr];
      T[ilyr] = theta[ilyr] * pow(plyr[ilyr]/p0, kappa);
      //printf("%d %f\n", ilyr, theta[ilyr]);
    }
  }
  return 0;
}

