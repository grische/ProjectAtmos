#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <plplot/plplot.h>
#include "schwarz.h"

#define FALSE 0
#define TRUE 1

#define R 8.31447     /* J/(mol K) */
#define M 0.0289644   /* kg / mol */
#define g 9.8065      /* m/s^2 */
#define EPSILON 0.99   /* break condition TMAX*EPS */

#define TMAX 255.0    /* K */
#define TMIN 100.0    /* K */

static int sortfunc (const void *a, const void *b) {             /* Defining sortfunc */
  
  if(*((double*)a)<*((double*)b)) {
    return 1;
  }
  if(*((double*)a)>*((double*)b)) {
    return -1;
  }  
  return 0;
  
}

double alpha (double T) {
  double x = -( R * T )/( g * M );
  return x;
}


int main() {                                                     /* Definition of Variables and Quantities  */
  
  double cp=1004;  /* J/kg K */
  double Esol=235;  /* W/m^2 = J*/
  double p0=1000;  /* hPa */
  double deltaTsurf=0; /* K */
  double deltat=100;  /* s */
  double Ra=287; /* J/kg K */
  double H=0; /* [Esol] */
  double Tsurf=0; /* K */

  int timesteps = 0;
  int ilev=0;
  int ilyr=0;
  int nlyr=10;  /* Number of Layers */
  int nlev=nlyr+1;
  int instabil=FALSE;
  
  double deltap=p0/nlyr;
  double kappa=Ra/cp;
  
  double *p=calloc(nlev,sizeof(double));
  double *T=calloc(nlyr,sizeof(double));
  double *theta=calloc(nlev,sizeof(double));
  double *plyr=calloc(nlyr,sizeof(double));
  double *z=calloc(nlyr,sizeof(double));
  double *deltaz=calloc(nlyr,sizeof(double));
  double *eup=calloc(nlev, sizeof(double));
  double *edn=calloc(nlev, sizeof(double));
  double *enet=calloc(nlyr, sizeof(double));
  double *deltaT=calloc(nlyr, sizeof(double));


  plscolbg (255, 255, 255);   /* background color white */
  plscol0  (15, 0, 0, 0);     /* set color 15 to black  */
  plsdev ("xwin"); /* if not called, the user is asked! */
  plinit ();
  
  plssub(2,1);


  for(ilev=0; ilev<nlev; ilev++) {                              /* Calculation of the Pressure at the Levels p[ilev] */
    p[ilev]=p0*ilev/(nlev-1);
    
    //printf("%d %f\n", ilev, p[ilev]);
   
  }
  
  for(ilyr=0; ilyr<nlyr; ilyr+=1) {                             /* Temperature for all Layers 255K */
    T[ilyr]=TMIN;
    printf("%f\n", T[ilyr]);
  }
 
  
  Tsurf=100;         /* Tsurf */

  while (Tsurf<400) {        /* Begrenzung */


    H=Esol-((5.67E-8)*(pow(Tsurf,4)));
  deltaTsurf=(H*g*deltat)/((deltap*100.0)*cp);                    /* Equation for Temperature of Tsurf */
  Tsurf += deltaTsurf;
  printf("deltaTsurf%f\n", deltaTsurf);

  printf("H %f, Tsurf %f\n",H, Tsurf);

  //schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup);

  schwarzschild(1.0, T, nlev, Tsurf, edn, eup);

    /* E-netto */
  for(ilyr=0; ilyr<nlyr; ilyr++) {
    enet[ilyr] = eup[ilyr+1] + edn[ilyr] - eup[ilyr] - edn[ilyr+1];
  printf("enet%f\n", enet[ilyr]);
  }



  for (ilyr=0; ilyr<nlyr; ilyr++)  {
    deltaT[ilyr]=(enet[ilyr]*g*deltat)/((deltap*100.0)*cp);        /* Temperature gain for each layer */
    T[ilyr] +=deltaT[ilyr];
    printf("T%f\n", T[ilyr]);
}

  
  
    for(ilyr=0; ilyr<nlyr; ilyr++) {                            /* Calculation of the Pressure in the Layers*/
      plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
      //printf("%d %f\n", ilyr, plyr[ilyr]);
    }
  
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {                           /* Conversion from T to theta */
      theta[ilyr]=pow(p0/plyr[ilyr], kappa)*T[ilyr];
      //printf("%d %f\n", ilyr, theta[ilyr]);
    }
    
    instabil = FALSE;
    
    for (ilyr=0; ilyr<nlyr;ilyr++) {                            /* Testing for Instability */
      if (theta[ilyr+1]>theta[ilyr]) {
	instabil = TRUE;
	break;
      }
    }

    if (instabil) {                                             /* Convection - Sorting of Layers according to theta */
      qsort (theta, nlyr, sizeof(double), sortfunc);
      //  printf ("Konvektion :-)\n");
    }
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {                           /* Conversion from theta to T */
      T[ilyr] = theta[ilyr] /  (pow(p0/plyr[ilyr], kappa));

      //printf("%d %f\n", ilyr, theta[ilyr]);
    }

    if(timesteps % (int)(10000000/deltat/nlyr) == 0) {
      printf ("timestep %d\n", timesteps);
    }

      pladv(1);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( 0, TMAX, p0, 0 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", (TMAX-TMIN)/4.0 , 0, "bcnst", 150.0, 0 );
      pllab ("temperature [K]", "p [hPa]", "");  /* axis labels     */

      plcol0 (9);                         /* color blue  */
      plline (nlyr, T, plyr);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */


      for (ilyr=0;ilyr<nlyr; ilyr++) {
	deltaz[ilyr]=(Ra*T[ilyr]*deltap)/(plyr[ilyr]*g);  
      }
      
      z[nlyr-1]=0;

      for (ilyr=nlyr-2; ilyr>-1; ilyr--) {
	z[ilyr]=z[ilyr+1]+deltaz[ilyr];
      }

      printf ("time: %d\n", (int)(timesteps*deltat));

      pladv(2);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( TMIN, TMAX, 0, 25000 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", (TMAX-TMIN)/4.0, 0, "bcnst", 5000.0, 0 );
      pllab ("temperature [K]", "z [m]", "");  /* axis labels     */

      plcol0 (9);                         /* color blue  */
      plline (nlyr, T, z);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */


      for (ilyr=0; ilyr<nlyr; ilyr++) {
	printf("z[%d]=%f,  plyr[%d]=%f, T=\n", ilyr, z[ilyr], ilyr, plyr[ilyr], T[ilyr]);
      }

      //sleep(1);
    

  }                                       /* End of Loop */

//  printf("\nTime %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
//for (ilyr=0; ilyr<nlyr; ilyr++) {
//printf("%d %f\n", ilyr, T[ilyr]);
//}

sleep(30);


 return 0;
}
