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

#define TSURF 288.0    /* K */

#define TIME_MAX 3600 * 24 * 360  /* seconds */

#define WAVELENGTH_MAX  1.0E-4  /* meters */
#define WAVELENGTH_STEP 1.0E-6  /* meters */

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
  double deltat=1000;  /* s */
  double Ra=287; /* J/kg K */
  double H=0; /* [Esol] */
  double Tsurf = TSURF;
 // double tau = 5.0;

  int timesteps = 0;
  int ilev=0;
  int ilyr=0;
  int nlyr=10;  /* Number of Layers */
  int nlev=nlyr+1;
  int instabil=FALSE;

  int iwvl=0;
  int nwvl=100;

  double dlambda=1.0E-6;  
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
  double *euptmp=calloc(nlev, sizeof(double));
  double *edntmp=calloc(nlev, sizeof(double));
  double *enet=calloc(nlyr, sizeof(double));

  double *deltaT=calloc(nlyr, sizeof(double));
  double *deltatau=calloc(nwvl,sizeof(double));
  double *lambda=calloc(nwvl,sizeof(double));


  plscolbg (255, 255, 255);   /* background color white */
  plscol0  (15, 0, 0, 0);     /* set color 15 to black  */
  plsdev ("xwin"); /* if not called, the user is asked! */
  plinit ();
  
  plssub(2,2);


  for(ilev=0; ilev<nlev; ilev++) {                              /* Calculation of the Pressure at the Levels p[ilev] */
    p[ilev]=p0*ilev/(nlev-1);
  }
    //printf("%d %f\n", ilev, p[ilev]);
    
    for(ilyr=0; ilyr<nlyr; ilyr++) {                            /* Calculation of the Pressure in the Layers*/
      plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
      //printf("%d %f\n", ilyr, plyr[ilyr]);
    }
   
  
  
  for(ilyr=0; ilyr<nlyr; ilyr++) {                             /* Temperature for all Layers 255K */
    T[ilyr]=Tsurf*pow((plyr[ilyr]/p0),kappa);
  }

  for (ilyr=0;ilyr<nlyr;ilyr++) {
    T[ilyr]=0;
  }
 
  for(ilyr=0;ilyr<nlyr; ilyr++){
    printf ("T=%f\n", T[ilyr]);
  }

  /* for (ilyr=0; ilyr<nlyr; ilyr++){ */
  /*   deltatau[ilyr]=0.1; */
  /*   printf ("tau=%f\n", deltatau[ilyr]); */
  /* } */

 

  lambda[0]=0;

  for (iwvl=0; iwvl<nwvl; iwvl++) {
    lambda[iwvl]= iwvl*dlambda;
    // printf("iwvl= %d, lambda= %g\n", iwvl, lambda[iwvl]);
  }

  for (ilyr=0; ilyr<nlyr; ilyr++) {

  for (iwvl=0; iwvl<8; iwvl++) {
    deltatau[ilyr][iwvl]=10.0/8.0;
  }

  for (iwvl=8; iwvl<13; iwvl++) {
    deltatau[ilyr][iwvl]=0.5/(13-8);
  }

  for (iwvl=13; iwvl<nwvl; iwvl++) {
    deltatau[ilyr][iwvl]=5.0/(nwvl-13);
  }

  for (iwvl=0; iwvl<nwvl; iwvl++) {
    printf("iwvl= %d, deltatau= %g, lambda= %g\n", iwvl, deltatau[iwvl], lambda[iwvl]);
  }
  }
      while (timesteps*deltat<TIME_MAX) {                                       /* Loop limited to 400K */
	// printf("\nNew time %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
	timesteps++;
  
	for(ilev=0; ilev<nlev; ilev++) {
	  edn[ilev] = 0;
	  eup[ilev] = 0;
	}

    //schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup)
    for(iwvl=0; iwvl<nwvl; iwvl++) {
      schwarzschild(deltatau, T, nlev, Tsurf, edntmp, euptmp, dlambda);
      //printf("dlambda %e: edn[4] = %e\n", dlambda, edntmp[4]*WAVELENGTH_STEP);
      for(ilev=0; ilev<nlev; ilev++) {
	edn[ilev] += edntmp[ilev]*WAVELENGTH_STEP;
	eup[ilev] += euptmp[ilev]*WAVELENGTH_STEP;
      }
    }
   
    /* E-netto */
    for(ilyr=0; ilyr<nlyr; ilyr++) {
      enet[ilyr] = eup[ilyr+1] + edn[ilyr] - eup[ilyr] - edn[ilyr+1];
      //printf("enet[%d] = %e - %e = %e\n", ilyr, eup[ilyr+1]+edn[ilyr], eup[ilyr] + edn[ilyr+1], enet[ilyr]);
    }

    /* Temperature gain for each layer */
    for (ilyr=0; ilyr<nlyr; ilyr++)  {
      deltaT[ilyr]=(enet[ilyr]*g*deltat)/((deltap*100.0)*cp);
      T[ilyr] +=deltaT[ilyr];
      //printf("dT[ilyr] = %f\n", ilyr, deltaT[ilyr]);
    }
  
   
    H=Esol-boltzmann(Tsurf)+edn[nlyr-1];
    deltaTsurf=(H*g*deltat)/((deltap*100.0)*cp);                    /* Equation for Temperature gain of the bottom Layer */
    Tsurf += deltaTsurf;
    // printf("%f\n", deltaT); */
  
    
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

    if(timesteps % (int)(50) == 0) {
      printf ("timestep %d\n", timesteps);

      /* Printing time in readable format
       * In order to remove days, hours or minutes
       * just add // in front of the appropriate line
       */
      int time = timesteps*deltat;
      printf("time ");
      printf("%dd = ", (int)(time / 3600 / 24));
      printf("%dh = ", (int)(time / 3600 ));
      printf("%dmin = ", (int)(time / 60 ));
      printf("%ds\n", time);

      printf(" Tsurf=%f\n", Tsurf);


      /* Plot T against p */

      pladv(1);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( 0, 400, p0, 0 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", 100, 0, "bcnst", 150.0, 0 );
      pllab ("temperature [K]", "p [hPa]", "");  /* axis labels     */

      plcol0 (14);                         /* color blue  */
      plline (nlyr, T, plyr);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */

      for (ilyr=0;ilyr<nlyr; ilyr++) {
	deltaz[ilyr]=(Ra*T[ilyr]*deltap)/(plyr[ilyr]*g);  
      }
      
      z[nlyr-1]=0;

      for (ilyr=nlyr-1; ilyr>-1; ilyr--) {
	z[ilyr]=z[ilyr+1]+deltaz[ilyr];
      }


      /* Plot T against z */

      pladv(3);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( 0, 400, 0, 24000 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", 100, 0, "bcnst", 4000.0, 0 );
      pllab ("temperature [K]", "z [m]", "");  /* axis labels     */

      plcol0 (14);                         /* color blue  */
      plline (nlyr, T, z);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */

      /* Plot Heating rate against p */

      pladv(2);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( -0.2, 0.2, p0, 0 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", 100, 0, "bcnst", 4000.0, 0 );
      pllab ("Heating Rate [T/timestep]", "p [hPa]", "");  /* axis labels     */

      plcol0 (12);                         /* color blue  */
      plline (nlyr, deltaT, plyr);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */


      /* Plot Heating rate against z */

      pladv(4);     /* select subpage 1  */
      plvsta();     /* standard viewport */
      plclear();    /* clear subpage     */
      plcol0 (15);  /* color black       */

      plwind( -0.2, 0.2, 0, 24000 );  /* xmin, xmax, ymin, ymax */
      plbox( "bcnst", 100, 0, "bcnst", 4000.0, 0 );
      pllab ("Heating Rate [T/timestep]", "z [m]", "");  /* axis labels     */

      plcol0 (12);                         /* color blue  */
      plline (nlyr, deltaT, z);  /* plot temperature profile  */

      plcol0 (15);                        /* color black */


      for (ilyr=0; ilyr<nlyr; ilyr++) {
        printf("ilyr %d, z=%f,  plyr=%f,theta=%f, T=%f\n", ilyr, z[ilyr], plyr[ilyr],theta[ilyr], T[ilyr]);
      }

      for (ilev=0; ilev<nlev; ilev++) {
	printf("p%f, edn%f, eup%f, enet%f\n", p[ilev], edn[ilev], eup[ilev], enet[ilev]);
      }
    }
    
  }                                       /* End of Loop */

  printf("\nTime %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
  for (ilyr=0; ilyr<nlyr; ilyr++) {
    printf("%d %f\n", ilyr, T[ilyr]);
  }

  printf("\nTsurf=%f\n", Tsurf);

  sleep(30);


  return 0;
    }
