#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <plplot/plplot.h>
#include "schwarz.h"
#include "ascii.h"

#define FALSE 0
#define TRUE 1

#define R 8.31447     /* J/(mol K) */
#define M 0.0289644   /* kg / mol */
#define g 9.8065      /* m/s^2 */

#define TSURF 100.0    /* K */
#define PSURF 1000     /* hPa */

#define TIME_MAX 3600 * 24 * 360  /* seconds */

#define WAVELENGTH_STEP_ROUGH 1.0E-7  /* meters */

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


int plotall(int nlyr, double* T, double* plyr, double* z, double* deltaTday) {
  /* Plot T against p */

  pladv(1);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */

  plwind( 0, 400, PSURF, 0 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 100, 0, "bcnst", 150.0, 0 );
  pllab ("temperature [K]", "p [hPa]", "");  /* axis labels     */

  plcol0 (14);                         /* color blue  */
  plline (nlyr, T, plyr);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */

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

  plwind( -20, 20, PSURF, 0 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 2, 0, "bcnst", 4000.0, 0 );
  pllab ("Heating Rate [T/day]", "p [hPa]", "");  /* axis labels     */

  plcol0 (12);                         /* color blue  */
  plline (nlyr, deltaTday, plyr);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */


  /* Plot Heating rate against z */

  pladv(4);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */


  plwind( -20, 20, 0, 24000 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 2, 0, "bcnst", 4000.0, 0 );
  pllab ("Heating Rate [T/day]", "z [m]", "");  /* axis labels     */

  plcol0 (12);                         /* color blue  */
  plline (nlyr, deltaTday, z);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */


  /* Tsurftime[timesteps]=T[nlyr-1]; */
  /* day[timesteps]=timesteps/86.4; */

  /* /\* Plot T against timestep *\/ */

  /* pladv(3);     /\* select subpage 1  *\/ */
  /* plvsta();     /\* standard viewport *\/ */
  /* plclear();    /\* clear subpage     *\/ */
  /* plcol0 (15);  /\* color black       *\/ */

  /* plwind( 0, 315360, 0, 400 );  /\* xmin, xmax, ymin, ymax *\/ */
  /* plbox( "bcnst", 100000000000, 0, "bcnst", 4000.0, 0 ); */
  /* pllab ("Heating Rate [T/day]", "z [m]", "");  /\* axis labels     *\/ */

  /* plcol0 (12);                         /\* color blue  *\/ */
  /* plline (timesteps, day, Tsurftime);  /\* plot temperature profile  *\/ */

  /* plcol0 (15);                        /\* color black *\/ */
}


int main() {                                                     /* Definition of Variables and Quantities  */
  
  double cp=1004;  /* J/kg K */
  double Esol=235;  /* W/m^2 = J*/
  double deltaTsurf=0; /* K */
  double deltat=1000;  /* s */
  double Ra=287; /* J/kg K */
  double H=0; /* [Esol] */
  double Tsurf = TSURF;
  // double tau = 5.0;

  int timesteps = 0;
  int ilev=0;
  int ilyr=0;
  int nlyr=20;  /* Number of Layers */
  int nlev=nlyr+1;
  int instabil=FALSE;

  int iwvl=0;
  int nwvl;

  double deltap=PSURF/nlyr;
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

  double *deltaTday=calloc(nlyr, sizeof(double));
  double *deltaT=calloc(nlyr, sizeof(double));
  double *lambda_rough;
  double *lambda_co2;
  double** deltatau_rough;
  double** deltatau_co2;

  /* initial values from file */
  int co2_lines;
  {
    /* the variables initialized inside brackets are only valid until the bracket closes */
    double *wvn;
    int status=0;
    int co2_nlyr;

    /* read file */
    if ( (status = ASCII_file2xy2D ("co2.dtau", &co2_lines, &co2_nlyr, &wvn, &deltatau_co2)) != 0 ) {
      printf("Error %d reading file, quitting!\n", status);
      return EXIT_FAILURE;
    }
    else {
      if ( co2_nlyr != nlyr) {
        printf("Error, file has different amount of columns (%d) than nylr (%d).\n", co2_nlyr, nlyr);
        return EXIT_FAILURE;
      }
      printf("co2.dtau: found %d lines in file\n", co2_lines);
    }

    lambda_co2 = calloc(co2_lines, sizeof(double));
    for(iwvl=0; iwvl<co2_lines; iwvl++) {
      lambda_co2[co2_lines-iwvl-1] = 1.0/(100.0*wvn[iwvl]);
    }
    printf("CO2 Profile: min=%e, max=%e\n", lambda_co2[0], lambda_co2[co2_lines-1]);
    free(wvn);

    /* nwvl goes only from: WAVELENGTH_STEP_ROUGH, ..., lambda_co2 */
    nwvl = lambda_co2[0] / WAVELENGTH_STEP_ROUGH;
    printf("nwvl = %e / %e = %d\n", lambda_co2[0], WAVELENGTH_STEP_ROUGH, nwvl);

    /* initialize tau matrix non-co2 iwvl */
    deltatau_rough = calloc(nwvl,sizeof(double*));
    for(iwvl=0; iwvl<nwvl; iwvl++) {
      deltatau_rough[iwvl] = calloc(nlyr,sizeof(double));
    }
   
    /* initialize lambda for each non-co2 iwvl */
    lambda_rough = calloc(nwvl,sizeof(double));
    for (iwvl=0; iwvl<nwvl; iwvl++) {
      lambda_rough[iwvl]= (iwvl+1)*WAVELENGTH_STEP_ROUGH;
    }

    for (ilyr=0; ilyr<nlyr; ilyr++) {
      for (iwvl=0; iwvl<nwvl; iwvl++) {
        if(lambda_rough[iwvl]<8.0E-6) {
          /* deltatau = 10.0/nlyr  for  0 ... 80 nm */
          deltatau_rough[iwvl][ilyr]=10.0/nlyr;
        } 
        else {
          /* deltatau =  0.5/nlyr  for  80 nm ... lambda_co2_min */
          deltatau_rough[iwvl][ilyr]=0.5/nlyr;
        }
      }
    }

  }


  plscolbg (255, 255, 255);   /* background color white */
  plscol0  (15, 0, 0, 0);     /* set color 15 to black  */
  plsdev ("xwin"); /* if not called, the user is asked! */
  plinit ();
  
  plssub(2,2);


  for(ilev=0; ilev<nlev; ilev++) {                              /* Calculation of the Pressure at the Levels p[ilev] */
    p[ilev]=PSURF*ilev/(nlev-1);
  }
  //printf("%d %f\n", ilev, p[ilev]);
    
  for(ilyr=0; ilyr<nlyr; ilyr++) {                            /* Calculation of the Pressure in the Layers*/
    plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
    //printf("%d %f\n", ilyr, plyr[ilyr]);
  }
   
  
  
  for(ilyr=0; ilyr<nlyr; ilyr++) {                             /* Temperature for all Layers 255K */
    T[ilyr]=Tsurf*pow((plyr[ilyr]/PSURF),kappa);
  }

  for (ilyr=0;ilyr<nlyr;ilyr++) {
    T[ilyr]=100;
  }
 
  for(ilyr=0;ilyr<nlyr; ilyr++){
    printf ("T=%f\n", T[ilyr]);
  }

  /* for (ilyr=0; ilyr<nlyr; ilyr++){ */
  /*   deltatau[ilyr]=0.1; */
  /*   printf ("tau=%f\n", deltatau[ilyr]); */
  /* } */

  while (timesteps*deltat<TIME_MAX) {                                       /* Loop limited to 400K */
    // printf("\nNew time %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
    timesteps++;
  
    for(ilev=0; ilev<nlev; ilev++) {
      edn[ilev] = 0;
      eup[ilev] = 0;
      enet[ilev]=0;
    }

    //schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup)
    /* Non-CO2 absorption */
    for(iwvl=0; iwvl<nwvl; iwvl++) {
      schwarzschild(deltatau_rough[iwvl], T, nlev, Tsurf, edntmp, euptmp, lambda_rough[iwvl]);
      //printf("dlambda %e: edn[4] = %e\n", dlambda, edntmp[4]*WAVELENGTH_STEP);
      for(ilev=0; ilev<nlev; ilev++) {
	edn[ilev] += edntmp[ilev]*WAVELENGTH_STEP_ROUGH;
	eup[ilev] += euptmp[ilev]*WAVELENGTH_STEP_ROUGH;
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
  
    
    /* Surface Temperature Gain */
    deltaTsurf=((Esol+edn[nlev-1]-eup[nlev-1])*g*deltat)/((deltap*100.0)*cp);
    Tsurf += deltaTsurf;

  
    
    for (ilyr=0; ilyr<nlyr; ilyr++) {                           /* Conversion from T to theta */
      theta[ilyr]=pow(PSURF/plyr[ilyr], kappa)*T[ilyr];
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
      T[ilyr] = theta[ilyr] /  (pow(PSURF/plyr[ilyr], kappa));

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

      for (ilyr=0;ilyr<nlyr; ilyr++) {
	deltaz[ilyr]=(Ra*T[ilyr]*deltap)/(plyr[ilyr]*g);  
      }
      
      z[nlyr-1]=0;

      for (ilyr=nlyr-1; ilyr>-1; ilyr--) {
	z[ilyr]=z[ilyr+1]+deltaz[ilyr];
      }

      for (ilyr=0; ilyr<nlyr; ilyr++) {
	deltaTday[ilyr]=deltaT[ilyr]*86.4/2;
      }

      plotall(nlyr, T, plyr, z, deltaTday);

      /* for (ilyr=0; ilyr<nlyr; ilyr++) { */
      /*   printf("ilyr %d, z=%f,  plyr=%f,theta=%f, T=%f\n", ilyr, z[ilyr], plyr[ilyr],theta[ilyr], T[ilyr]); */
      /* } */

      /* for (ilev=0; ilev<nlev; ilev++) { */
      /* 	printf("p%f, edn%f, eup%f, enet%f\n", p[ilev], edn[ilev], eup[ilev], enet[ilev]); */
      /* } */
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
