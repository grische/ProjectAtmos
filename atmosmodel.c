#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <plplot/plplot.h>
#include "schwarz.h"
#include "mstrnx/ascii.h"

#define FALSE 0
#define TRUE 1

#define R 8.31447     /* J/(mol K) */
#define M 0.0289644   /* kg / mol */
#define g 9.8065      /* m/s^2 */

#define TSURF 100.0    /* K */
#define PSURF 1000     /* hPa */

#define TIME_MAX 3600 * 24 * 360  /* seconds */

#define WAVELENGTH_STEP_ROUGH 1.0E-8  /* meters */
#define WAVELENGTH_INC_STEP   1     /* read every nth value */

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
#ifndef _NOPLOT
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

  plwind( 0, 400, 0, 40000 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 100, 0, "bcnst", 5000.0, 0 );
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
  plbox( "bcnst", 2, 0, "bcnst", 150.0, 0 );
  pllab ("Heating Rate [T/day]", "p [hPa]", "");  /* axis labels     */

  plcol0 (12);                         /* color blue  */
  plline (nlyr, deltaTday, plyr);  /* plot temperature profile  */

  plcol0 (15);                        /* color black */


  /* Plot Heating rate against z */

  pladv(4);     /* select subpage 1  */
  plvsta();     /* standard viewport */
  plclear();    /* clear subpage     */
  plcol0 (15);  /* color black       */


  plwind( -20, 20, 0, 40000 );  /* xmin, xmax, ymin, ymax */
  plbox( "bcnst", 2, 0, "bcnst", 5000.0, 0 );
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
#endif
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
  int nwvl_small;

  double deltap=PSURF/nlyr;
  double kappa=Ra/cp;
  
  double *p=calloc(nlev,sizeof(double));
  double *T=calloc(nlyr,sizeof(double));
  double *theta=calloc(nlev,sizeof(double));
  double *plyr=calloc(nlyr,sizeof(double));
  double *z=calloc(nlyr,sizeof(double));
  double *deltaz=calloc(nlyr,sizeof(double));

  double *tmplev=calloc(nlev, sizeof(double)); /* temporary vector with length lev */
  double *tmplyr=calloc(nlyr, sizeof(double)); /* temporary vector with length lyr */

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
  int co2_nwvl;
  {
    /* the variables initialized inside brackets are only valid until the bracket closes */
    double *wvn;
    int status=0;
    int co2_nlyr;

    /* read file */
    if ( (status = ASCII_file2xy2D ("co2.dtau", &co2_nwvl, &co2_nlyr, &wvn, &deltatau_co2)) != 0 ) {
      printf("Error %d reading file, quitting!\n", status);
      return EXIT_FAILURE;
    }
    else {
      if ( co2_nlyr != nlyr) {
        printf("Error, file has different amount of columns (%d) than nylr (%d).\n", co2_nlyr, nlyr);
        return EXIT_FAILURE;
      }
      printf("co2.dtau: found %d lines in file\n", co2_nwvl);
    }

    lambda_co2 = calloc(co2_nwvl, sizeof(double));
    for(iwvl=0; iwvl<co2_nwvl; iwvl++) {
      lambda_co2[co2_nwvl-iwvl-1] = 1.0/(100.0*wvn[iwvl]);
    }
    printf("CO2 Profile: min=%e, max=%e\n", lambda_co2[0], lambda_co2[co2_nwvl-1]);
    free(wvn);

    printf("Flipping deltatau_co2\n");
    double* temp = calloc(nlyr, sizeof(double));
    for (iwvl=0; iwvl<co2_nwvl; iwvl++) {
        for(ilyr=0; ilyr<nlyr; ilyr++) 
           temp[nlyr-1-ilyr] = deltatau_co2[iwvl][ilyr];
        for(ilyr=0; ilyr<nlyr; ilyr++)
           deltatau_co2[iwvl][ilyr] = temp[ilyr];
    }
    free(temp);

    /* nwvl goes only from: WAVELENGTH_STEP_ROUGH, ..., lambda_co2 */
    nwvl_small=(lambda_co2[0]/WAVELENGTH_STEP_ROUGH);
    nwvl = nwvl_small+((1000.0E-6 - lambda_co2[co2_nwvl-1])/WAVELENGTH_STEP_ROUGH);
    printf("nwvl = %e / %e = %d\n", lambda_co2[0], WAVELENGTH_STEP_ROUGH, nwvl);

    /* initialize tau matrix non-co2 iwvl */
    deltatau_rough = calloc(nwvl,sizeof(double*));
    for(iwvl=0; iwvl<nwvl; iwvl++) {
      deltatau_rough[iwvl] = calloc(nlyr,sizeof(double));
    }
   
    /* initialize lambda for each non-co2 iwvl */
    lambda_rough = calloc(nwvl,sizeof(double));
    for (iwvl=0; iwvl<nwvl_small; iwvl++) {
      lambda_rough[iwvl]= (iwvl+1)*WAVELENGTH_STEP_ROUGH;
    }
    for (iwvl=nwvl_small; iwvl<nwvl; iwvl++) {
      lambda_rough[iwvl]=(20E-6 - lambda_rough[nwvl_small-1])+((iwvl+1)*WAVELENGTH_STEP_ROUGH);
    }

    for (ilyr=0; ilyr<nlyr; ilyr++) {
      for (iwvl=0; iwvl<nwvl; iwvl++) {
        if(lambda_rough[iwvl]<8.0E-6) {
          /* deltatau = 10.0/nlyr  for  0 ... 80 nm */
          deltatau_rough[iwvl][ilyr]=10.0/nlyr;
        } 
        if(lambda_rough[iwvl] >= 8.0E-6 && lambda_rough[iwvl] < lambda_co2[0]) {
          /* deltatau =  0.5/nlyr  for  80 nm ... lambda_co2_min */
          deltatau_rough[iwvl][ilyr]=0.5/nlyr;
        }
	
	if(lambda_rough[iwvl] > lambda_co2[co2_nwvl-1]) {
	  deltatau_rough[iwvl][ilyr]=10.0/nlyr;
	}
      }
      
    }

  }

#ifndef _NOPLOT
  plscolbg (255, 255, 255);   /* background color white */
  plscol0  (15, 0, 0, 0);     /* set color 15 to black  */
  plsdev ("xwin"); /* if not called, the user is asked! */
  plinit ();
  
  plssub(2,2);
#endif

  for(ilev=0; ilev<nlev; ilev++) {                              /* Calculation of the Pressure at the Levels p[ilev] */
    p[ilev]=PSURF*ilev/(nlev-1);
  }
  //printf("%d %f\n", ilev, p[ilev]);
    
  for(ilyr=0; ilyr<nlyr; ilyr++) {                            /* Calculation of the Pressure in the Layers*/
    plyr[ilyr]= 0.5*(p[ilyr]+p[ilyr+1]);
    //printf("%d %f\n", ilyr, plyr[ilyr]);
  }
   
  
  
  {
    double *Tlev;
    double *temp;
    int status;
    int ntemplev;
   
    
    status = read_3c_file("fpda.atm", &temp, &temp, &Tlev, &ntemplev);
    if (status !=0) {
      printf("Error reading Temperature profile\n");
      return EXIT_FAILURE;
    }
    else if (ntemplev != nlev) {
      printf("Error! nlev fault!!!!!\n");
      return EXIT_FAILURE;
    }

    for (ilev=0; ilev<nlev; ilev++) {
      printf("ilev = %d, Tlev = %f\n", ilev, Tlev[ilev]);
    }


    for (ilyr=0;ilyr<nlyr;ilyr++) {
      T[ilyr]=(Tlev[ilyr]+Tlev[ilyr+1])/2;
      printf("ilyr = %d, T = %f\n", ilyr, T[ilyr]);
    }

    Tsurf=Tlev[nlev-1];
    
    printf("Tsurf = %f\n", Tsurf);
    
    free(temp);
    free(Tlev);

  }

  while (timesteps*deltat<TIME_MAX) {                                       /* Loop limited to 400K */
    // printf("\nNew time %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
    timesteps++;
  
    for(ilev=0; ilev<nlev; ilev++) {
      edn[ilev] = 0;
      eup[ilev] = 0;
      if(ilev < nlyr)
        enet[ilev]=0;
    }

    //schwarzschild(const double tau, const double *T, const int nlev, const double Ts, double *edn, double *eup)
    /* Non-CO2 absorption */
    for(iwvl=0; iwvl<nwvl; iwvl++) {
      schwarzschild(deltatau_rough[iwvl], T, nlev, Tsurf, edntmp, euptmp, lambda_rough[iwvl], tmplev, tmplyr);
      //printf("dlambda %e: edn[4] = %e\n", dlambda, edntmp[4]*WAVELENGTH_STEP);
      for(ilev=0; ilev<nlev; ilev++) {
	edn[ilev] += edntmp[ilev]*WAVELENGTH_STEP_ROUGH;
	eup[ilev] += euptmp[ilev]*WAVELENGTH_STEP_ROUGH;
      }
    }
    /* between non-co2 and co2 */
    schwarzschild(deltatau_rough[nwvl_small-1], T, nlev, Tsurf, edntmp, euptmp, (lambda_co2[0] + lambda_rough[nwvl_small-1])/2.0, tmplev, tmplyr);
    for(ilev=0; ilev<nlev; ilev++) {
      edn[ilev] += edntmp[ilev]*(lambda_co2[0]-lambda_rough[nwvl_small-1]);
      eup[ilev] += euptmp[ilev]*(lambda_co2[0]-lambda_rough[nwvl_small-1]);
    }

   
    /* CO2 absorption */
    for(iwvl=0; iwvl<co2_nwvl-WAVELENGTH_INC_STEP-1; iwvl+=WAVELENGTH_INC_STEP) {
      schwarzschild(deltatau_co2[iwvl], T, nlev, Tsurf, edntmp, euptmp, lambda_co2[iwvl], tmplev, tmplyr);
      //printf("dlambda %e: edn[4] = %e\n", dlambda, edntmp[4]*WAVELENGTH_STEP);
      const double dstep = lambda_co2[iwvl+WAVELENGTH_INC_STEP]-lambda_co2[iwvl];
      for(ilev=0; ilev<nlev; ilev++) {
	edn[ilev] += edntmp[ilev]*dstep;
	eup[ilev] += euptmp[ilev]*dstep;
      }
    }  
    /* E-netto */
    for(ilyr=0; ilyr<nlyr-1; ilyr++) {
      enet[ilyr] = eup[ilyr+1] + edn[ilyr] - eup[ilyr] - edn[ilyr+1];
     
    }

    for (ilev=0; ilev<nlev; ilev++) {
      printf("ilev = %d, eup = %f, edn = %f\n", ilev, eup[ilev], edn[ilev]);
    }

    /* Temperature gain for each layer */
    for (ilyr=0; ilyr<nlyr-1; ilyr++)  {
      deltaT[ilyr]=(enet[ilyr]*g*deltat)/((deltap*100.0)*cp);
      T[ilyr] +=deltaT[ilyr];
      //printf("dT[ilyr] = %f\n", ilyr, deltaT[ilyr]);
    }
  
    
    /* Surface Temperature Gain */
    deltaT[nlyr-1]=((Esol+edn[nlyr-1]-eup[nlyr-1])*g*deltat)/((deltap*100.0)*cp);
    T[nlyr-1] += deltaT[nlyr-1];


    Tsurf=T[nlyr-1];

  
    
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

    if(timesteps % (int)(5) == 0) {
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
      for (ilyr=nlyr-2; ilyr >= 0; ilyr--) {
	z[ilyr]=z[ilyr+1]+deltaz[ilyr];
      }

      for (ilyr=0; ilyr<nlyr; ilyr++) {
	deltaTday[ilyr]=deltaT[ilyr]*86.4/2;
      }

      plotall(nlyr, T, plyr, z, deltaTday);

      /* for (ilyr=0; ilyr<nlyr; ilyr++) { */
      /*   printf("ilyr %d, z=%f,  plyr=%f,theta=%f, T=%f\n", ilyr, z[ilyr], plyr[ilyr],theta[ilyr], T[ilyr]); */
      /* } */

      for (ilev=0; ilev<nlev; ilev++) {
      	printf("p%f, edn%f, eup%f, enet%f\n", p[ilev], edn[ilev], eup[ilev], enet[ilev]);
      }
    }
    
  }                                       /* End of Loop */

  free(tmplev);
  free(tmplyr);

  printf("\nTime %d: T = %f\n", (int)(timesteps*deltat), T[nlyr-1]);
  for (ilyr=0; ilyr<nlyr; ilyr++) {
    printf("%d %f\n", ilyr, T[ilyr]);
  }

  printf("\nTsurf=%f\n", Tsurf);

  sleep(30);

  return 0;
}
