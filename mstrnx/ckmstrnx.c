#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "ascii.h"
#include "ckmstrnx.h"
#include "hydrostat.h"

void ck_(int *nln, double *zb, double *pb, double *tb, int *lyrflag,
	 double *h2o, double *co2, double *o3, double *n2o, double *co,
	 double *ch4, double *o2,
	 double *wvn, double *wgt, double *dtau);

int ck_mstrnx (double *zb, double *pb, double *tb, double *h2o, double *o3,
	       int nlev, int lyrflag,
	       double c_co2, double c_n2o, double c_co, double c_ch4,
	       double c_o2,
	       double ****dtau, double ***wgt, double **wavelength,
	       int *nbnd, int **nch)
{
  static int first=1;
  int ilyr=0, ilev=0, iv=0, iq=0, j=0, status=0;

  double *co2_SI=NULL, *n2o_SI=NULL, *co_SI=NULL, *ch4_SI=NULL, *o2_SI=NULL;
  double *h2o_SI=NULL, *o3_SI=NULL;
  double *h2o_l=NULL, *co2_l=NULL, *o3_l=NULL, *n2o_l=NULL, *co_l=NULL, *ch4_l=NULL, *o2_l=NULL;
  double *denscol=NULL;

  double taumax=100.0; /* maximum optical thickness per layer */

  double *xtmp=NULL, *ytmp=NULL;
  int ntmp=0;

  double *wgt_fortran=NULL, *dtau_fortran=NULL;

  double *wvn=NULL;

  char filename[FILENAME_MAX]="";

  int nlyr = nlev-1;
  
  /* parameters to be used for Fortran and C */
  #include "param.h"

  wvn  = calloc (kbnd+1, sizeof(double));
  

  if (first) {
    if (lyrflag)
      fprintf (stderr, " ... interpreting pressure, temperature, and densities as layer quantities\n");
    else 
      fprintf (stderr, " ... interpreting pressure, temperature, and densities as level quantities\n");

    if (nlev>kln+1) {
      fprintf (stderr, "Error, increase kln in param.h and param.inc to at least %d\n", nlev-1);
      return -1;
    }

    *nbnd   =  kbnd;   

    strcpy(filename, "mstrnx.data/nch.29");
    if ((status=read_2c_file (filename, &xtmp, &ytmp, &ntmp))!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
  
    if (ntmp!=*nbnd) {
      fprintf (stderr, "Error, expected %d entries in %s, found %d.\n", *nbnd, filename, ntmp); 
      return -1;
    }
  
    *nch = calloc (*nbnd, sizeof(int));
    *wavelength  = calloc (*nbnd+1, sizeof(double));

    for (iv=0; iv<*nbnd; iv++) 
      (*nch)[*nbnd-1-iv] = (int) (ytmp[iv]+0.5);

    free (xtmp); free (ytmp);
  }
  
  if (lyrflag) { /* interpret pressure, temperature, and densities as layer quantities */
    h2o_l = calloc (nlev-1, sizeof(double));
    co2_l = calloc (nlev-1, sizeof(double));
    o3_l  = calloc (nlev-1, sizeof(double));
    n2o_l = calloc (nlev-1, sizeof(double));
    co_l  = calloc (nlev-1, sizeof(double));
    ch4_l = calloc (nlev-1, sizeof(double));
    o2_l  = calloc (nlev-1, sizeof(double));
    
    for (ilev=0; ilev<nlev-1; ilev++)  {
      h2o_l[ilev] = h2o[ilev];
      o3_l [ilev] = o3 [ilev];
      co2_l[ilev] = c_co2;  
      n2o_l[ilev] = c_n2o;
      co_l [ilev] = c_co;
      ch4_l[ilev] = c_ch4;
      o2_l [ilev] = c_o2;
    }
  }
  else {         /* interpret pressure, temperature, and densities as level quantities */
#if 0
    /* allocate memory for trace gas profiles */
    h2o_SI = calloc (nlev, sizeof(double));
    o3_SI  = calloc (nlev, sizeof(double));
    co2_SI = calloc (nlev, sizeof(double));
    n2o_SI = calloc (nlev, sizeof(double));
    co_SI  = calloc (nlev, sizeof(double));
    ch4_SI = calloc (nlev, sizeof(double));
    o2_SI  = calloc (nlev, sizeof(double));

    /* convert from ppm to mol/mol */
    for (ilev=0; ilev<nlev; ilev++)  {
      h2o_SI[ilev] = h2o[ilev]*1E-6;
      o3_SI [ilev] = o3 [ilev]*1E-6;
      co2_SI[ilev] = c_co2    *1E-6;  
      n2o_SI[ilev] = c_n2o    *1E-6;
      co_SI [ilev] = c_co     *1E-6;
      ch4_SI[ilev] = c_ch4    *1E-6;
      o2_SI [ilev] = c_o2     *1E-6;
    }
    
    /* now convert level densities to layer-averaged quantities assuming         */
    /* hydrostatic pressure, and linear temperature and concentration variations */
    
    calc_denscol_hydrostatic_profile (nlev, pb, tb, h2o_SI, 0, &h2o_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb, co2_SI, 0, &co2_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb,  o3_SI, 0, &o3_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb, n2o_SI, 0, &n2o_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb,  co_SI, 0, &co_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb, ch4_SI, 0, &ch4_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb,  o2_SI, 0, &o2_l);
    calc_denscol_hydrostatic_profile (nlev, pb, tb,   NULL, 0, &denscol);    /* air */
    
    
    /* convert to layer-averaged concentration [ppm] by dividing by the air column */
    for (ilyr=0;ilyr<nlyr;ilyr++)  {
      h2o_l [ilyr] *= 1E6 / denscol[ilyr];
      co2_l [ilyr] *= 1E6 / denscol[ilyr];
      o3_l  [ilyr] *= 1E6 / denscol[ilyr];
      n2o_l [ilyr] *= 1E6 / denscol[ilyr];
      co_l  [ilyr] *= 1E6 / denscol[ilyr];
      ch4_l [ilyr] *= 1E6 / denscol[ilyr];
      o2_l  [ilyr] *= 1E6 / denscol[ilyr];
    }
#else
    fprintf(stderr,"Error, this mode not implemented\n");
    return -1;
#endif
  }
    
//  fprintf (stderr, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", h2o_l[nlyr-1], co2_l[nlyr-1], o3_l[nlyr-1], n2o_l[nlyr-1], co_l[nlyr-1], ch4_l[nlyr-1], o2_l[nlyr-1]);

  wgt_fortran  = calloc (*nbnd*kch, sizeof(double));
  dtau_fortran = calloc (*nbnd*kch*kln, sizeof(double));

  /* get correlated-k dtau's and weights */
  ck_(&nlyr,zb,pb,tb,&lyrflag,h2o_l,co2_l,o3_l,n2o_l,co_l,ch4_l,o2_l,
      wvn,wgt_fortran,dtau_fortran);

  /* convert wavenumbers [cm-1] to wavelengths [nm] */
  for (iv=0; iv<=*nbnd; iv++) 
    (*wavelength)[*nbnd-iv] = 1.0E7 / wvn[iv];

  if (first) {
    *wgt = calloc (*nbnd, sizeof(double *));
    for (iv=0;iv<*nbnd;iv++)
      (*wgt)[iv] = calloc ((*nch)[iv], sizeof(double));

    *dtau = calloc (*nbnd, sizeof(double **));
    for (iv=0;iv<*nbnd;iv++) {
      (*dtau)[iv] = calloc ((*nch)[iv], sizeof(double *));
      for (iq=0;iq<(*nch)[iv];iq++)
	(*dtau)[iv][iq] = calloc (nlyr, sizeof(double));
    }
  }

  /* convert 1D Fortran array to 2D and 3D C arrays */
  j=0;
  for (iv=0; iv<*nbnd; iv++)
    for (iq=0; iq<kch; iq++) {
      if (iq<(*nch)[*nbnd-1-iv])
	(*wgt)[*nbnd-1-iv][iq] = wgt_fortran[j];
      j++;
    }

  j=0;
  for (iv=0;iv<*nbnd;iv++)
    for (iq=0;iq<kch;iq++)
      for (ilyr=0;ilyr<kln;ilyr++) {
	if (iq<(*nch)[*nbnd-1-iv] && ilyr<nlyr)
	  (*dtau)[*nbnd-1-iv][iq][ilyr] = dtau_fortran[j];
	j++;
      }

  /* restrict layer optical thickness to taumax, */
  /* otherwise the twostream solution is likely  */
  /* to become unstable; the solution did not    */
  /* change significantly when taumax was varied */
  /* between 10 and 300 (above we obtained NANs) */
  for (iv=0;iv<*nbnd;iv++)
    for (iq=0;iq<(*nch)[iv];iq++)
      for (ilyr=0;ilyr<nlyr;ilyr++)
	if ((*dtau)[iv][iq][ilyr]>taumax)
	  (*dtau)[iv][iq][ilyr]=taumax;
  
  free (wvn);
  free (wgt_fortran);
  free (dtau_fortran);

  free(h2o_SI);
  free(o3_SI);
  free(co2_SI);
  free(n2o_SI);
  free(co_SI);
  free(ch4_SI);
  free(o2_SI);

  free(denscol);

  free(h2o_l);
  free(co2_l);
  free(o3_l);
  free(n2o_l);
  free(co_l);
  free(ch4_l);
  free(o2_l);

  first=0;

  return 0;
}


void ck_mstrnx_fortran_ (double *zb, double *pb, double *tb, double *h2o, double *o3, int *nlev, int *lyrflag, 
			 double *c_co2, double *c_n2o, double *c_co, double *c_ch4, double *c_o2,
			 double *dtau, double *wgt, double *wavelength, int *nbnd, int *nch, int *status)
{
  double ***dtau_d;
  double **wgt_d;  
  double *wavelength_d;  
  int *nch_d;  
  int iv=0, iq=0,j=0,ilyr=0;

  *status = ck_mstrnx (zb, pb, tb, h2o, o3, *nlev, *lyrflag,
		       *c_co2, *c_n2o, *c_co, *c_ch4, *c_o2,
		       &dtau_d, &wgt_d, &wavelength_d, nbnd, &nch_d);
  if (*status!=0)
    fprintf(stderr, "Error %d returned by ck_mstrnx()\n", *status);

  for (iv=0; iv<*nbnd; iv++)
    nch[iv] = nch_d[iv];

  for (iv=0; iv<*nbnd; iv++)
    wavelength[iv] = wavelength_d[iv];

  j=0;
  for (iv=0; iv<*nbnd; iv++)
    for (iq=0; iq<10; iq++) {
      if (iq<nch[iv])
	wgt[j] = wgt_d[iv][iq];
      j++;
    }

  j=0;
  for (iv=0 ; iv<*nbnd; iv++)
    for (iq=0; iq<10; iq++)
      for (ilyr=0;ilyr<*nlev-1;ilyr++) {
	if (iq<nch[iv])
	  dtau[j] = dtau_d[iv][iq][ilyr];
	j++;
      }

}
