/* @(#)ckmstrnx.h
 */

#ifndef _ckmstrnx_h
#define _ckmstrnx_h 1


#if defined (__cplusplus)
extern "C" {
#endif

int ck_mstrnx (double *zb, double *pb, double *tb, double *h2o, double *o3,
	       int nlev, int lyrflag,
	       double c_co2, double c_n2o, double c_co, double c_ch4,
	       double c_o2,
	       double ****dtau, double ***wgt, double **wavelength,
	       int *nbnd, int **nch);

#if defined (__cplusplus)
}
#endif

#endif

