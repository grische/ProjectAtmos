/*--------------------------------------------------------------------
 * ask Robert Buras
 * ulrike 16.06.2010
 *--------------------------------------------------------------------*/

#ifndef __bandec_h
#define __bandec_h

#if defined (__cplusplus)
extern "C" {
#endif

void bandec(double **a, long n, double **al,
	    long indx[], double *d);

void banbks(double **a, long n, double **al,
	    long indx[], double b[]);


#endif
