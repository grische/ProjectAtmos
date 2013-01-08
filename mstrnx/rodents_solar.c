/*--------------------------------------------------------------------
 * $Id$
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
 * RODENTS: ROberts' Delta-EddingtoN Two-Stream
 *
 * by Robert Buras
 *
 * development in joint cooperation with Bernhard Mayer
 *
 * implementation partly by Ulrike Wissmeier
 *
 * Programmed using the book by Zdunkowski, Trautmann and Bott,
 * "Radiation in the Atmosphere", chapters 6.1-6.4
 *
 * Note that in the 2007 print there are misprints in Eqs. 6.50 and 6.88
 * Also, we implemented thermal emission differently than in chapter 6.5
 * of the book.
 *
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "rodents_solar.h"
#include "ascii.h"
#include "bandec.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

static inline int err_out ( char *out, int status );
static void free_input(rodents_input_struct input);
static int setup_factors ( /* INPUT */
			  int     nlyr,
			  double *dtau, 
			  double *omega_0, 
			  double *g, 
			  double *f, 
			  double  S_0,
			  double  mu_0,
			  double  albedo,
			  /* OUTPUT */
			  rodents_input_struct *input);

static int multi_layer_solver ( rodents_input_struct   input,
				rodents_output_struct  output);

static void solve_direct(double *a33, double *a33_u, int nlyr, double S_0,
			 double *S, double *S_u);

static void setup_equations( /* input */
			    rodents_input_struct  input,
			    double               *s_direct,
			    /* output */
			    double              **matr,
			    double               *rhs );

static void cp_to_result(double *b,
			 int     nlyr,
			 double *e_minus,
			 double *e_plus);

static int process_output ( /* input */
			   rodents_input_struct   input,
			   rodents_output_struct  output,
			   /* output */
			   double                *fldn,
			   double                *flup,
			   double                *fldir );


/*************************************************************/
/* rodents, former tsm                                       */
/*************************************************************/

int rodents_solar ( /* INPUT */
		   int     nlyr,
		   double *dtau,
		   double *omega_0,
		   double *g, 
		   double *f,
		   double  S_0,
		   double  mu_0,
		   double  albedo,
		   /* OUTPUT */
		   double *fldn,
		   double *flup,
		   double *fldir )
{
  int status=0;
  rodents_input_struct input;
  rodents_output_struct output;

  output.e_minus       = calloc((size_t) nlyr+1, sizeof(double));
  output.e_plus        = calloc((size_t) nlyr+1, sizeof(double));
  output.s_direct      = calloc((size_t) nlyr+1, sizeof(double));
  output.s_direct_unsc = calloc((size_t) nlyr+1, sizeof(double));
    
  /* initialization */

  status = setup_factors ( nlyr, dtau, omega_0, g, f,
			   S_0, mu_0, albedo,
			   &input);
  if (status)
    return err_out("Error %d returned from setup_factors!\n",status);

  /* solvers */
  status = multi_layer_solver (input, output);
  if (status)
    return err_out("Error %d returned from multi_layer_solver!\n",status);

  status = process_output (input, output, fldn, flup, fldir);
  if (status)
    return err_out("Error %d returned from process_output!\n",status);

  /* free memory */  

  free_input (input);

  free(output.e_minus);
  free(output.e_plus);
  free(output.s_direct);
  free(output.s_direct_unsc);
  
  return 0;
}


/*************************************************************/
/* setup_factors                                             */
/* number of flops per layer: (needs update)                 */
/* add : 14                                                  */
/* mult: 23                                                  */
/* div :  5                                                  */
/* exp :  2                                                  */
/* sqrt:  1                                                  */
/*************************************************************/

static int setup_factors ( /* INPUT */
			  int     nlyr,
			  double *dtau, 
			  double *omega_0, 
			  double *g, 
			  double *f, 
			  double  S_0,
			  double  mu_0,
			  double  albedo,
			  /* OUTPUT */
			  rodents_input_struct *input)
{
  int i=0;
  double bscr=0, b_mmu_0 =0;
  double dtau_d=0.0, g_d=0.0, omega_0_d=0.0;
  double alpha_1=0.0, alpha_2=0.0, alpha_3=0.0, alpha_4=0.0;
  double alpha_5=0.0, alpha_6=0.0;
  double A=0.0, lambda=0.0;
  double mu_0_inv=0.0;
  double den1=0.0, exp1 = 0.0, term1=0.0;

  double epsilon = 1e-6;

  input->nlyr     = nlyr;
  input->S_0      = S_0;
  input->mu_0     = mu_0;
  input->albedo   = albedo;

  /* allocate memory */
  input->a11 = calloc( nlyr, sizeof(double) );
  input->a12 = calloc( nlyr, sizeof(double) );

  input->a13 = calloc( nlyr, sizeof(double) );
  input->a23 = calloc( nlyr, sizeof(double) );
  input->a33 = calloc( nlyr, sizeof(double) );

  /* unscaled extinction for direct radiation */
  input->a33_unsc = calloc( nlyr, sizeof(double) );

  /* choose delta scaling factor using HG */
  /* for (i=0;i<nlyr;i++) */
  /*   f[i] = g[i] * g[i]; */

  /* calculate parameters */

  mu_0_inv = 1./mu_0;

  for (i=0;i<nlyr;i++) {

    /* catch some singularities */ 
    if ( omega_0[i] >= 1.0 )
      omega_0[i] = 1.0 - epsilon;

    if ( omega_0[i] * g[i] == 1.0 )
      omega_0[i] *= 1.0 - epsilon;

    /* do delta scaling */
    dtau_d = dtau[i] * ( 1. - omega_0[i] * f[i] );
    g_d    = ( g[i] - f[i] ) / ( 1. - f[i] );
    omega_0_d = omega_0[i] * ( 1. - f[i] )
      / ( 1. - f[i] * omega_0[i] );

    b_mmu_0 = 0.5 - 0.75 * g_d * mu_0;

    /* Eq. 6.64 */
    bscr = 0.5 - 0.375 * g_d;
    alpha_1 = 2. * ( 1. - omega_0_d * ( 1. - bscr ) ) - 0.25;
    alpha_2 = 2. * omega_0_d * bscr - 0.25;

    /* Eq. 6.69 */
    lambda = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 );

    if ( lambda * dtau_d > 1e2 ) {
      input->a11[i] = 0.0;
      input->a12[i] = ( alpha_1 - lambda ) / alpha_2;
    }
    else { 
      /* exp1, term1 needed for A, a_matrix */
      exp1  = exp( lambda * dtau_d );
      term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

      /* Eq. 6.82 */
      A = 1.0 / ( term1 - 1. / term1 );

      /* Eq. 6.88 and own calculations */
      input->a11[i] = A * 2.0 * lambda / alpha_2;
      input->a12[i] = A * ( exp1 - 1. / exp1 );
    }

    /* den1 needed for alpha_5 and alpha_6 */
    den1 = 1. / ( mu_0_inv * mu_0_inv - lambda * lambda );

    /* Eqs. 6.64 and 6.77 */
    alpha_3 = - omega_0_d * b_mmu_0;
    alpha_4 = omega_0_d + alpha_3;
    alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
    alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

    input->a33[i]      = exp ( - dtau_d  * mu_0_inv );   
    input->a33_unsc[i] = exp ( - dtau[i] * mu_0_inv );   

    input->a13[i] =
      + alpha_5 * ( 1.0 - input->a33[i] * input->a11[i] )
      - alpha_6 * input->a12[i];

    input->a23[i] =
      - alpha_5 * input->a33[i] * input->a12[i]
      + alpha_6 * ( input->a33[i] - input->a11[i] );

  } /* end do ilyr */

  return 0;
}


/*************************************************************/
/* multi_layer_solver                                        */
/*************************************************************/

int multi_layer_solver (rodents_input_struct   input,
			rodents_output_struct  output)
{
  double **matr, **mlu, *d, *b;
  long *indx;
  int i=0, n=2*input.nlyr+2;

  /* calloc memory */
  matr = calloc(n,sizeof(double *));
  for (i=0;i<n;i++)
    matr[i] = calloc(5,sizeof(double));

  mlu = calloc(n,sizeof(double *));
  for (i=0;i<n;i++)
    mlu[i] = calloc(2,sizeof(double));

  indx = calloc(n,sizeof(long));
  d    = calloc(n,sizeof(double));
  b    = calloc(n,sizeof(double));

  /******** RBUW start of important stuff ********/

  /* solve direct radiation */
  if (input.S_0 != 0.0) {
    solve_direct(input.a33, input.a33_unsc, input.nlyr, input.S_0,
		 output.s_direct, output.s_direct_unsc);
  }

  /* setup matrix and vector for diffuse radiation */
  setup_equations (input, output.s_direct, matr, b);

  /* solve diffuse radiation */
  bandec(matr, n, mlu, indx, d);
  banbks(matr, n, mlu, indx, b);
  
  /* put result */
  cp_to_result(b, input.nlyr, output.e_minus, output.e_plus);

  /******** RBUW end of important stuff ********/

  for (i=0;i<n;i++) {
    free(mlu[i]);
    free(matr[i]);
  }
  free(b);
  free(d);
  free(indx);
  free(mlu);
  free(matr);

  return 0;
}


/*************************************************************/
/* solve_direct                                              */
/*************************************************************/

static void solve_direct(double *a33, double *a33_unsc, int nlyr, double S_0,
			 double *S, double *S_unsc)
{
  int i=0;
  
  S[0]      = S_0;
  S_unsc[0] = S_0;
  
  for (i=0;i<nlyr;i++) {
    S[i+1]      = a33[i]      * S[i];
    S_unsc[i+1] = a33_unsc[i] * S_unsc[i];
  }
}


/*************************************************************/
/* setup_equations                                           */
/*************************************************************/

static void setup_equations ( /* input */
			     rodents_input_struct  input,
			     double               *s_direct,
			     /* output */
			     double              **matr,
			     double               *rhs )
{
  int i=0, i2=0;

  /* the matrix: */

  matr[0][2] = -1.0;

  for (i=0;i<input.nlyr;i++) {
    i2=2*i+1;
    matr[i2  ][2] = -1.0;
    matr[i2  ][4] = input.a11[i];
    matr[i2  ][1] = input.a12[i];
    matr[i2+1][2] = -1.0;
    matr[i2+1][0] = input.a11[i];
    matr[i2+1][3] = input.a12[i];
  }

  i=input.nlyr;
  i2=2*i+1;
  matr[i2  ][2] = -1.0;
  matr[i2  ][1] = input.albedo;

  /* the RHS */

  for (i=0;i<input.nlyr;i++) {
    i2=2*i+1;
    rhs[i2  ] += - input.a13[i] * s_direct[i];
    rhs[i2+1] += - input.a23[i] * s_direct[i];
  }

  i=input.nlyr;
  i2=2*i;
  rhs[i2+1] += - input.albedo * s_direct[i] * input.mu_0;

}


/*************************************************************/
/* cp_to_result                                              */
/*************************************************************/

static void cp_to_result ( double *b,
			   int     nlyr,
			   double *e_minus,
			   double *e_plus )
{
  int i=0;
  for (i=0;i<nlyr+1;i++) {
    e_minus[i] = b[2*i];
    e_plus [i] = b[2*i+1];
  }
}

/*************************************************************/
/* process_output                                            */
/*************************************************************/

static int process_output ( /* input */
			   rodents_input_struct   input,
			   rodents_output_struct  output,
			   /* output */
			   double                *fldn,
			   double                *flup,
			   double                *fldir )
{
  int i=0;

  for (i=0;i<input.nlyr+1;i++)  {
    fldn [i] = output.e_minus [i]
      + ( output.s_direct[i] - output.s_direct_unsc[i] ) * input.mu_0;
    flup [i] = output.e_plus  [i];
    fldir[i] = output.s_direct_unsc[i] * input.mu_0;
  }

  return 0;
}

/*************************************************************/
/* free_input                                                */
/*************************************************************/

static void free_input(rodents_input_struct input)
{
  free(input.a11);
  free(input.a12);
  free(input.a13);
  free(input.a23);
  free(input.a33);
  free(input.a33_unsc);
}

/*************************************************************/
/* err_out                                                   */
/*************************************************************/

static inline int err_out ( char *out, int status )
{
  fprintf(stderr,out,status);
  return status;
}
