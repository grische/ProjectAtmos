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

#ifndef __rodents_h
#define __rodents_h

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct {
  int    thermal;
  int    solar;
  int    nlyr;
  double S_0;
  double mu_0;
  double albedo;
  double *a11, *a12, *a13, *a23, *a33, *a33_unsc;
} rodents_input_struct;

typedef struct {
  double *e_minus, *e_plus, *s_direct, *s_direct_unsc;
} rodents_output_struct;


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
		   double *fldir );


#endif
