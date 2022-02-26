/*
 * multimin.cxx
 * 
 * Copyright 2018 Anirban Pal <anirban@ZeroPointEnergy>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

/* Paraboloid centered on (p[0],p[1]), with  
   scale factors (p[2],p[3]) and minimum p[4] */

double my_f (const gsl_vector *v, void *params)
{
	double x, y;
	double *p = (double *)params;
	
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);
	
	return p[2] * (x - p[0]) * (x - p[0]) +
		   p[3] * (y - p[1]) * (y - p[1]) + p[4]; 
}

/* The gradient of f, df = (df/dx, df/dy). */
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
	double x, y;
	double *p = (double *)params;
	
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);
	
	gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
	gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
	*f = my_f(x, params); 
	my_df(x, params, df);
}

int main(int argc, char **argv)
{
	double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};

	const gsl_multimin_fminimizer_type *T = 
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;
	
	/* Starting point */
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, 5.0);
	gsl_vector_set (x, 1, 7.0);
	
	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 1.0);
	
	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = my_f;
	minex_func.params = par;
	
	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	
	do
	{
	  iter++;
	  status = gsl_multimin_fminimizer_iterate(s);
	  
	  if (status) 
		break;
	
	  size = gsl_multimin_fminimizer_size (s);
	  status = gsl_multimin_test_size (size, 1e-2);
	
	  if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}
	
	  printf ("%5zu %10.3e %10.3e f() = %7.3f size = %.3f\n", 
			  iter,
			  gsl_vector_get (s->x, 0), 
			  gsl_vector_get (s->x, 1), 
			  s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);
	
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	
	return status;
}

