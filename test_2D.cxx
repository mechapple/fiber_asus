/*
 * num_int_2D_test.cxx
 * 
 * Copyright 2017 Anirban <anirban@ZeroPointEnergy>
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


#include "cuba.h"

#include <iostream>
#include <cmath>
#include <ctime>
#include "cubature.h"
#define PI 3.14159265358979323846264338327950288419716939937510582

double theta = 90.0*PI/180.0;
double costheta = cos(theta);
double sintheta = sin(theta);

double ar = 1.0e3;
double sig = 1.0e-3;
//double x20 = 0.5*(1-costheta);
//double y20 = 0.5*(-sintheta);

double x20 = 0.0;
double y20 = 0.0;

double func(double x1,double x2) { 
	//double y = sqrt(1 + x1*x1 + x2*x2);
	double dx = x2*costheta+x20-x1;
	double dy = x2*sintheta+y20;
	double rn = exp(-ar*( sqrt(dx*dx+dy*dy)-sig ));
	return 1.0*(rn*rn-2*rn);
	//return 1.0/y; 
}

int f1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    fval[0] = func(x[0],x[1]);
    return 0; // success
}

static int f2(const int *ndim, const cubareal x[],
  const int *ncomp, cubareal f[], void *userdata) {
	
	f[0] = func(x[0],x[1]);

  return 0;
}

int main(int argc, char **argv)
{
	int start_s=clock();
	int Nsim = 1000;
	double h = 1/(2.0*Nsim), sum = 0.0, lsum[3], mid1, mid2;
	double lam1 = 3.0*(35.0+sqrt(385.0))/140.0, lam2 = 3.0*(35.0-sqrt(385.0))/140.0;
	double b1 = (77.0-3.0*sqrt(385.0))/891.0, b2 = (77.0+3.0*sqrt(385.0))/891.0;
	double mu = sqrt(3.0/5.0), c = 25.0/324.0;
	double m[2][2] = {{4,8},{8,16}};
		
	for(int i=0;i<Nsim;i++) 
		for(int j=0;j<Nsim;j++) {
			mid1 = i*2*h + h; mid2 = j*2*h + h;
			
			lsum[0] = func(mid1 + lam1*h,mid2) + func(mid1 - lam1*h,mid2) + func(mid1,mid2 + lam1*h) + func(mid1,mid2 - lam1*h);
			lsum[1] = func(mid1 + lam2*h,mid2) + func(mid1 - lam2*h,mid2) + func(mid1,mid2 + lam2*h) + func(mid1,mid2 - lam2*h);
			lsum[2] = func(mid1 - mu*h,mid2 - mu*h) + func(mid1 - mu*h,mid2 + mu*h) + func(mid1 + mu*h,mid2 - mu*h) + func(mid1 + mu*h,mid2 + mu*h);
			sum += b1*lsum[0] + b2*lsum[1] + c*lsum[2];
	}
	
	sum *= 4.0*h*h;
	//double error = fabs(sum - 0);
	//std::cout << "Integral = " << sum*h*h << "\n";
	printf("Computed Integral with Simpson's rule %.6f %.15f %.15g\n",h,sum,sum);
	
	int stop_s=clock();
	std::cout << "time (s): " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	
	double xmin[2] = {0,0}, xmax[2] = {1,1}, sigma = 0.5, val[1], err[1];
    hcubature(1, f1, &sigma, 2, xmin, xmax, 0, 0, 1e-8, ERROR_INDIVIDUAL, val, err);
    printf("Computed integral with hcubature = %0.15f +/- %.15g\n", val[0], err[0]);
	
	int stop_s1=clock();
	std::cout << "time (s): " << (stop_s1-stop_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	
	pcubature(1, f1, &sigma, 2, xmin, xmax, 0, 0, 1e-8, ERROR_INDIVIDUAL, val, err);
    printf("Computed integral with pcubature = %0.15f +/- %.15g\n", val[0], err[0]);
	
	int stop_s2=clock();
	std::cout << "time (s): " << (stop_s2-stop_s1)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	
	int comp, nregions, neval, fail;
	cubareal integral[1], error[1], prob[1];
	
	Cuhre(2, 1, f2, NULL, 1, 1e-8, 1e-12, 0 | 4, 0, 50000, 0, NULL, NULL,
	      &nregions, &neval, &fail, integral, error, prob);
	printf("Computed integral with cuba-cuhre = %0.15f +/- %.15g\n", (double)integral[0], (double)error[0]);
    
    int stop_s3=clock();
	std::cout << "time (s): " << (stop_s3-stop_s2)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	
	return 0;
}

