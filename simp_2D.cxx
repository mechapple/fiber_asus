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

#include <iostream>
#include <cmath>
#include <ctime>

#include <lj.h>

double func(double x1,double x2) { 
	double y = sqrt(1 + x1*x1 + x2*x2);
	return 1.0/y; 
}

int main(int argc, char **argv)
{
	int start_s=clock();
	int Nsim = 1000;
	//double L = atof(argv[1]);
	
	for(double L = 0.09; L<0.4; L+=0.001) {
		printf("x,LJ_energy %.6f %.15f\n",L,pair_energy(L));
	}
	printf("\n");
	
	for(double L = 0.09; L<0.4; L+=0.001) {
		double h = 1.0/(Nsim), sum = 0.0,x1,x2;
		int m[2][2] = {{4,8},{8,16}};
		
		for(int i=0;i<=Nsim;i++) 
			for(int j=0;j<=Nsim;j++) {
				//x1 = i*h; x2 = j*h;
				//sum += m[i%2][j%2]*func(i*h,j*h);
				sum += m[i%2][j%2]*pair_energy(sqrt(pow((i*h-j*h),2)+L*L));
		}
		
		int n[2] = {2,4};
		for(int i=0;i<=Nsim;i++) 
			for(double t=0;t<=1.0;t+=1.0) {
				//x = i*h; y = j*h;
				sum -= n[i%2]*pair_energy(sqrt(pow((i*h-t),2)+L*L));
				sum -= n[i%2]*pair_energy(sqrt(pow((t-i*h),2)+L*L));
				//sum -= n[i%2]*func(i*h,t);
				//sum -= n[i%2]*func(t,i*h);
		}
		
		sum += pair_energy(sqrt(pow((0.0),2)+L*L))+pair_energy(sqrt(pow((-1.0),2)+L*L))+pair_energy(sqrt(pow((1.0),2)+L*L))+pair_energy(sqrt(pow((0.0),2)+L*L));
		
		sum *= h*h/9.0;
		//double error = fabs(sum - 0.793359121);
		//std::cout << "Integral = " << sum*h*h << "\n";
		printf("L,h,Integral %.6f %.6f %.15f\n",L,h,sum);
	}
	
	int stop_s=clock();
	std::cout << "time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << std::endl;
	return 0;
}

