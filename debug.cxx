/*
 * main.cxx 
 * MPASS - Massively Parallel Adaptive String Simulator
 * 
 * //(0,1),(c,1),(1,c),(1,0) c=0.55191502449 for unit circle
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

double yfac=10.0,xfac=0.1;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <cubature.h>
#include <bernstein.h>
#include <eigen2/Eigen/Dense>

#include <lj.h>
#include <point.h>
#include <bezier.h>

#include <min.h>

int main(int argc, char **argv)
{
	//Bezier b,b_min;
	std::vector<Bezier> b(Nbez),b_min(Nbez);
	double sep = sig*pow(2.0,1.0/6);
	if(argc==2) sep = atof(argv[1]);
	//FILE *fcom = fopen("Bezier.lammpstrj","w");
	//FILE *fcom2 = fopen("Bezier_cps.lammpstrj","w");
	
	std::ifstream inputfile;
	inputfile.open("InitialCPs.dat");
	
	double dump[7];
	for(int i=0;i<Nbez;i++) 
		for(int j=0;j<4;j++) {
			inputfile >> dump[0] >> dump[1] >> dump[2] >> dump[3];
			inputfile >> b[i].p[j].x >> b[i].p[j].y >> b[i].p[j].z;
	}
	
	
	//b[0].p[0].x = 0.0; b[0].p[0].y = -0.5*sep; b[0].p[0].z = 0.0;
	//b[0].p[1].x = 0.2; b[0].p[1].y = -0.5*sep; b[0].p[1].z = 0.0;
	//b[0].p[2].x = 0.8; b[0].p[2].y = -2.0*sep; b[0].p[2].z = 0.0;
	//b[0].p[3].x = 1.0; b[0].p[3].y = -5.0*sep; b[0].p[3].z = 0.0;
	
	//b[1].p[0].x = 0.0; b[1].p[0].y = 0.5*sep; b[1].p[0].z = 0.0;
	//b[1].p[1].x = 0.2; b[1].p[1].y = 0.5*sep; b[1].p[1].z = 0.0;
	//b[1].p[2].x = 0.8; b[1].p[2].y = 2.0*sep; b[1].p[2].z = 0.0;
	//b[1].p[3].x = 1.0; b[1].p[3].y = 5.0*sep; b[1].p[3].z = 0.0;
	
	//double yfac=10.0,xfac=0.1;
	//b[0].p[0].x = 0.0*xfac; b[0].p[0].y = 0.0*yfac; b[0].p[0].z = 0.0;
	//b[0].p[1].x = 0.2*xfac; b[0].p[1].y = 0.2*yfac; b[0].p[1].z = 0.0;
	//b[0].p[2].x = 0.8*xfac; b[0].p[2].y = 0.8*yfac; b[0].p[2].z = 0.0;
	//b[0].p[3].x = 1.0*xfac-0.5*sep; b[0].p[3].y = 1.0*yfac; b[0].p[3].z = 0.0;
	
	//b[1].p[0].x = 2.0*xfac; b[1].p[0].y = 0.0*yfac; b[1].p[0].z = 0.0;
	//b[1].p[1].x = 1.8*xfac; b[1].p[1].y = 0.2*yfac; b[1].p[1].z = 0.0;
	//b[1].p[2].x = 1.2*xfac; b[1].p[2].y = 0.8*yfac; b[1].p[2].z = 0.0;
	//b[1].p[3].x = 1.0*xfac+0.5*sep; b[1].p[3].y = 1.0*yfac; b[1].p[3].z = 0.0;
	
	
	printf("Initial constraints %.12f %.12f %.12f\n",b[0].cons[0].x,b[0].cons[0].y,b[0].cons[0].z);
	printf("Initial constraints %.12f %.12f %.12f\n",b[1].cons[0].x,b[1].cons[0].y,b[1].cons[0].z);
	
	b[0].cons[0].x = 0.0; b[0].cons[0].y = 0.0; b[0].cons[0].z = 0.0; //fix left end
	//for(int j=0;j<4;j++) {b[0].cons[j].x = 0.0; b[0].cons[j].y = 0.0; b[0].cons[j].z = 0.0; } // fix bottom fiber
	
	b[1].cons[0].x = 0.0; b[1].cons[0].y = 0.0; b[1].cons[0].z = 0.0; //fix left end
	
	printf("Modified constraints %.12f %.12f %.12f\n",b[0].cons[0].x,b[0].cons[0].y,b[0].cons[0].z);
	printf("Modified constraints %.12f %.12f %.12f\n",b[1].cons[0].x,b[1].cons[0].y,b[1].cons[0].z);
	
	//b.p[0].x = 0.0; b.p[0].y = 1.0; b.p[0].z = 0.0;
	//b.p[1].x = 0.55191502449; b.p[1].y = 1.0; b.p[1].z = 0.0;
	//b.p[2].x = 1.0; b.p[2].y = 0.55191502449; b.p[2].z = 0.0;
	//b.p[3].x = 1.0; b.p[3].y = 0.0; b.p[3].z = 0.0;
	
	//b[0].p[0].x = 0.0; b[0].p[0].y = 1.0; b[0].p[0].z = 0.0;
	//b[0].p[1].x = 0.55191502449; b[0].p[1].y = 1.0; b[0].p[1].z = 0.0;
	//b[0].p[2].x = 1.0; b[0].p[2].y = 0.55191502449; b[0].p[2].z = 0.0;
	//b[0].p[3].x = 1.0; b[0].p[3].y = 0.0; b[0].p[3].z = 0.0;
	
	int fac = 1;
	//b[1].p[0].x = 0.0; b[1].p[0].y = 1.0*fac; b[1].p[0].z = sep;
	//b[1].p[1].x = 0.55191502449*fac; b[1].p[1].y = 1.0*fac; b[1].p[1].z = sep;
	//b[1].p[2].x = 1.0*fac; b[1].p[2].y = 0.55191502449*fac; b[1].p[2].z = sep;
	//b[1].p[3].x = 1.0*fac; b[1].p[3].y = 0.0; b[1].p[3].z = sep;
	
	
	
	printf("Bezier length %.15f %.15f\n",b[0].length(),b[1].length());
	printf("Bezier axial energy %.15f\n",b[0].axial_energy()+b[1].axial_energy());
	printf("Bezier bending energy %.15f\n",b[0].bending_energy()+b[1].bending_energy());
	printf("Bezier inter energy %.15f\n",inter_energy(b[0],b[1]));
	
	b[0].intra_force();
	b[0].compute_mass();
	
	vecBez bz(2); bz = inter_force(b[0],b[1]);
	
	//b_min[0] = fire_min0(b[0],fcom);
	//b_min = fire_min(b,fcom,fcom2);
	
	//b.axial_force();
	//b.bending_force();
	//std::cout << "Bezier length " << b.length() << "\n";
	
	//fclose(fcom); fclose(fcom2);
	
	return 0;
}

