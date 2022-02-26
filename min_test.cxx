/*
 * cg_test.cxx
 * test code for conjugate minimizaton with backtracking line search
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
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>

#include <random>
#include <algorithm>
#include <iterator>
#include <functional>

#define PI 3.14159265359
#define TOL 1e-8
#define MAX 1e8

typedef std::vector<double> vec;

struct ext_vector {
	vec x,v,f;
};

class quartic {
	public:
		int N; double mass;
		vec a,b;
		double value(vec);
		vec gradient(vec);
		vec normalize(vec);
		vec minimize(vec);
		vec fire_min(vec);
		ext_vector integrate_verlet(vec,vec,vec,double);
		vec scalar_multiply(vec,double);
		vec add(vec,vec);
		double dot(vec,vec);
		double mag(vec);
		void initialize_vectors();
		void display_function();
		void print_vector(vec);
};

void quartic::initialize_vectors() {
	for(int i=1;i<=N;i++)
	{
		a.push_back(fabs(sin(i*PI*0.5/N)));
		b.push_back(fabs(cos(i*PI*0.5/N)));
	}
}

void quartic::print_vector(vec x) {
	for(int i=0;i<N;i++)
	{
		if(i!=0) printf(" ");
		printf(" %.6f",x[i]);
	}
	printf("\n");
}


void quartic::display_function() {
	for(int i=0;i<N;i++)
	{
		if(i!=0) printf("+");
		printf(" %.2f*(x%d-%.2f)^4 ",a[i],i+1,b[i]);
	}
	printf("\n");
}

double quartic::value(vec x) {
	double f=0;
	for(int i=0;i<N;i++) {
		f += a[i]*pow((x[i]-b[i]),4);
	}
	return f;
}

vec quartic::gradient(vec x) {
	vec df;
	for(int i=0;i<N;i++) {
		df.push_back(4.0*a[i]*pow((x[i]-b[i]),3));
	}
	return(df);
}

vec quartic::normalize(vec x) {
	vec xn;
	double mag=0;
	for(int i=0;i<N;i++) {
		mag += x[i]*x[i];
	}
	xn = scalar_multiply(x,1.0/sqrt(mag));
	return(xn);
}


vec quartic::scalar_multiply(vec x, double a) {
	vec mx;
	for(int i=0;i<N;i++) {
		mx.push_back(x[i]*a);
	}
	return(mx);
}


vec quartic::add(vec x1, vec x2) {
	vec y=x1;
	for(int i=0;i<N;i++) {
		y[i] = x1[i] + x2[i];
	}
	return y;
}

double quartic::dot(vec x1, vec x2) {
	double sum=0;
	for(int i=0;i<N;i++) {
		sum += x1[i]*x2[i];
	}
	return sum;
}

double quartic::mag(vec x1) {
	double sum=0;
	for(int i=0;i<N;i++) {
		sum += x1[i]*x1[i];
	}
	return sqrt(sum);
}

vec quartic::minimize(vec x0) {
	
	vec xk,pk,gk,xl,xk1,pk1,gk1;
	xk1 = x0;
	double c1=0.0001,c2=0.1,tau=0.5,alpha_k,alpha,beta_k;
	
	alpha_k=1.0/tau;
	gk1 = normalize(gradient(xk1)); pk1 = scalar_multiply(gk1,-1);
	
	do {
		alpha_k=1.0/tau;
		
		pk = pk1; xk = xk1; gk = gk1;
		
		do { //backtracking Wolfe conditions
			alpha_k *= tau;
			xl = add(xk,scalar_multiply(pk,alpha_k));
			printf("\n\tBacktracking %lf %lf %lf",alpha_k,value(xl),dot(pk,pk));
		}while( value(xl)>(value(xk)+alpha_k*c1*dot(gk,pk)) || dot(gradient(xl),pk)<c2*dot(gk,pk) );
		
		xk1 = add(xk,scalar_multiply(pk,alpha_k));
		gk1 = gradient(xk1);
		//beta_k = dot(gk1,gk1)/dot(gk,gk);
		beta_k = dot(gk1,add(gk1,pk))/dot(gk,gk);
		//beta_k = dot(gk1,add(gk1,pk))/dot(pk,add(gk1,pk));
		
		pk1 = add(scalar_multiply(gk1,-1),scalar_multiply(pk,beta_k));
		pk1 = normalize(pk1);
		printf("\nIteration %lf %lf %lf %lf",value(xk1),beta_k,dot(gk1,gk1),dot(pk1,pk1));
	}while(fabs(dot(gk1,gk1))>TOL);
	
	return xk1;
}

ext_vector quartic::integrate_verlet(vec xk,vec vk,vec fk,double dt) {
	ext_vector Xk1;
	vec dxk = add(scalar_multiply(vk,dt),scalar_multiply(fk,0.5*dt*dt/mass));
	Xk1.x = add(xk,dxk);
	Xk1.f = scalar_multiply(gradient(Xk1.x),-1);
	vec dvk = scalar_multiply(add(fk,Xk1.f),0.5*dt/mass);
	Xk1.v = add(vk,dvk);
	return Xk1;
}

vec quartic::fire_min(vec x0) {
	mass = 1.0;
	vec xk(N,0.0),vk(N,0.0),fk(N,0.0),xk1(N,0.0),vk1(N,0.0),fk1(N,0.0);
	ext_vector Xk,Xk1;
	
	double dt = 0.01; int count=0;
	double alpha0 = 0.1,alpha = alpha0, fdec = 0.5, finc = 1.1, dtmax = 1.0, falpha = 0.99;
	
	xk = x0; fk = scalar_multiply(gradient(xk),-1);
	std::fill(vk.begin(), vk.end(), 0.0);
	
	do {
		double P = dot(fk,vk);
		if(P<0) {
			std::fill(vk.begin(), vk.end(), 0.0);
			dt *= fdec; alpha = alpha0;
		}else {
			double cf = alpha*sqrt(dot(vk,vk)/dot(fk,fk));
			double cv = 1.0 - alpha;
			vk = add(scalar_multiply(vk,cv),scalar_multiply(fk,cf));
			dt = std::min(dt*finc,dtmax); alpha *= falpha;
		}
		
		//printf("\n	Iteration prior %d %lf %lf",count,value(xk),mag(xk));
		Xk1 = integrate_verlet(xk,vk,fk,dt); //Integration via velocity verlet
		xk = Xk1.x; vk = Xk1.v; fk = Xk1.f;
		if(count%100==0) printf("\nIteration %d %.12f %.12f",count,value(xk),mag(fk));
		count++ ;
	}while(mag(fk)>TOL && count<MAX);
	
	return xk;
}

int main(int argc, char **argv)
{
	std::clock_t start;
    double duration;
    start = std::clock();
    std::cout<<"Starting\n";
    
	quartic q;
	if(argc>1) q.N=atoi(argv[1]);
	else q.N=10;
	
	q.initialize_vectors();
	//q.display_function();
	
	std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine(rnd_device());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    auto gen = std::bind(dist, mersenne_engine);
    vec x(q.N),x_min;
    std::generate(begin(x), end(x), gen);
	
	//printf(" Initial Vector\n"); q.print_vector(x);
	//printf(" Initial Reverse Vector\n"); q.print_vector(q.scalar_multiply(x,-1));
	//printf(" %lf\n",q.value(x));
	x_min = q.fire_min(x);
	
	printf("\n Initial Vector "); q.print_vector(x);
	printf("\n Function form "); q.display_function();
	printf("\n Final Vector "); q.print_vector(x_min);
	//printf(" %lf\n",q.value(x_min));
	//q.print_vector(q.gradient(x_min));
	
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Duration "<< duration <<"seconds \n";
	
	return 0;
}

