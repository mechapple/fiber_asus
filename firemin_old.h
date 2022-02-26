#define Nsimp2D 1e3
#define xscale 50.0

typedef std::vector<Bezier> vecBez;

void display_allcurves(vecBez bk,int time,FILE *fcom) {
	double min=-100.0,max=100.0;
	fprintf(fcom,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp",time,(Nres+1)*Nbez);
	fprintf(fcom,"\n%.3f %.3f \n%.3f %.3f \n%.3f %.3f \nITEM: ATOMS id mol type q xu yu zu\n",min,max,min,max,min,max);
	double h = 1.0/Nres;
	int count=0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0; j<=Nres; j++) {
			double t = j*h; point R = bk[i].r(t);
			fprintf(fcom,"%d 1 1 0.00 %lf %lf %lf\n",++count,R.x*xscale,R.y,R.z);
		}
	}
}

void display_all_cps(vecBez bk,int time,FILE *fcom) {
	double min=-100.0,max=100.0;
	fprintf(fcom,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp",time,4*Nbez);
	fprintf(fcom,"\n%.3f %.3f \n%.3f %.3f \n%.3f %.3f \nITEM: ATOMS id mol type q xu yu zu\n",min,max,min,max,min,max);
	int count=0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0; j<4; j++) {
			fprintf(fcom,"%d 1 1 0.00 %lf %lf %lf\n",++count,bk[i].p[j].x*xscale,bk[i].p[j].y,bk[i].p[j].z);
		}
	}
}

double inter_energy(Bezier b1,Bezier b2) {
	int Nsim = Nsimp2D;
	
	int m[2][2] = {{4,8},{8,16}};
	
	double tmin = 0.9,tmax=1.0;
	double h = (tmax-tmin)/(Nsim), sum = 0.0,x,y;
	
	for(int i=0;i<=Nsim;i++) 
		for(int j=0;j<=Nsim;j++) {
			double t1 = tmin+i*h,t2 = tmin+j*h; point R1 = b1.r(t1), R2 = b2.r(t2);
			double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
			sum += m[i%2][j%2]*pair_energy(d)*normpoint(b1.dr(t1))*normpoint(b2.dr(t2));
			//printf("\n value %.15f %.15f %.15f %.15f %.15f",t1,t2,d,pair_energy(d),pair_energy(d)*normpoint(b1.dr(t1))*normpoint(b2.dr(t2)));
	}
	
	int n[2] = {2,4};
	for(int i=0;i<=Nsim;i++) 
		for(double t=tmin;t<=tmax;t+=(tmax-tmin)) {
			double t1 = tmin+i*h; 
			
			point R1 = b1.r(t1), R2 = b2.r(t);
			double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
			sum -= n[i%2]*pair_energy(d)*normpoint(b1.dr(t1))*normpoint(b2.dr(t));
			
			R1 = b1.r(t); R2 = b2.r(t1);
			d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
			sum -= n[i%2]*pair_energy(d)*normpoint(b1.dr(t))*normpoint(b2.dr(t1));
			
			//sum -= n[i%2]*func(i*h,t);
			//sum -= n[i%2]*func(t,i*h);
	}
	
	for(double t1=tmin;t1<=tmax;t1+=(tmax-tmin)) for(double t2=tmin;t2<=tmax;t2+=(tmax-tmin)) {
		point R1 = b1.r(t1), R2 = b2.r(t2);
		double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
		sum += pair_energy(d)*normpoint(b1.dr(t1))*normpoint(b2.dr(t2));
	}
	sum *= h*h/9.0;
	
	//printpoint(b2.p[3]);
	printf("\nOrig %.15f \n", sum);
		
	double xmin[2] = {tmin,tmin}, xmax[2] = {tmax,tmax}, val[1], err[1];
	point rp[8] = {b1.p[0],b1.p[1],b1.p[2],b1.p[3],b2.p[0],b2.p[1],b2.p[2],b2.p[3]};
	
	//printpoint(rp[7]);
	pcubature(1, f_InterEnergy, rp, 2, xmin, xmax, 0, 0, TOL2, ERROR_INDIVIDUAL, val, err);
	printf("Computed integral = %0.15f +/- %.15g\n", val[0], err[0]);
	
	
	return sum;
}

vecBez inter_force(Bezier b1,Bezier b2) {
	vecBez bz(2);
	int Nsim = Nsimp2D;
	
	double tmin = 0.7,tmax=1.0;
	double h = (tmax-tmin)/(Nsim), sum = 0.0,x,y;
	int m[2][2] = {{4,8},{8,16}};
	point f[2][4];
	for(int i=0;i<2;i++) for(int j=0;j<4;j++) { f[i][j].x = 0.0; f[i][j].y = 0.0; f[i][j].z = 0.0; }
	
	for(int i=0;i<=Nsim;i++) 
		for(int j=0;j<=Nsim;j++) {
			double t1 = tmin+i*h,t2 = tmin+j*h; 
			point R1 = b1.r(t1), R2 = b2.r(t2);
			point dR1 = b1.dr(t1), dR2 = b2.dr(t2);
			double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
			
			if(d<=rcut) {
				
				double ds[2] = { normpoint(dR1),normpoint(dR2) };
				double pg = pair_gradient(d), pe = pair_energy(d);
				double Bt[4][2] = { {B[0](t1),B[0](t2)}, {B[1](t1),B[1](t2)}, {B[2](t1),B[2](t2)}, {B[3](t1),B[3](t2)} };
				double dBt[4][2] = { {dB[0](t1),dB[0](t2)}, {dB[1](t1),dB[1](t2)}, {dB[2](t1),dB[2](t2)}, {dB[3](t1),dB[3](t2)} };
				
				for(int k=0;k<4;k++) {
					
					f[0][k].x += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R1.x-R2.x)*Bt[k][0]/d  +  pe*dR1.x*dBt[k][0]*ds[1]/ds[0]);
					f[0][k].y += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R1.y-R2.y)*Bt[k][0]/d  +  pe*dR1.y*dBt[k][0]*ds[1]/ds[0]);
					f[0][k].z += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R1.z-R2.z)*Bt[k][0]/d  +  pe*dR1.z*dBt[k][0]*ds[1]/ds[0]);
					
					f[1][k].x += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R2.x-R1.x)*Bt[k][1]/d  +  pe*dR2.x*dBt[k][1]*ds[0]/ds[1]);
					f[1][k].y += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R2.y-R1.y)*Bt[k][1]/d  +  pe*dR2.y*dBt[k][1]*ds[0]/ds[1]);
					f[1][k].z += m[i%2][j%2]*(ds[0]*ds[1]*pg*(R2.z-R1.z)*Bt[k][1]/d  +  pe*dR2.z*dBt[k][1]*ds[0]/ds[1]);
				}
			}
			
			//printf("(%lf,%lf) (%lf,%lf) %lf %lf %lf %lf\n",R1.x,R1.y,R2.x,R2.y,t1,t2,d,pair_gradient(d));
	}
	
	int n[2] = {2,4};
	for(int i=0;i<=Nsim;i++) 
		for(double t=tmin;t<=tmax;t+=(tmax-tmin)) {
			 
			{	double t1 = tmin+i*h, t2 = t;
				
				point R1 = b1.r(t1), R2 = b2.r(t2);
				point dR1 = b1.dr(t1), dR2 = b2.dr(t2);
				double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
				
				if(d<=rcut) {
					
					double ds[2] = { normpoint(dR1),normpoint(dR2) };
					double pg = pair_gradient(d), pe = pair_energy(d);
					double Bt[4][2] = { {B[0](t1),B[0](t2)}, {B[1](t1),B[1](t2)}, {B[2](t1),B[2](t2)}, {B[3](t1),B[3](t2)} };
					double dBt[4][2] = { {dB[0](t1),dB[0](t2)}, {dB[1](t1),dB[1](t2)}, {dB[2](t1),dB[2](t2)}, {dB[3](t1),dB[3](t2)} };
				
					for(int k=0;k<4;k++) {
						
						f[0][k].x -= n[i%2]*(ds[0]*ds[1]*pg*(R1.x-R2.x)*Bt[k][0]/d  +  pe*dR1.x*dBt[k][0]*ds[1]/ds[0]);
						f[0][k].y -= n[i%2]*(ds[0]*ds[1]*pg*(R1.y-R2.y)*Bt[k][0]/d  +  pe*dR1.y*dBt[k][0]*ds[1]/ds[0]);
						f[0][k].z -= n[i%2]*(ds[0]*ds[1]*pg*(R1.z-R2.z)*Bt[k][0]/d  +  pe*dR1.z*dBt[k][0]*ds[1]/ds[0]);
						
						f[1][k].x -= n[i%2]*(ds[0]*ds[1]*pg*(R2.x-R1.x)*Bt[k][1]/d  +  pe*dR2.x*dBt[k][1]*ds[0]/ds[1]);
						f[1][k].y -= n[i%2]*(ds[0]*ds[1]*pg*(R2.y-R1.y)*Bt[k][1]/d  +  pe*dR2.y*dBt[k][1]*ds[0]/ds[1]);
						f[1][k].z -= n[i%2]*(ds[0]*ds[1]*pg*(R2.z-R1.z)*Bt[k][1]/d  +  pe*dR2.z*dBt[k][1]*ds[0]/ds[1]);
					}
				}
			}
			
			{	double t1 = t, t2 = tmin+i*h;
				
				point R1 = b1.r(t1), R2 = b2.r(t2);
				point dR1 = b1.dr(t1), dR2 = b2.dr(t2);
				double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
				
				if(d<=rcut) {
					
					double ds[2] = { normpoint(dR1),normpoint(dR2) };
					double pg = pair_gradient(d), pe = pair_energy(d);
					double Bt[4][2] = { {B[0](t1),B[0](t2)}, {B[1](t1),B[1](t2)}, {B[2](t1),B[2](t2)}, {B[3](t1),B[3](t2)} };
					double dBt[4][2] = { {dB[0](t1),dB[0](t2)}, {dB[1](t1),dB[1](t2)}, {dB[2](t1),dB[2](t2)}, {dB[3](t1),dB[3](t2)} };
					
					for(int k=0;k<4;k++) {
						
						f[0][k].x -= n[i%2]*(ds[0]*ds[1]*pg*(R1.x-R2.x)*Bt[k][0]/d  +  pe*dR1.x*dBt[k][0]*ds[1]/ds[0]);
						f[0][k].y -= n[i%2]*(ds[0]*ds[1]*pg*(R1.y-R2.y)*Bt[k][0]/d  +  pe*dR1.y*dBt[k][0]*ds[1]/ds[0]);
						f[0][k].z -= n[i%2]*(ds[0]*ds[1]*pg*(R1.z-R2.z)*Bt[k][0]/d  +  pe*dR1.z*dBt[k][0]*ds[1]/ds[0]);
						
						f[1][k].x -= n[i%2]*(ds[0]*ds[1]*pg*(R2.x-R1.x)*Bt[k][1]/d  +  pe*dR2.x*dBt[k][1]*ds[0]/ds[1]);
						f[1][k].y -= n[i%2]*(ds[0]*ds[1]*pg*(R2.y-R1.y)*Bt[k][1]/d  +  pe*dR2.y*dBt[k][1]*ds[0]/ds[1]);
						f[1][k].z -= n[i%2]*(ds[0]*ds[1]*pg*(R2.z-R1.z)*Bt[k][1]/d  +  pe*dR2.z*dBt[k][1]*ds[0]/ds[1]);
					}
				}
			}
	}
	
	for(double t1=tmin;t1<=tmax;t1+=(tmax-tmin)) for(double t2=tmin;t2<=tmax;t2+=(tmax-tmin)) {
		point R1 = b1.r(t1), R2 = b2.r(t2);
		point dR1 = b1.dr(t1), dR2 = b2.dr(t2);
		double d = sqrt(pow(R1.x-R2.x,2)+pow(R1.y-R2.y,2)+pow(R1.z-R2.z,2));
		
		double ds[2] = { normpoint(dR1),normpoint(dR2) };
		double pg = pair_gradient(d), pe = pair_energy(d);
		double Bt[4][2] = { {B[0](t1),B[0](t2)}, {B[1](t1),B[1](t2)}, {B[2](t1),B[2](t2)}, {B[3](t1),B[3](t2)} };
		double dBt[4][2] = { {dB[0](t1),dB[0](t2)}, {dB[1](t1),dB[1](t2)}, {dB[2](t1),dB[2](t2)}, {dB[3](t1),dB[3](t2)} };
		
		if(d<=rcut) {
			
			double ds[2] = { normpoint(dR1),normpoint(dR2) };
			double pg = pair_gradient(d), pe = pair_energy(d);
			double Bt[4][2] = { {B[0](t1),B[0](t2)}, {B[1](t1),B[1](t2)}, {B[2](t1),B[2](t2)}, {B[3](t1),B[3](t2)} };
			double dBt[4][2] = { {dB[0](t1),dB[0](t2)}, {dB[1](t1),dB[1](t2)}, {dB[2](t1),dB[2](t2)}, {dB[3](t1),dB[3](t2)} };
				
			for(int k=0;k<4;k++) {
				
				f[0][k].x += 1.0*(ds[0]*ds[1]*pg*(R1.x-R2.x)*Bt[k][0]/d  +  pe*dR1.x*dBt[k][0]*ds[1]/ds[0]);
				f[0][k].y += 1.0*(ds[0]*ds[1]*pg*(R1.y-R2.y)*Bt[k][0]/d  +  pe*dR1.y*dBt[k][0]*ds[1]/ds[0]);
				f[0][k].z += 1.0*(ds[0]*ds[1]*pg*(R1.z-R2.z)*Bt[k][0]/d  +  pe*dR1.z*dBt[k][0]*ds[1]/ds[0]);
				
				f[1][k].x += 1.0*(ds[0]*ds[1]*pg*(R2.x-R1.x)*Bt[k][1]/d  +  pe*dR2.x*dBt[k][1]*ds[0]/ds[1]);
				f[1][k].y += 1.0*(ds[0]*ds[1]*pg*(R2.y-R1.y)*Bt[k][1]/d  +  pe*dR2.y*dBt[k][1]*ds[0]/ds[1]);
				f[1][k].z += 1.0*(ds[0]*ds[1]*pg*(R2.z-R1.z)*Bt[k][1]/d  +  pe*dR2.z*dBt[k][1]*ds[0]/ds[1]);
			}
		}
	}
	
	//printf("Printing orig forces\n");
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",b1.f[k].x,b1.f[k].y,b1.f[k].z);
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",b2.f[k].x,b2.f[k].y,b2.f[k].z);
		
	double factor = -1.0*h*h/9.0;
	//printf("Printing inter forces\n");
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",f[0][k].x*factor,f[0][k].y*factor,f[0][k].z*factor);
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",f[1][k].x*factor,f[1][k].y*factor,f[1][k].z*factor);
		
	//for(int k=0;k<4;k++) {
		//b1.f[k].x += f[0][k].x*factor;	b1.f[k].y += f[0][k].y*factor; b1.f[k].z += f[0][k].z*factor;
		//b2.f[k].x += f[1][k].x*factor;	b2.f[k].y += f[1][k].y*factor; b2.f[k].z += f[1][k].z*factor;
	//}
	
	for(int k=0;k<4;k++) {
		bz[0].f[k].x = f[0][k].x*factor;	bz[0].f[k].y = f[0][k].y*factor; bz[0].f[k].z = f[0][k].z*factor;
		bz[1].f[k].x = f[1][k].x*factor;	bz[1].f[k].y = f[1][k].y*factor; bz[1].f[k].z = f[1][k].z*factor;
	}
	
	//printf("Printing modf. forces\n");
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",b1.f[k].x,b1.f[k].y,b1.f[k].z);
	//for(int k=0;k<4;k++) printf("%lf %lf %lf\n",b2.f[k].x,b2.f[k].y,b2.f[k].z);
	
	return bz;
}


Bezier fire_min0(Bezier b0, FILE *fcom) {
	Bezier bk;
	
	double dt = 0.001; int count=0;
	double alpha0 = 0.1,alpha = alpha0, fdec = 0.5, finc = 1.1, dtmax = 1.0, falpha = 0.99;
	b0.intra_force(); 
	
	copy(b0.p,bk.p); copy(b0.f,bk.f); 
	bk.compute_mass();
	
	for(int j=0;j<4;j++) { bk.v[j].x = 0.0; bk.v[j].y = 0.0; bk.v[j].z = 0.0; }
	
	do {
		double P = dot(bk.f,bk.v);
		//printf("\n P %lf",P);
		if(P<0) {
			for(int j=0;j<4;j++) { bk.v[j].x = 0.0; bk.v[j].y = 0.0; bk.v[j].z = 0.0; }
			dt *= fdec; alpha = alpha0;
		}else {
			double cf = alpha*sqrt(dot(bk.v,bk.v)/dot(bk.f,bk.f));
			double cv = 1.0 - alpha;
			if( dot(bk.f,bk.f) == 0) cf = 0.0;
			//printf("\n P cf cv %lf %lf %lf",P,cf,cv);
			for(int j=0;j<4;j++) { 
				bk.v[j].x = cv*bk.v[j].x + cf*bk.f[j].x; 
				bk.v[j].y = cv*bk.v[j].y + cf*bk.f[j].y; 
				bk.v[j].z = cv*bk.v[j].z + cf*bk.f[j].z; 
			}
			dt = std::min(dt*finc,dtmax); alpha *= falpha;
		}
		
		bk.integrate_verlet(dt);
		bk.display_curve(count,fcom);
		////printf("\n	Iteration prior %d %lf %lf",count,value(xk),mag(xk));
		//bk1 = integrate_verlet(bk,dt); //Integration via velocity verlet
		//xk = Xk1.x; vk = Xk1.v; fk = Xk1.f;
		printf("Iteration %d %.12f %.12f\n",count,bk.intra_energy(),norm(bk.f));
		count++;
	}while(norm(bk.f)>TOL && count<MAX);
	
	return bk;
}

vecBez fire_min(vecBez b0, FILE *fcom, FILE *fcom2) {
	vecBez bk(Nbez),bz(2);
	
	double dt = 1.0e-4; int count=0;
	double alpha0 = 0.1,alpha = alpha0, fdec = 0.5, finc = 1.1, dtmax = 1e-2, falpha = 0.99, total_norm;
	
	for(int i=0;i<Nbez;i++) b0[i].intra_force();
		
	for(int i=0;i<Nbez;i++) {
		copy(b0[i].p,bk[i].p); copy(b0[i].f,bk[i].f); copy(b0[i].cons,bk[i].cons);
		bk[i].compute_mass();
		for(int j=0;j<4;j++) { bk[i].v[j].x = 0.0; bk[i].v[j].y = 0.0; bk[i].v[j].z = 0.0; }
	}
	
	bz = inter_force(bk[0],bk[1]);
	for(int j=0;j<4;j++) {  
		bk[0].f[j].x += bz[0].f[j].x; bk[0].f[j].y += bz[0].f[j].y; bk[0].f[j].z += bz[0].f[j].z;
		bk[1].f[j].x += bz[1].f[j].x; bk[1].f[j].y += bz[1].f[j].y; bk[1].f[j].z += bz[1].f[j].z;
	}	
	
	do {
		double P = 0; for(int i=0;i<Nbez;i++) P += dot(bk[i].f,bk[i].v);
		//printf("\n P %lf",P);
		if(P<0) {
			for(int i=0;i<Nbez;i++) for(int j=0;j<4;j++) { bk[i].v[j].x = 0.0; bk[i].v[j].y = 0.0; bk[i].v[j].z = 0.0; }
			dt *= fdec; alpha = alpha0;
		} else {
			double vv=0.0,ff=0.0;
			for(int i=0;i<Nbez;i++) { vv += dot(bk[i].v,bk[i].v); ff += dot(bk[i].f,bk[i].f); }
			double cf = alpha*sqrt(vv/ff);
			double cv = 1.0 - alpha;
			if( ff == 0) cf = 0.0;
			//printf(" P cf cv %lf %lf %lf\n",P,cf,cv);
			for(int i=0;i<Nbez;i++) {
				for(int j=0;j<4;j++) { 
					bk[i].v[j].x = cv*bk[i].v[j].x + cf*bk[i].f[j].x; 
					bk[i].v[j].y = cv*bk[i].v[j].y + cf*bk[i].f[j].y; 
					bk[i].v[j].z = cv*bk[i].v[j].z + cf*bk[i].f[j].z; 
				}
			}
			dt = std::min(dt*finc,dtmax); alpha *= falpha;
		}
		
		for(int i=0;i<Nbez;i++)
		{ 
			bk[i].compute_mass();
			bk[i].integrate_verlet1(dt);
		}
		
		for(int i=0;i<Nbez;i++) bk[i].intra_force();
		
		//printf("\n\nFire_min: After intra_force\n");
		//for(int i=0;i<Nbez;i++) {printpoint(bk[i].f[0]); printpoint(bk[i].f[1]); printpoint(bk[i].f[2]); printpoint(bk[i].f[3]);}
		
		bz = inter_force(bk[0],bk[1]);
		
		for(int j=0;j<4;j++) {  
			bk[0].f[j].x += bz[0].f[j].x; bk[0].f[j].y += bz[0].f[j].y; bk[0].f[j].z += bz[0].f[j].z;
			bk[1].f[j].x += bz[1].f[j].x; bk[1].f[j].y += bz[1].f[j].y; bk[1].f[j].z += bz[1].f[j].z;
		}		
		
		//printf("Fire_min: After inter_force\n");
		//for(int i=0;i<Nbez;i++) {printpoint(bk[i].f[0]); printpoint(bk[i].f[1]); printpoint(bk[i].f[2]); printpoint(bk[i].f[3]);}
		
		for(int i=0;i<Nbez;i++) bk[i].integrate_verlet2(dt);
		
		if(count%2==0) display_allcurves(bk,count,fcom);
		if(count%2==0) display_all_cps(bk,count,fcom2);
		
		double total_energy=0.0;total_norm=0.0;
		for(int i=0;i<Nbez;i++) {
			total_energy += bk[i].intra_energy();
			total_norm += dot(bk[i].f,bk[i].f);
		}
		
		double interE = inter_energy(bk[0],bk[1]);
		total_energy += interE;
		total_norm = sqrt(total_norm);
		if(count%2==0) printf("Iteration %d %lf %.12f %.12f %.12f\n",count,dt,total_energy,total_norm,interE);
		count++;
	}while(total_norm>TOL && count<MAX);
	//while(norm(bk[0].f)>TOL && count<MAX);
	
	return bk;
}
