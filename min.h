#define Nsimp2D 1e3
#define xscale 50.0

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

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
	double tmin = 0.0,tmax=1.0;
	double xmin[2] = {tmin,tmin}, xmax[2] = {tmax,tmax}, val[1], err[1];
	point rp[8] = {b1.p[0],b1.p[1],b1.p[2],b1.p[3],b2.p[0],b2.p[1],b2.p[2],b2.p[3]};
	
	hcubature(1, f_InterEnergy, rp, 2, xmin, xmax, maxEval, ABSTOL, TOL2, ERROR_INDIVIDUAL, val, err);
	printf("Computed integral_e = %0.15f +/- %.15g\n", val[0], err[0]);
	return val[0];
}

double constrain(Bezier b1,Bezier b2) {
	double cE1 = pow(b1.p[0].x,2)+pow(b1.p[0].y,2)+pow(b1.p[0].z,2);
	double cE2 = pow(b2.p[0].x-2.0*xfac,2)+pow(b2.p[0].y,2)+pow(b2.p[0].z,2);
	
	return 100.0*(cE1 + cE2);
}

double inter_area(Bezier b1,Bezier b2) {
	double tmin = 0.0,tmax=1.0;
	double xmin[2] = {tmin,tmin}, xmax[2] = {tmax,tmax}, val[1], err[1];
	point rp[8] = {b1.p[0],b1.p[1],b1.p[2],b1.p[3],b2.p[0],b2.p[1],b2.p[2],b2.p[3]};
	
	hcubature(1, f_InterArea, rp, 2, xmin, xmax, maxEval, ABSTOL, TOL2, ERROR_INDIVIDUAL, val, err);
	printf("Computed area_e = %0.15f +/- %.15g\n", val[0], err[0]);
	return val[0];
}

vecBez inter_force(Bezier b1,Bezier b2) {
	vecBez bz(2);
	double tmin = 0.0,tmax=1.0;
	double xmin[2] = {tmin,tmin}, xmax[2] = {tmax,tmax}, val[24], err[24];
	
	point rp[8] = {b1.p[0],b1.p[1],b1.p[2],b1.p[3],b2.p[0],b2.p[1],b2.p[2],b2.p[3]};
	hcubature(24, f_InterForce, rp, 2, xmin, xmax, maxEval, ABSTOL, TOL3, ERROR_INDIVIDUAL, val, err);
	printf("Computed integral_f0 = %0.15f +/- %.15g\n", val[0], err[0]);
	
	//printf("Printing inter cubature forces\n");
	//for(int k=0;k<4;k++) printf("%.15f %.15f %.15f\n",val[3*k],val[3*k+1],val[3*k+2]);
	//for(int k=0;k<4;k++) printf("%.15f %.15f %.15f\n",val[3*k+12],val[3*k+13],val[3*k+14]);

	for(int k=0;k<4;k++) {
		bz[0].f[k].x = val[3*k];	bz[0].f[k].y = val[3*k+1]; bz[0].f[k].z = val[3*k+2];
		bz[1].f[k].x = val[3*k+12];	bz[1].f[k].y = val[3*k+13]; bz[1].f[k].z = val[3*k+14];
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
	double alpha0 = 0.1,alpha = alpha0, fdec = 0.5, finc = 1.1, dtmax = 1e-3, falpha = 0.99, total_norm;
	
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

double my_f (const gsl_vector *v, void *params)
{
	vecBez bk(Nbez);
	double *p = (double *)params;
	
	int count_x = 0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0;j<4;j++) { 
			bk[i].p[j].x = gsl_vector_get (v, count_x);
			bk[i].p[j].y = gsl_vector_get (v, count_x+1);
			bk[i].p[j].z = gsl_vector_get (v, count_x+2);
			count_x += 3;
		}
	}
	
	bk[0].orig_length = p[0]; bk[0].init = 1;
	bk[1].orig_length = p[1]; bk[1].init = 1;
	
	double larea = inter_area(bk[0],bk[1]);
	
	double sum = bk[0].intra_energy() + bk[1].intra_energy() + 1.0*inter_energy(bk[0],bk[1]) + constrain(bk[0],bk[1]);
	return sum;
}

/* The gradient of f, df = (df/dx, df/dy). */
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
	vecBez bk(Nbez),bz;
	double *p = (double *)params;
	
	int count_x = 0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0;j<4;j++) { 
			bk[i].p[j].x = gsl_vector_get (v, count_x);
			bk[i].p[j].y = gsl_vector_get (v, count_x+1);
			bk[i].p[j].z = gsl_vector_get (v, count_x+2);
			count_x += 3;
		}
	}
	
	bk[0].orig_length = p[0]; bk[0].init = 1;
	bk[1].orig_length = p[1]; bk[1].init = 1;
	
	for(int i=0;i<Nbez;i++) bk[i].intra_force();
	bz = inter_force(bk[0],bk[1]);
	
	for(int j=0;j<4;j++) {  
		bk[0].f[j].x += bz[0].f[j].x; bk[0].f[j].y += bz[0].f[j].y; bk[0].f[j].z += bz[0].f[j].z;
		bk[1].f[j].x += bz[1].f[j].x; bk[1].f[j].y += bz[1].f[j].y; bk[1].f[j].z += bz[1].f[j].z;
	}
	
	count_x = 0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0;j<4;j++) { 
			
			gsl_vector_set(df, count_x, -bk[i].f[j].x);
			gsl_vector_set(df, count_x+1, -bk[i].f[j].y);
			gsl_vector_set(df, count_x+2, -bk[i].f[j].z);
			
			count_x += 3;
		}
	}
		
	//gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
	//gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
}

/* Compute both f and df together. */
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
	*f = my_f(x, params); 
	my_df(x, params, df);
}

vecBez bfgs2_min(vecBez b0, FILE *fcom, FILE *fcom2, FILE *flog) 
{	
	vecBez bk(Nbez);
	size_t iter = 0;
	int status;
	
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	/* Position of the minimum (1,2), scale factors 
	 10,20, height 30. */
	double par[2] = { b0[0].length(), b0[1].length() };
	
	gsl_vector *x;
	gsl_multimin_function_fdf my_func;
	
	my_func.n = 24;
	my_func.f = my_f;
	my_func.df = my_df;
	my_func.fdf = my_fdf;
	my_func.params = par;
	
	/* Starting point, x = b0 */
	x = gsl_vector_alloc (24);
	
	int count_x = 0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0;j<4;j++) { 
			gsl_vector_set (x, count_x, b0[i].p[j].x);
			gsl_vector_set (x, count_x+1, b0[i].p[j].y);
			gsl_vector_set (x, count_x+2, b0[i].p[j].z);
			count_x += 3;
		}
	}
	
	T = gsl_multimin_fdfminimizer_conjugate_fr;
	//T = gsl_multimin_fdfminimizer_vector_bfgs2;
	s = gsl_multimin_fdfminimizer_alloc (T, 24);
	
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.1, 1e-2);
	
	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		
		if(status) break;
		
		status = gsl_multimin_test_gradient (s->gradient, 1e-4);
		gsl_vector *df = gsl_multimin_fdfminimizer_gradient(s);
		
		if (status == GSL_SUCCESS) printf("Minimum found at:\n");
		
		printf("%5zu %10.5f %10.5f\n", iter, s->f, gsl_blas_dnrm2(df));
		
		//if(iter%1==0) 
		{
			
			count_x = 0;
			for(int i=0;i<Nbez;i++) {
				for(int j=0;j<4;j++) { 
					
					bk[i].p[j].x = gsl_vector_get (s->x, count_x);
					bk[i].p[j].y = gsl_vector_get (s->x, count_x+1);
					bk[i].p[j].z = gsl_vector_get (s->x, count_x+2);
					
					count_x += 3;
				}
			}
			
			display_allcurves(bk,iter,fcom);
			display_all_cps(bk,iter,fcom2);
		}
		
	
	}
	while (status == GSL_CONTINUE && iter < ITERMAX);	
	
	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);
	
	return bk;
}

vecBez simplex_min(vecBez b0, FILE *fcom, FILE *fcom2, FILE *flog) 
{	
	vecBez bk(Nbez);
	
	const gsl_multimin_fminimizer_type *T = 
	gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;
	
	/* Position of the minimum (1,2), scale factors 
	 10,20, height 30. */
	double par[2] = { b0[0].length(), b0[1].length() };
	
	/* Starting point, x = b0 */
	x = gsl_vector_alloc (24);
	
	int count_x = 0;
	for(int i=0;i<Nbez;i++) {
		for(int j=0;j<4;j++) { 
			gsl_vector_set (x, count_x, b0[i].p[j].x);
			gsl_vector_set (x, count_x+1, b0[i].p[j].y);
			gsl_vector_set (x, count_x+2, b0[i].p[j].z);
			count_x += 3;
		}
	}
	
	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (24);
	gsl_vector_set_all (ss, 1.0e-2);
	
	/* Initialize method and iterate */
	minex_func.n = 24;
	minex_func.f = my_f;
	minex_func.params = par;
	
	s = gsl_multimin_fminimizer_alloc (T, 24);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate (s);
		
		if(status) break;
		
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-6	);
				
		if (status == GSL_SUCCESS) printf("Minimum found at:\n");
		
		//if(iter%1==0) 
		{
			
			count_x = 0;
			for(int i=0;i<Nbez;i++) {
				for(int j=0;j<4;j++) { 
					
					bk[i].p[j].x = gsl_vector_get (s->x, count_x);
					bk[i].p[j].y = gsl_vector_get (s->x, count_x+1);
					bk[i].p[j].z = gsl_vector_get (s->x, count_x+2);
					
					count_x += 3;
				}
			}
			
			display_allcurves(bk,iter,fcom);
			display_all_cps(bk,iter,fcom2);
		}
		
		printf("%5zu %10.5f %10.5f\n", iter, s->fval, inter_energy(bk[0],bk[1]) );
					
	}
	while (status == GSL_CONTINUE && iter < ITERMAX);	
	
	gsl_vector_free (x);
	gsl_vector_free (ss);
	gsl_multimin_fminimizer_free (s);
	
	return bk;
}
