class Bezier {
		
	public: 
		int init;
		double orig_length;
		point p[4],pd[4],v[4],a[4],a1[4],f[4],fax[4],fbx[4],cons[4],p0[4];
		Eigen::MatrixXf mass,mass_inv,mass2;
		point r(double);
		point dr(double);
		point ddr(double);
		double length();
		double axial_energy();
		double bending_energy();
		void axial_force();
		void bending_force();
		void compute_mass();
		void intra_force();
		double intra_energy();
		void integrate_verlet(double);
		void integrate_verlet1(double);
		void integrate_verlet2(double);
		void display_curve(int,FILE *);
		Bezier();
};

Bezier :: Bezier() {
	init = 0;
	for (int i=0;i<4;i++) { cons[i].x = 1.0; cons[i].y = 1.0; cons[i].z = 1.0;}
}

point Bezier::r(double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*B[i](t);
		pt.y += p[i].y*B[i](t);
		pt.z += p[i].z*B[i](t);
	}
	return pt;
}

point Bezier::dr(double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*dB[i](t);
		pt.y += p[i].y*dB[i](t);
		pt.z += p[i].z*dB[i](t);
	}
	return pt;
}

point Bezier::ddr(double t) {
	point pt = { .x=0.0, .y=0.0, .z=0.0 };
	for (int i=0;i<4;i++) {
		pt.x += p[i].x*ddB[i](t);
		pt.y += p[i].y*ddB[i](t);
		pt.z += p[i].z*ddB[i](t);
	}
	return pt;
}

double Bezier::length() {
	double xmin[1] = {0}, xmax[1] = {1}, val[1], err[1];
	hcubature(1, f_Length, p, 1, xmin, xmax, 0, 0, TOL1, ERROR_INDIVIDUAL, val, err);
	return val[0];
}

double Bezier::axial_energy(){ 
	double L = length();
	
	if(init!=0) return 0.5*EA*pow(L-orig_length,2);
	else {
		orig_length = L; init=1.0;
		return 0.0;
	}
}

void Bezier::axial_force() {
	double xmin[1] = {0}, xmax[1] = {1}, val[12], err[12];
	hcubature(12, f_AxForce, p, 1, xmin, xmax, 0, 0, TOL1, ERROR_INDIVIDUAL, val, err);
	//printf("fax_curbature\n");
	//for(int i=0;i<3;i++) {for(int j=0;j<4;j++) printf("%lf ",val[i*4+j]); printf("\n");}
		
	double factor = -1.0*EA*(length()-orig_length);
	for(int i=0;i<3;i++) for(int j=0;j<4;j++) val[i*4+j] *= factor;
	
	fax[0].x = val[0]; fax[0].y = val[4]; fax[0].z = val[8];
	fax[1].x = val[1]; fax[1].y = val[5]; fax[1].z = val[9];
	fax[2].x = val[2]; fax[2].y = val[6]; fax[2].z = val[10];
	fax[3].x = val[3]; fax[3].y = val[7]; fax[3].z = val[11];
}

double Bezier::bending_energy() { 
	double xmin[1] = {0}, xmax[1] = {1}, val[1], err[1];
	hcubature(1, f_BeEnergy, p, 1, xmin, xmax, 0, 0, TOL1, ERROR_INDIVIDUAL, val, err);
	return val[0]*0.5*EI;
}

void Bezier::bending_force() {
	double xmin[1] = {0}, xmax[1] = {1}, val[12], err[12];
	hcubature(12, f_BeForce, p, 1, xmin, xmax, 0, 0, TOL1, ERROR_INDIVIDUAL, val, err);
	//printf("fbx_curbature\n");
	//for(int i=0;i<3;i++) {for(int j=0;j<4;j++) printf("%.15f ",val[i*4+j]); printf("\n");}
	
	double factor = -0.5*EI;
	for(int i=0;i<3;i++) for(int j=0;j<4;j++) val[i*4+j] *= factor;
	
	fbx[0].x = val[0]; fbx[0].y = val[4]; fbx[0].z = val[8];
	fbx[1].x = val[1]; fbx[1].y = val[5]; fbx[1].z = val[9];
	fbx[2].x = val[2]; fbx[2].y = val[6]; fbx[2].z = val[10];
	fbx[3].x = val[3]; fbx[3].y = val[7]; fbx[3].z = val[11];
}

void Bezier::intra_force() {
	bending_force();
	axial_force();
	for(int j=0;j<4;j++) {
		f[j].x = fax[j].x + fbx[j].x; 	f[j].y = fax[j].y + fbx[j].y;	f[j].z = fax[j].z + fbx[j].z;
	}
}

double Bezier::intra_energy() {
	double Eb = bending_energy();
	double Ef = axial_energy();
	return (Eb+Ef);
}

void Bezier::compute_mass() {
	mass = Eigen::MatrixXf::Zero(4,4);
	double xmin[1] = {0}, xmax[1] = {1}, val[16], err[16];
	hcubature(16, f_Mass, p, 1, xmin, xmax, 0, 0, TOL1, ERROR_INDIVIDUAL, val, err);
	for(int j=0;j<4;j++) for(int k=0;k<4;k++) mass(j,k) = rho*val[j*4+k];
	
	mass_inv = mass.inverse();
	//std::cout << "Here is the curbature mass matrix :\n" << mass << std::endl;
	//std::cout << "Here is the inverse mass matrix :\n" << mass_inv << std::endl;
}

void Bezier::display_curve(int time, FILE *fcom) {
	double min=-100.0,max=100.0;
	fprintf(fcom,"ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n%d\nITEM: BOX BOUNDS pp pp pp",time,Nres+1);
	fprintf(fcom,"\n%.3f %.3f \n%.3f %.3f \n%.3f %.3f \nITEM: ATOMS id mol type q xu yu zu\n",min,max,min,max,min,max);
	double h = 1.0/Nres;
	for(int i=0; i<=Nres; i++) {
		double t = i*h; point R = r(t);
		fprintf(fcom,"%d 1 1 0.00 %lf %lf %lf\n",i+1,R.x,R.y,R.z);
	}
	
}

void copy(point src[], point dest[]) {
	for(int j=0;j<4;j++) {
		dest[j].x = src[j].x;	dest[j].y = src[j].y;	dest[j].z = src[j].z;
	}
}

double norm(point src[]) {
	double sum=0.0;
	for(int j=0;j<4;j++) {
		sum += src[j].x*src[j].x;
		sum += src[j].y*src[j].y;
		sum += src[j].z*src[j].z;
	}
	return sqrt(sum);
}

double dot(point p1[],point p2[]) {
	double sum=0.0;
	for(int j=0;j<4;j++) {
		sum += p1[j].x*p2[j].x;
		sum += p1[j].y*p2[j].y;
		sum += p1[j].z*p2[j].z;
	}
	return sum;
}

void Bezier::integrate_verlet(double dt) {
	Eigen::VectorXf force(4),accl(4);
	force << f[0].x, f[1].x, f[2].x, f[3].x; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].x = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	force << f[0].y, f[1].y, f[2].y, f[3].y; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].y = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	force << f[0].z, f[1].z, f[2].z, f[3].z; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].z = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	for(int j=0;j<4;j++) {
		p[j].x += (v[j].x*dt + a[j].x*dt*dt*0.5)*cons[j].x;
		p[j].y += (v[j].y*dt + a[j].y*dt*dt*0.5)*cons[j].y;
		p[j].z += (v[j].z*dt + a[j].z*dt*dt*0.5)*cons[j].z;
	}
	
	intra_force();
	
	Eigen::VectorXf force1(4),accl1(4);
	force1 << f[0].x, f[1].x, f[2].x, f[3].x; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].x = accl1[j];
	
	force1 << f[0].y, f[1].y, f[2].y, f[3].y; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].y = accl1[j];
	
	force1 << f[0].z, f[1].z, f[2].z, f[3].z; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].z = accl1[j];
	
	for(int j=0;j<4;j++) {
		v[j].x += (a[j].x+a1[j].x)*dt*0.5*cons[j].x;
		v[j].y += (a[j].y+a1[j].y)*dt*0.5*cons[j].y;
		v[j].z += (a[j].z+a1[j].z)*dt*0.5*cons[j].z;
	}
}

void Bezier::integrate_verlet1(double dt) {
	Eigen::VectorXf force(4),accl(4);
	force << f[0].x, f[1].x, f[2].x, f[3].x; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].x = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	force << f[0].y, f[1].y, f[2].y, f[3].y; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].y = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	force << f[0].z, f[1].z, f[2].z, f[3].z; accl = mass_inv*force;
	for(int j=0;j<4;j++) a[j].z = accl[j];
	//std::cout << "Here is the accl vector :\n" << accl << std::endl;
	
	for(int j=0;j<4;j++) {
		p[j].x += (v[j].x*dt + a[j].x*dt*dt*0.5)*cons[j].x;
		p[j].y += (v[j].y*dt + a[j].y*dt*dt*0.5)*cons[j].y;
		p[j].z += (v[j].z*dt + a[j].z*dt*dt*0.5)*cons[j].z;
	}
}

void Bezier::integrate_verlet2(double dt) {
	Eigen::VectorXf force1(4),accl1(4);
	force1 << f[0].x, f[1].x, f[2].x, f[3].x; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].x = accl1[j];
	
	force1 << f[0].y, f[1].y, f[2].y, f[3].y; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].y = accl1[j];
	
	force1 << f[0].z, f[1].z, f[2].z, f[3].z; accl1 = mass_inv*force1;
	for(int j=0;j<4;j++) a1[j].z = accl1[j];
	
	for(int j=0;j<4;j++) {
		v[j].x += (a[j].x+a1[j].x)*dt*0.5*cons[j].x;
		v[j].y += (a[j].y+a1[j].y)*dt*0.5*cons[j].y;
		v[j].z += (a[j].z+a1[j].z)*dt*0.5*cons[j].z;
	}
}
