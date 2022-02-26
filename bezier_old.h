class Bezier {
		
	public: 
		int init;
		double orig_length;
		point p[4],pd[4],v[4],a[4],a1[4],f[4],fax[4],fbx[4],cons[4];
		Eigen::MatrixXf mass,mass_inv;
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
	double sum = 0.0,h = 1.0/Nsimp;
	int m[3] = {2,3,3};
	for(int i=0; i<=Nsimp; i++) sum += normpoint(dr(i*h))*m[i%3]; //Composite Simpson's 3/8 rule
	sum -= (normpoint(dr(0.0))+normpoint(dr(1.0))); //correct for end points
	
	return sum*(3*h/8);
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
	double sum[3][4]={0},h = 1.0/Nsimp;
	int m[3] = {2,3,3};
	for(int i=0; i<=Nsimp; i++) {
		double t = i*h,ds = normpoint(dr(t)); point dR = dr(t);
		for(int j=0;j<4;j++) sum[0][j] += m[i%3]*(1.0/ds)*dR.x*dB[j](t);
		for(int j=0;j<4;j++) sum[1][j] += m[i%3]*(1.0/ds)*dR.y*dB[j](t);
		for(int j=0;j<4;j++) sum[2][j] += m[i%3]*(1.0/ds)*dR.z*dB[j](t);
	}
	for(double t=0; t<=1.0; t+=1.0) { //correct for end points 0 and 1
		double ds = normpoint(dr(t)); point dR = dr(t);
		for(int j=0;j<4;j++) sum[0][j] -= (1.0/ds)*dR.x*dB[j](t);
		for(int j=0;j<4;j++) sum[1][j] -= (1.0/ds)*dR.y*dB[j](t);
		for(int j=0;j<4;j++) sum[2][j] -= (1.0/ds)*dR.z*dB[j](t);
	}
	
	double factor = -1.0*EA*(length()-orig_length)*(3*h/8);
	for(int i=0;i<3;i++) for(int j=0;j<4;j++) sum[i][j] *= factor;

	//printf("fax\n");
	//for(int i=0;i<3;i++) {for(int j=0;j<4;j++) printf("%lf ",sum[i][j]); printf("\n");}
	
	fax[0].x = sum[0][0]; fax[0].y = sum[1][0]; fax[0].z = sum[2][0];
	fax[1].x = sum[0][1]; fax[1].y = sum[1][1]; fax[1].z = sum[2][1];
	fax[2].x = sum[0][2]; fax[2].y = sum[1][2]; fax[2].z = sum[2][2];
	fax[3].x = sum[0][3]; fax[3].y = sum[1][3]; fax[3].z = sum[2][3];
}

double Bezier::bending_energy() { 
	double sum = 0.0,h = 1.0/Nsimp;
	int m[3] = {2,3,3};
	for(int i=0; i<=Nsimp; i++) {
		double t = i*h;
		double num = pow((ddr(t).z*dr(t).y - ddr(t).y*dr(t).z),2) + pow((ddr(t).x*dr(t).z - ddr(t).z*dr(t).x),2) + pow((ddr(t).y*dr(t).x - ddr(t).x*dr(t).y),2);
		double den = pow(normpoint(dr(t)),5);
		sum += m[i%3]*num/den;  //Composite Simpson's 3/8 rule
	}
	
	for(double t=0; t<=1.0; t+=1.0) { //correct for end points 0 and 1
		double num = pow((ddr(t).z*dr(t).y - ddr(t).y*dr(t).z),2) + pow((ddr(t).x*dr(t).z - ddr(t).z*dr(t).x),2) + pow((ddr(t).y*dr(t).x - ddr(t).x*dr(t).y),2);
		double den = pow(normpoint(dr(t)),5);
		sum -= num/den;
	}
	return sum*(3*h/8)*0.5*EI;
}

void Bezier::bending_force() {
	
	double sum[3][4]={0},h = 1.0/Nsimp;
	int m[3] = {2,3,3};
	for(int i=0; i<=Nsimp; i++) { //Composite Simpson's 3/8 rule
		double t = i*h; point dR = dr(t), ddR = ddr(t);
		double num = pow((ddR.z*dR.y - ddR.y*dR.z),2) + pow((ddR.x*dR.z - ddR.z*dR.x),2) + pow((ddR.y*dR.x - ddR.x*dR.y),2);
		double Nm[3] = {(ddR.z*dR.y - ddR.y*dR.z), (ddR.x*dR.z - ddR.z*dR.x), (ddR.y*dR.x - ddR.x*dR.y)};
		double ds = normpoint(dR); 
		for(int j=0;j<4;j++) {
			double dNdxj = 2*Nm[1]*(ddB[j](t)*dR.z - ddR.z*dB[j](t)) + 2*Nm[2]*(ddR.y*dB[j](t) - ddB[j](t)*dR.y);
			double dsdxj = (dR.x/ds)*dB[j](t);
			sum[0][j] += m[i%3]*(ds*dNdxj - 5*num*dsdxj)/pow(ds,6);
		}
		
		for(int j=0;j<4;j++) {
			double dNdyj = 2*Nm[0]*(ddR.z*dB[j](t) - ddB[j](t)*dR.z) + 2*Nm[2]*(ddB[j](t)*dR.x - ddR.x*dB[j](t));
			double dsdyj = (dR.y/ds)*dB[j](t);
			sum[1][j] += m[i%3]*(ds*dNdyj - 5*num*dsdyj)/pow(ds,6);
		}
		
		for(int j=0;j<4;j++) {
			double dNdzj = 2*Nm[0]*(ddB[j](t)*dR.y - ddR.y*dB[j](t)) + 2*Nm[1]*(ddR.x*dB[j](t) - ddB[j](t)*dR.x);
			double dsdzj = (dR.z/ds)*dB[j](t);
			sum[2][j] += m[i%3]*(ds*dNdzj - 5*num*dsdzj)/pow(ds,6);
		}
	}
	
	for(double t=0; t<=1.0; t+=1.0) { //correct for end points 0 and 1
		point dR = dr(t), ddR = ddr(t);
		double num = pow((ddR.z*dR.y - ddR.y*dR.z),2) + pow((ddR.x*dR.z - ddR.z*dR.x),2) + pow((ddR.y*dR.x - ddR.x*dR.y),2);
		double Nm[3] = {(ddR.z*dR.y - ddR.y*dR.z), (ddR.x*dR.z - ddR.z*dR.x), (ddR.y*dR.x - ddR.x*dR.y)};
		double ds = normpoint(dR); 
		for(int j=0;j<4;j++) {
			double dNdxj = 2*Nm[1]*(ddB[j](t)*dR.z - ddR.z*dB[j](t)) + 2*Nm[2]*(ddR.y*dB[j](t) - ddB[j](t)*dR.y);
			double dsdxj = (dR.x/ds)*dB[j](t);
			sum[0][j] -= 1.0*(ds*dNdxj - 5*num*dsdxj)/pow(ds,6);
		}
		
		for(int j=0;j<4;j++) {
			double dNdyj = 2*Nm[0]*(ddR.z*dB[j](t) - ddB[j](t)*dR.z) + 2*Nm[2]*(ddB[j](t)*dR.x - ddR.x*dB[j](t));
			double dsdyj = (dR.y/ds)*dB[j](t);
			sum[1][j] -= 1.0*(ds*dNdyj - 5*num*dsdyj)/pow(ds,6);
		}
		
		for(int j=0;j<4;j++) {
			double dNdzj = 2*Nm[0]*(ddB[j](t)*dR.y - ddR.y*dB[j](t)) + 2*Nm[1]*(ddR.x*dB[j](t) - ddB[j](t)*dR.x);
			double dsdzj = (dR.z/ds)*dB[j](t);
			sum[2][j] -= 1.0*(ds*dNdzj - 5*num*dsdzj)/pow(ds,6);
		}
	}
	
	double factor = -0.5*EI*(3*h/8);
	for(int i=0;i<3;i++) for(int j=0;j<4;j++) sum[i][j] *= factor;
	
	//printf("fbx\n");
	//for(int i=0;i<3;i++) {for(int j=0;j<4;j++) printf("%lf ",sum[i][j]); printf("\n");}
	
	fbx[0].x = sum[0][0]; fbx[0].y = sum[1][0]; fbx[0].z = sum[2][0];
	fbx[1].x = sum[0][1]; fbx[1].y = sum[1][1]; fbx[1].z = sum[2][1];
	fbx[2].x = sum[0][2]; fbx[2].y = sum[1][2]; fbx[2].z = sum[2][2];
	fbx[3].x = sum[0][3]; fbx[3].y = sum[1][3]; fbx[3].z = sum[2][3];
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
	double sum = 0.0,h = 1.0/Nsimp;
	int m[3] = {2,3,3};
	//for(int j=0;j<4;j++) printpoint(p[j]);
	for(int i=0; i<=Nsimp; i++) {
		double t = i*h; point dR = dr(t);
		double ds = normpoint(dR);
		double Bt[4]; for(int j=0;j<4;j++) Bt[j] = B[j](t);
		for(int j=0;j<4;j++) for(int k=0;k<4;k++) {
			mass(j,k) += m[i%3]*ds*Bt[j]*Bt[k];
		}
	}
	
	for(double t=0; t<=1.0; t+=1.0) { //correct for end points 0 and 1
		point dR = dr(t); double ds = normpoint(dR);
		double Bt[4]; for(int j=0;j<4;j++) Bt[j] = B[j](t);
		for(int j=0;j<4;j++) for(int k=0;k<4;k++) {
			mass(j,k) -= ds*Bt[j]*Bt[k];
		}
	}
	
	double factor = rho*(3*h/8);
	for(int j=0;j<4;j++) for(int k=0;k<4;k++) mass(j,k) *= factor;
	
	mass_inv = mass.inverse();
	//std::cout << "Here is the mass matrix :\n" << mass << std::endl;
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
