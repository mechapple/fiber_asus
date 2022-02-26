#define sig 0.01
#define eps 1.0e-2
#define rcut 1.05

double LJ_energy(double r) {return 4*eps*(pow(sig/r,12) - pow(sig/r,6)); }
double LJ_force(double r) {return 48*(eps/sig)*( pow(sig/r,13) - 0.5*pow(sig/r,7) ); }

const double lj1 = LJ_force(rcut);
const double lj2 = -rcut*lj1 - LJ_energy(rcut);

double pair_energy(double r) {
	if(r>rcut) return 0;
	double rn = pow(sig/r,6);
	//double rn = pow(sig/r,3);
	
	//return 4*eps*(r-sig)*(r-sig);
	return (4*eps*(rn*rn-rn) + lj1*r + lj2); //shifted LJ
}

double pair_force(double r) {
	if(r>rcut) return 0;
	double rn = pow(sig/r,6);
	return (48*eps*(rn*rn-0.5*rn)/r - lj1); //shifted LJ
}

double pair_gradient(double r) {
	if(r>rcut) return 0;
	double rn = pow(sig/r,6);
	return (lj1 - 48*eps*(rn*rn-0.5*rn)/r); //shifted LJ
}
