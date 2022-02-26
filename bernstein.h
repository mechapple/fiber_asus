
double Bern0(double t) { return pow(1.0-t,3); }
double dBern0(double t) { return -3.0*pow(1.0-t,2); }
double ddBern0(double t) { return 6.0*(1.0-t); }

double Bern1(double t) { return 3.0*t*pow(1.0-t,2); }
double dBern1(double t) { return 3.0*(1.0-t)*(1.0-3*t); }
double ddBern1(double t) { return 6.0*(3.0*t-2); }

double Bern2(double t) { return 3.0*pow(t,2)*(1.0-t); }
double dBern2(double t) { return 3.0*t*(2.0-3*t); }
double ddBern2(double t) { return 6.0*(1-3.0*t); }

double Bern3(double t) { return pow(t,3); }
double dBern3(double t) { return 3.0*pow(t,2); }
double ddBern3(double t) { return 6.0*t; }
