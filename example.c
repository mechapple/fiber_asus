#include <stdio.h>
#include <math.h>
#include "cubature.h"

int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    double sigma[2] = {((double *) fdata)[0],((double *) fdata)[1]}; // we can pass σ via fdata argument
    double sum = 0;
    unsigned i;
    for (i = 0; i < ndim; ++i) sum += x[i] * x[i];
    // compute the output value: note that fdim should == 1 from below
    fval[0] = exp(-sigma[0] * sum);
    fval[1] = exp(-sigma[1] * 2 * sum);
    return 0; // success
}

int main(int argc, char *argv[]) 
{
    //double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma[2] = {0.5,0.25}, val[2], err[2];
    double xmin[1] = {-2}, xmax[1] = {2}, sigma[2] = {0.5,0.25}, val[2], err[2];
    hcubature(2, f, sigma, 1, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, val, err);
    printf("Computed integral = %0.10g +/- %g\n", val[0], err[0]);
    printf("Computed integral = %0.10g +/- %g\n", val[1], err[1]);
    //hcubature_v(1, f1, &sigma, 3, xmin, xmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &val, &err);
    return 0;
}
