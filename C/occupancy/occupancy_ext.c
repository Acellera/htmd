#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

const int dim = 3;

inline double exp2(double x) {
//ok for inputs less than 5 if precision is not so important
  x = 1.0 + x / 1024;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
}


void occupancy_ext(double* centers1D, float* coords, double* sigmas, double* occus,
                   double* center, double maxlength, int ncenters, int natoms) {
    int a,c;
    double length2 = (maxlength+6)*(maxlength+6);
    for (a=0;a<natoms;a++) {
        const float* coo = &coords[a*dim];
        double d[3];
        d[0] = coo[0]-center[0];
        d[1] = coo[1]-center[1];
        d[2] = coo[2]-center[2];
        double D2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
        if (D2<length2) 
        {
            for (c=0;c<ncenters;c++) if (occus[c]<1.0) {
                const double* cent = &centers1D[c*dim];
                double d[3];
                d[0] = coo[0]-cent[0];
                d[1] = coo[1]-cent[1];
                d[2] = coo[2]-cent[2];
                float D2 = (float) d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                if (D2<25.0f) {//at 5 is already very small
                    float D_1 = 1.0f/sqrt(D2);
                    //float D_1 = rsqrt(D2);
                    double x = sigmas[a]*D_1;
                    double x3= x*x*x;
                    double x12= x3*x3*x3*x3;
                    double value = 1.0 - exp(-x12);
                    occus[c] = MAX(occus[c], value);
                }
            }
        }
    }

}

#define PI (3.141592653589793)
#define SQRT2PI_1  (1.0/sqrt(2.0*PI))

void descriptor_ext(double* centers1D, float* coords, double* sigmas, double* occus,
                   double* center, double maxlength, int ncenters, int natoms, int nchannels) {
    int a,c,h;
    double length2 = (maxlength+6)*(maxlength+6);
    for (a=0;a<natoms;a++) {
        const float* coo = &coords[a*dim];
        const double* atomsigmas = &sigmas[a*nchannels];
        double d[3];
        d[0] = coo[0]-center[0];
        d[1] = coo[1]-center[1];
        d[2] = coo[2]-center[2];
        double D2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
        if (D2<length2) {
            for (c=0;c<ncenters;c++) {
                const double* cent = &centers1D[c*dim];
                double d[3];
                d[0] = coo[0]-cent[0];
                d[1] = coo[1]-cent[1];
                d[2] = coo[2]-cent[2];
                float D2 =  d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                if (D2<25.0) {//at 5 is already very smalli
                    for (h=0; h < nchannels; h++) {
                        //OLD
                    	//double sigma_1 = 1.0/atomsigmas[h];
                   	    // //double value = exp(-0.5*D2*sigma_1*sigma_1); //unnormalized
                   	    //double value = SQRT2PI_1*sigma_1*exp(-0.5*D2*sigma_1*sigma_1); //normalized
                    	//occus[c*nchannels + h] = MAX(occus[c*nchannels + h], value);
                    	//NEW
                    	float D_1 = 1.0f/sqrt(D2);
                        double x = atomsigmas[h]*D_1;
                        double x3= x*x*x;
                        double x12= x3*x3*x3*x3;
                        double value = 1.0 - exp(-x12);
                        occus[c*nchannels + h] = MAX(occus[c*nchannels + h], value);
                    }
                }
            }
        }
    }
}


void descriptor_weight_ext(double* centers1D, float* coords, double* sigma,  double* weights, double* occus,
                   double* center, double maxlength, int ncenters, int natoms, int nchannels) {
    int a,c,h;
    double length2 = (maxlength+6)*(maxlength+6);
    for (a=0;a<natoms;a++) {
        const float* coo = &coords[a*dim];
        const double* atomweights = &weights[a*nchannels];
        double d[3];
        d[0] = coo[0]-center[0];
        d[1] = coo[1]-center[1];
        d[2] = coo[2]-center[2];
        double D2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
        if (D2<length2) {
            for (c=0;c<ncenters;c++) {
                const double* cent = &centers1D[c*dim];
                double d[3];
                d[0] = coo[0]-cent[0];
                d[1] = coo[1]-cent[1];
                d[2] = coo[2]-cent[2];
                float D2 =  d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                if (D2<25.0) {//at 5 is already very smalli
                    for (h=0; h < nchannels; h++) {
                        //OLD
                    	//double sigma_1 = 1.0/atomsigmas[h];
                   	    // //double value = exp(-0.5*D2*sigma_1*sigma_1); //unnormalized
                   	    //double value = SQRT2PI_1*sigma_1*exp(-0.5*D2*sigma_1*sigma_1); //normalized
                    	//occus[c*nchannels + h] = MAX(occus[c*nchannels + h], value);
                    	//NEW
                    	float D_1 = 1.0f/sqrt(D2);
                        double x = sigma[a]*D_1;
                        double x3= x*x*x;
                        double x12= x3*x3*x3*x3;
                        double value = atomweights[h]*(1.0 - exp(-x12));
                        occus[c*nchannels + h] = MAX(occus[c*nchannels + h], value);
                    }
                }
            }
        }
    }

}


