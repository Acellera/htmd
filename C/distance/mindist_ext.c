#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define X(atom, frame, nframes, nf3) atom*nf3
#define Y(atom, frame, nframes, nf3) atom*nf3+1
#define Z(atom, frame, nframes, nf3) atom*nf3+2

// I redefine here for frames to optimize calculations in the single frame use case
#define Xf(atom, frame, nframes, nf3) atom*nf3+frame
#define Yf(atom, frame, nframes, nf3) atom*nf3+nframes+frame
#define Zf(atom, frame, nframes, nf3) atom*nf3+2*nframes+frame


const int dim = 3;

// Too slow calling a function
//void get_coords(float* coords, int atom, int frame, int nframes, float* c){
//    c[0] = coords[atom*3*nframes+0*nframes+frame];
//    c[1] = coords[atom*3*nframes+1*nframes+frame];
//    c[2] = coords[atom*3*nframes+2*nframes+frame];
//}

void mindist_single_frame(float* coords, int* groups1, int* groups2, int gn1, int gn2, int na, float* dist) {
    float mindist;
    int a,b,g1,g2,n1,n2,g1atm,g2atm;
    float coo1[3], coo2[3];
    int f=0, nf=1;  // current frame and total frames
    int nf3 = nf*3; // Precalculate 3 * nframes for the coordinate lookup macro

    // Iterate over the two group sets
    for (g1=0; g1 < gn1; g1++){
        for (g2=0; g2 < gn2; g2++){
            mindist = -1;
            // Iterate over atoms in the two groups
            for (a = 0; a < na; a++) {
                g1atm = groups1[g1 * na + a];
                if (g1atm == -1) break;
                //get_coords(coords, g1atm, f, nf, coo1);
                //const float* coo1 = &coords[g1atm*dim];

                for (b=0; b < na; b++) {
                    g2atm = groups2[g2 * na + b];
                    if (g2atm == -1) break;
                    //get_coords(coords, g2atm, f, nf, coo2);
                    //const float* coo2 = &coords[g2atm*dim];

                    float d[3];
                    d[0] = coords[X(g1atm,f,nf,nf3)]-coords[X(g2atm,f,nf,nf3)];
                    d[1] = coords[Y(g1atm,f,nf,nf3)]-coords[Y(g2atm,f,nf,nf3)];
                    d[2] = coords[Z(g1atm,f,nf,nf3)]-coords[Z(g2atm,f,nf,nf3)];
                    //printf("coor: %f %f %f / %f %f %f\n", coo1[0], coo1[1], coo1[2], coo2[0], coo2[1], coo2[2]);

                    float D = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                    //printf("index: %d/%d dist: %f\n", g1atm, g2atm, sqrt(D));

                    if (D < mindist || mindist < 0) {
                        mindist = D;
                    }
                }
            }
            //printf("%d %d %d %d: %f\n", g1, gn2, g2, g1*gn2+g2, sqrt(mindist));
            dist[g1*gn2+g2] = sqrt(mindist);
        }
    }
}


void mindist_trajectory(float* coords, float* box, int* groups1, int* groups2, int gn1, int gn2, int na, int nf, int pbc, float* dist) {
    float mindist;
    int a,b,g1,g2,n1,n2,g1atm,g2atm,f;
    float coo1[3], coo2[3];
    int nf3 = nf*3; // Precalculate 3 * nframes for the coordinate lookup macro
    int groupprod = gn1*gn2;

    // Iterate over all frames
    for (f=0; f < nf; f++){
        // Iterate over the two group sets
        for (g1=0; g1 < gn1; g1++){
            for (g2=0; g2 < gn2; g2++){
                mindist = -1;
                // Iterate over atoms in the two groups
                for (a = 0; a < na; a++) {
                    g1atm = groups1[g1 * na + a];
                    if (g1atm == -1) break;
                    //get_coords(coords, g1atm, f, nf, coo1);
                    //const float* coo1 = &coords[g1atm*dim];

                    for (b=0; b < na; b++) {
                        g2atm = groups2[g2 * na + b];
                        if (g2atm == -1) break;
                        //get_coords(coords, g2atm, f, nf, coo2);
                        //const float* coo2 = &coords[g2atm*dim];

                        float d[3];
                        d[0] = coords[Xf(g1atm,f,nf,nf3)]-coords[Xf(g2atm,f,nf,nf3)];
                        d[1] = coords[Yf(g1atm,f,nf,nf3)]-coords[Yf(g2atm,f,nf,nf3)];
                        d[2] = coords[Zf(g1atm,f,nf,nf3)]-coords[Zf(g2atm,f,nf,nf3)];

                        if (pbc){
                            d[0] = d[0] - box[0] * round(d[0] / box[0]);
                            d[1] = d[1] - box[1] * round(d[1] / box[1]);
                            d[2] = d[2] - box[2] * round(d[2] / box[2]);
                        }
                        //printf("coor: %f %f %f / %f %f %f\n", coo1[0], coo1[1], coo1[2], coo2[0], coo2[1], coo2[2]);

                        float D = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
                        //printf("index: %d/%d dist: %f\n", g1atm, g2atm, sqrt(D));

                        if (D < mindist || mindist < 0) {
                            mindist = D;
                        }
                    }
                }
                //printf("%d %d %d %d: %f\n", g1, gn2, g2, g1*gn2+g2, sqrt(mindist));
                dist[g1*gn2+g2+(f*groupprod)] = sqrt(mindist);  // Add here the f
            }
        }
    }
}
