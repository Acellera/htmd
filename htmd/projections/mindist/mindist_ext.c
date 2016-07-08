#include <stdio.h>
#include <stdlib.h>
#include <math.h>


const int dim = 3;

void mindist_ext(float* coords, int* groups1, int* groups2, int gn1, int gn2, int na, float* dist) {
    float mindist;
    int a,b,g1,g2,n1,n2,g1atm,g2atm;

    //for (a=0; a<4; a++){
    //    const float* coo2 = &coords[a*dim];
    //    printf("%f %f %f\n", coo2[0], coo2[1], coo2[2]);
    //}

    // Iterate over the two group sets
    for (g1=0; g1 < gn1; g1++){
        for (g2=0; g2 < gn2; g2++){
            mindist = -1;
            // Iterate over atoms in the two groups
            for (a = 0; a < na; a++) {
                g1atm = groups1[g1 * na + a];
                if (g1atm == -1) break;
                const float* coo1 = &coords[g1atm*dim];

                for (b=0; b < na; b++) {
                    g2atm = groups2[g2 * na + b];
                    if (g2atm == -1) break;
                    const float* coo2 = &coords[g2atm*dim];

                    float d[3];
                    d[0] = coo1[0]-coo2[0];
                    d[1] = coo1[1]-coo2[1];
                    d[2] = coo1[2]-coo2[2];
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