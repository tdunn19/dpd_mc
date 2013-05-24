/* function.c */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "dpd.h"

Vector vdist(Vector v, Vector u) {
    Vector dr;

    dr.x = v.x - u.x;
    dr.y = v.y - u.y;
    dr.z = v.z - u.z;

    return dr;
}

double vmag(Vector v) {
        return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

int mod(int a, int b) {
    int c = a % b;

    if (c < 0) {
        c += b;
    }

    return c;
}
