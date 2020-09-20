//
// Created by mvrd on 13.11.2018.
//

#include "functions_bvp.h"

// var1
boundary_set bvp1 = {
        .alpha = { 0, 1 },
        .betta = { 1, 0 },
        .bigA = 0,
        .bigB = 1
};
double p1(double x) { return 1; }
double q1(double x) { return 0;}
double f1(double x) { return 1; }
double analitic1(double x) { return -1 * exp(-1.) + exp(-1 * x) + x; }


//var2
boundary_set bvp2 = {
        .alpha = { 4.5, -1 },
        .betta = { 0, 1 },
        .bigA = 2.5 / 1.5,
        .bigB = sqrt(2.5)
};
double p2(double x) { return 2.25 / (1.5 * x + 1); }
double q2(double x) { return 0; }
double f2(double x) { return 3 / sqrt(1.5 * x + 1); }
double analitic2(double x) { return 0.666667 * sqrt(1.5 * x + 1) * (x + 0.666667) + 0.148148; }
