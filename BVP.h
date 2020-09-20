#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

#include "functions_bvp.h"

using std::cout;
using std::endl;
using std::vector;
using std::setw;
using std::setprecision;

typedef double (*func)(double);
typedef vector<double> Vector;

struct boundary_set {
    Vector alpha;
    Vector betta;
    double bigA;
    double bigB;
};

class BVP {

public:
    BVP();
    ~BVP();

    static Vector solve_bvp2(double x0, double xn, double h, func p, func q, func f, boundary_set &set, Vector &coeff);
    static Vector bvp_coeff(double h, double step, func p, func q, Vector &s);
    static Vector boundary_firstapp(double &h, boundary_set &set);
    static Vector boundary_secondapp(double &x0, double &xn, double &h, func p, func q, func f, boundary_set &set);
    static void show_bvp(double &x0, double &xn, double &h, Vector &sol1, Vector &sol2, func analitic);
};

