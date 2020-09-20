#include "BVP.h"

// https://old.math.tsu.ru/EEResources/cm/

void BVP::show_bvp(double &x0, double &xn, double &h, Vector &sol1, Vector &sol2, func analitic) {
    cout
            << setw(5) << "x"
            << setw(20) << "first approx"
            << setw(20) << "second approx"
            << setw(20) << "ananlitic solut"
            << endl;

    int n = int((xn - x0) / h);
    double deviation_fa = 0;
    double deviation_sa = 0;
    for (int i = 0; i <= n; ++i) {
        cout
                << setw(5) << setprecision(4) << x0
                << setw(20) << setprecision(5) << sol1[i]
                << setw(20) << setprecision(5) << sol2[i]
                << setw(20) << setprecision(5) << analitic(x0)
                << endl;
        double err_first_app = fabs(sol1[i] - analitic(x0));
        double err_second_app = fabs(sol2[i] - analitic(x0));
        if (err_first_app > deviation_fa )
            deviation_fa = err_first_app;
        if (err_second_app > deviation_sa )
            deviation_sa = err_second_app;
        x0 += h;
    }
    if (analitic != nullptr) {
        cout << setw(5) << setprecision(4) << "Error first approximation: "  << deviation_fa << endl;
        cout << setw(5) << setprecision(4) << "Error second approximation: " << deviation_sa << endl;
    }
}

Vector BVP::bvp_coeff(double h, double x, func p, func q, Vector &s){
    double a = 1 / (h * h) * (1 - h / 2 * p(x));
    double b = 1 / (h * h) * (-2 + h * h * q(x));
    double c = 1 / (h * h) * (1 + h / 2 * p(x));
    return s = { a, b, c};
}

Vector BVP::boundary_firstapp(double &h, boundary_set &set){

    Vector s(6); s.reserve(6);

    double b0 = set.alpha[0] - set.alpha[1] / h;
    double c0 = set.alpha[1] / h;
    double d0 = set.bigA;
    double bn = set.betta[0] + set.betta[1] / h;
    double an = -set.betta[1] / h;
    double dn = set.bigB;

    s = { b0, c0, d0, bn, an, dn };
    return s;
}

Vector BVP::boundary_secondapp(double &x0, double &xn, double &h, func p, func q, func f, boundary_set &set){

    Vector s(6); s.reserve(6);
    Vector coeff(3); coeff.reserve(3);

    int n = int((xn - x0) / h) + 1;

    double x = h;
    BVP::bvp_coeff(h, x, p, q, coeff);
    double b0 = set.alpha[0] + set.alpha[1] * (coeff[0] - 3. * coeff[2]) / coeff[2] / 2. / h;
    double c0 = (coeff[1] + 4. * coeff[2]) * set.alpha[1] / coeff[2] / 2. / h;
    double d0 = set.bigA + set.alpha[1] * f(x) / coeff[2] / 2. / h;

    x = (n - 2) * h;
    BVP::bvp_coeff(h, x, p, q, coeff);
    double bn = set.betta[0] + set.betta[1] * (3. * coeff[0] - coeff[2]) / coeff[0] / 2. / h;
    double an = -set.betta[1] * (4. * coeff[0] + coeff[1]) / coeff[0] / 2. / h;
    double dn = set.bigB - set.betta[1] * f(x) / coeff[0] / 2. / h;

    s = { b0, c0, d0, bn, an, dn };
    return s;
}

Vector BVP::solve_bvp2 (double x0, double xn, double h, func p, func q, func f, boundary_set &set, Vector &coeff) {

    double x = x0 + h;
    unsigned int n = (unsigned int) ((xn - x0) / h) + 1;

    Vector a(n); a.reserve(n);
    Vector b(n); b.reserve(n);
    Vector c(n); c.reserve(n);
    Vector d(n); d.reserve(n);
    Vector P(n); P.reserve(n);
    Vector Q(n); Q.reserve(n);
    Vector y(n); y.reserve(n);

    for (int i = 0; i < n; i++) {
        a[i] = (1 / (h * h) * (1 - h / 2 * p(x)));
        b[i] = (1 / (h * h) * (-2 + h * h * q(x)));
        c[i] = (1 / (h * h) * (1 + h / 2 * (p(x))));
        d[i] = (f(x));
        x += h;
    }

    b[0] = coeff[0];
    c[0] = coeff[1];
    d[0] = coeff[2];

    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    b[n - 1] = coeff[3];
    a[n - 1] = coeff[4];
    d[n - 1] = coeff[5];

    for (int i = 1; i < n; ++i) {
        P[i] = -c[i] / (b[i] + a[i] * P[i - 1]);
        Q[i] = (d[i] - a[i] * Q[i - 1]) / (b[i] + a[i] * P[i - 1]);
    }

    for (int i = 0; i < n; i++) {
        y[n - 1 - i] = P[n - 1 - i] * y[n - i] + Q[n - 1 - i];
    }

    return y;
}

BVP::BVP() = default;

BVP::~BVP() = default;

