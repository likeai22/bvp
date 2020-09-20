//
// Created by userpc on 14.11.2018.
//

#include "test.h"

void test_equation (double x0, double xn, double h, func p, func q, func f, func analitic, boundary_set &bc)
{
    Vector boundary_first_approximation = BVP::boundary_firstapp(h, bc);
    Vector boundary_second_approximation = BVP::boundary_secondapp(x0, xn, h, p, q, f, bc);

    Vector sol1 = BVP::solve_bvp2(x0, xn, h, p, q, f, bc, boundary_first_approximation);
    Vector sol2 = BVP::solve_bvp2(x0, xn, h, p, q, f, bc, boundary_second_approximation);

    BVP::show_bvp(x0, xn, h, sol1, sol2, analitic);

    sol1.clear();
    sol2.clear();
}