//
// Created by userpc on 14.11.2018.
//

#pragma once

#include "functions_bvp.h"
#include "BVP.h"

void test_equation(double x0, double xn, double h, func p, func q, func f, func analitic, boundary_set &set);
