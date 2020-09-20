#include "BVP.h"
#include "functions_bvp.h"
#include "test.h"

void test1() {
    cout << "y'' + y' = 1" << endl;
    test_equation(0, 1, 0.05, p1, q1, f1, analitic1, bvp1);
}

void test2() {
    cout << "y'' + 2.25 / (1.5x + 1) y' = 3 / sqrt(1.5x + 1)" << endl;
    test_equation(0, 1, 0.05, p2, q2, f2, analitic2, bvp2);
}


int main()
{
    test1();
    test2();

    return 0;
}
