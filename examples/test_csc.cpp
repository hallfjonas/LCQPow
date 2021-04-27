#include "osqp.h"

int main() {
    // Load problem data
    c_float P_x[3] = {4.0, 1.0, 2.0, };
    c_int P_nnz = 3;
    c_int P_i[3] = {0, 0, 1, };
    c_int P_p[3] = {0, 1, 3, };
    c_int n = 2;

    csc* P = csc_matrix(n, n, P_nnz, P_x, P_i, P_p);
    print_csc_matrix(P, "P");

    return 0;
}