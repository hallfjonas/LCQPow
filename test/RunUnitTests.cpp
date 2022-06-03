/*
 *	This file is part of LCQPow.
 *
 *	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2022 by Jonas Hall et al.
 *
 *	LCQPow is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPow is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPow; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "Utilities.hpp"
#include "Options.hpp"
#include "LCQProblem.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

// Testing standard matrix multiplications
TEST(UtilitiesTest, MatrixMultiplicationTest) {
    // A = [1 0 2; 3 1 1]
    // B = [2 0 0 2; 1 0 0 1; 0 -1 -1 0]
    // C = A*B = [2 -2 -2 2; 7 -1 -1 7]
    int m = 2;
    int n = 3;
    int p = 4;

    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1 };
    double* B = new double[n*p] { 2, 0, 0, 2, 1, 0, 0, 1, 0, -1, -1, 0 };
    double* C = new double[m*p];

    LCQPow::Utilities::MatrixMultiplication(A, B, C, m, n, p);

    ASSERT_EQ(C[0], 2);
    ASSERT_EQ(C[1], -2);
    ASSERT_EQ(C[2], -2);
    ASSERT_EQ(C[3], 2);
    ASSERT_EQ(C[4], 7);
    ASSERT_EQ(C[5], -1);
    ASSERT_EQ(C[6], -1);
    ASSERT_EQ(C[7], 7);
}

// Testing transposed matrix multiplications
TEST(UtilitiesTest, TransposedMatrixMultiplicationTest) {
    // A = [1 0 2; 3 1 1]
    // B = [98 -10]
    // C = A*B = [68 -10 186]
    int m = 2;
    int n = 3;
    int p = 1;

    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1 };
    double* B = new double[n*m] { 98, -10 };
    double* C = new double[n*p];

    LCQPow::Utilities::TransponsedMatrixMultiplication(A, B, C, m, n, p);
    ASSERT_EQ(C[0], 68);
    ASSERT_EQ(C[1], -10);
    ASSERT_EQ(C[2], 186);
}

// Testing the matrix symmetrization product
TEST(UtilitiesTest, MatrixSymmetrization) {
    // A = [1 0 2; 3 1 1]
    // B = [2 0 1; 0 0 -1]
    // C = A'*B + B'*A = [4 1 -2; 1 0 0; -2 0 -4]
    int m = 2;
    int n = 3;

    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1 };
    double* B = new double[n*m] { 2, 0, 1, 0, 0, -1 };
    double* C = new double[n*n];
    LCQPow::Utilities::MatrixSymmetrizationProduct(A, B, C, m, n);

    ASSERT_EQ(C[0], 4);
    ASSERT_EQ(C[1], 0);
    ASSERT_EQ(C[2], 2);
    ASSERT_EQ(C[3], 0);
    ASSERT_EQ(C[4], 0);
    ASSERT_EQ(C[5], -1);
    ASSERT_EQ(C[6], 2);
    ASSERT_EQ(C[7], -1);
    ASSERT_EQ(C[8], 2);
}

// Testing standard and symmetrization matrix multiplications
TEST(UtilitiesTest, AffineTransformation) {
    // alpha = 2;
    // A = [1 0 2; 3 1 1]
    // b = [2 0 1]
    // c = [-3; -3]
    // d = [2; 8] = alpha*A*b + c

    int m = 2;
    int n = 3;

    double alpha = 2;
    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1 };
    double* b = new double[n*m] { 2, 0, 1 };
    double* c = new double[m] { -3, -3 };
    double* d = new double[m];

    LCQPow::Utilities::AffineLinearTransformation(alpha, A, b, c, d, m, n);
    ASSERT_EQ(d[0], 5);
    ASSERT_EQ(d[1], 11);
}

// Testing matrix add
TEST(UtilitiesTest, MatrixAdd) {
    // alpha = -1;
    // A = [0 1; 3 1; 10 1]
    // beta = 0.5;
    // B = [2 0; 0 4; 2 2]
    // C = [1 -1; -3 1; -9 0]

    int m = 3;
    int n = 2;

    double alpha = -1;
    double* A = new double[m*n] { 0, 1, 3, 1, 10, 1 };

    double beta = 0.5;
    double* B = new double[m*n] { 2, 0, 0, 4, 2, 2 };
    double* C = new double[m*n];

    LCQPow::Utilities::WeightedMatrixAdd(alpha, A, beta, B, C, m, n);
    ASSERT_EQ(C[0], 1);
    ASSERT_EQ(C[1], -1);
    ASSERT_EQ(C[2], -3);
    ASSERT_EQ(C[3], 1);
    ASSERT_EQ(C[4], -9);
    ASSERT_EQ(C[5], 0);
}

// Testing vector add
TEST(UtilitiesTest, VectorAdd) {
    // alpha = 2;
    // a = [0; 1; 2; 3]
    // beta = -1;
    // b = [10; 2; 0; 3]
    // d = [-10; 0; 4; -3]

    int m = 4;

    double alpha = 2;
    double* a = new double[m] { 0, 1, 2, 3 };

    double beta = -1;
    double* b = new double[m] { 10, 2, 0, 3 };
    double* d = new double[m];

    LCQPow::Utilities::WeightedVectorAdd(alpha, a, beta, b, d, m);
    ASSERT_EQ(d[0], -10);
    ASSERT_EQ(d[1], 0);
    ASSERT_EQ(d[2], 4);
    ASSERT_EQ(d[3], 3);
}

// Testing Quadratic Form Product
TEST(UtilitiesTest, QuadraticFormProduct) {
    // p = [1; 2; 3]
    // Q = [0 1 0; 1 2 1; 0 1 0]
    // ret = [1 2 3]*[2; 8; 2] = 24

    int m = 3;
    double* p = new double[m] { 1, 2, 3 };
    double* Q = new double[m*m] { 0, 1, 0, 1, 2, 1, 0, 1, 0 };

    double ret = LCQPow::Utilities::QuadraticFormProduct(Q, p, m);
    ASSERT_EQ(ret, 24);
}

// Testing dot product
TEST(UtilitiesTest, DotProduct) {
    // a = [0; 1; 2; 3]
    // b = [10; 2; 0; 3]
    // ret = a'*b

    int m = 4;
    double* a = new double[m] { 0, 1, 2, 3 };
    double* b = new double[m] { 10, 2, 0, 3 };

    double ret = LCQPow::Utilities::DotProduct(a, b, m);
    ASSERT_EQ(ret, 11);
}

// Testing 1 norm
TEST(UtilitiesTest, MaxAbs) {
    int m = 4;

    // a = [0; 1; 2; 3]
    // ret = 3
    double* a = new double[m] { 0, 1, 2, 3 };
    double ret = LCQPow::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 3);

    // a = [0; -1; 2; 0]
    // ret = 2
    a = new double[m] { 0, -1, 2, 0 };
    ret = LCQPow::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 2);

    // a = [0; -4; 2; 0]
    // ret = 4
    a = new double[m] { 0, -4, 2, 0 };
    ret = LCQPow::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 4);
}

// Testing Options constructors, default settings, consistency
TEST(UtilitiesTest, Options) {
    LCQPow::Options opts;

    // Check changed values
    opts.setInitialPenaltyParameter( 100 );
    opts.setPenaltyUpdateFactor( 100 );
    ASSERT_EQ(opts.getInitialPenaltyParameter(), 100);
    ASSERT_EQ(opts.getPenaltyUpdateFactor(), 100);

    // Check copy constructor
    LCQPow::Options opts2(opts);
    ASSERT_EQ(opts2.getInitialPenaltyParameter(), 100);
    ASSERT_EQ(opts2.getPenaltyUpdateFactor(), 100);
}

// Testing csc to dns
TEST(UtilitiesTest, CSCtoDNS) {
    // First test matrix
    // | 2 1 0 |
    // | 0 2 0 |

    int m = 2;
    int n = 3;
    int Q_nnx = 3;
    double Q_data[3] = { 2.0, 1.0, 2.0 };
    int Q_i[3] = {0, 0, 1};
    int Q_p[4] = {0, 1, 3, 4};

    csc* Q = csc_matrix(m, n, Q_nnx, Q_data, Q_i, Q_p);

    double* Q_full = LCQPow::Utilities::csc_to_dns(Q);
    ASSERT_TRUE(Q_full != 0);

    ASSERT_DOUBLE_EQ(Q_full[0], 2);
    ASSERT_DOUBLE_EQ(Q_full[1], 1);
    ASSERT_DOUBLE_EQ(Q_full[2], 0);
    ASSERT_DOUBLE_EQ(Q_full[3], 0);
    ASSERT_DOUBLE_EQ(Q_full[4], 2);
    ASSERT_DOUBLE_EQ(Q_full[5], 0);
    delete[] Q_full;

    // Modify some values: Second test matrix
    // | 2 0 0  |
    // | 1 0 10 |
    Q_data[2] = 10;
    Q_p[1] = 2;
    Q_p[2] = 2;
    Q_i[1] = 1;
    Q = csc_matrix(m, n, Q_nnx, Q_data, Q_i, Q_p);

    Q_full = LCQPow::Utilities::csc_to_dns(Q);
    ASSERT_TRUE(Q_full != 0);

    ASSERT_EQ(Q_full[0], 2);
    ASSERT_EQ(Q_full[1], 0);
    ASSERT_EQ(Q_full[2], 0);
    ASSERT_EQ(Q_full[3], 1);
    ASSERT_EQ(Q_full[4], 0);
    ASSERT_EQ(Q_full[5], 10);
    delete[] Q_full;

    // Transpose: Third test matrix
    // | 2  1 |
    // | 0  0 |
    // | 10 0 |
    m = 3;
    n = 2;
    double T_data[3] = { 2.0, 10.0, 1.0 };
    int T_i[3] = {0, 2, 0};
    int T_p[3] = {0, 2, 3};
    int T_nnx = 3;
    csc* T = csc_matrix(m, n, T_nnx, T_data, T_i, T_p);

    double* T_full = LCQPow::Utilities::csc_to_dns(T);
    ASSERT_TRUE(T_full != 0);

    ASSERT_EQ(T_full[0], 2.0);
    ASSERT_EQ(T_full[1], 1.0);
    ASSERT_EQ(T_full[2], 0.0);
    ASSERT_EQ(T_full[3], 0.0);
    ASSERT_EQ(T_full[4], 10.0);
    ASSERT_EQ(T_full[5], 0.0);
    delete[] T_full;
}

// Testing matrix format change (csc to dns and vice versa)
TEST(UtilitiesTest, SparseDenseBackAndForth) {

    int numExp = 100;
    int m = 2;
    int n = 5;

    srand((unsigned int)time(NULL));

    for (int i = 0; i < numExp; i++) {
        double* Q = new double[m*n]();

        // Monitor number of nonzeros
        int nnx = 0;

        // Randomly fill values
        for (int j=0; j < m*n; j++)
        {
            // Get random integers
            int rd = std::rand();

            // Only fill even integers (i.e should be rougly 25% filled).
            if (rd % 4 == 0) {
                Q[j] = rd % 9;
                nnx++;
            }
        }

        // Convert to sparse
        csc* Q_sparse = LCQPow::Utilities::dns_to_csc(Q, m, n);

        double* Q_control = LCQPow::Utilities::csc_to_dns(Q_sparse);

        for (int j = 0; j < m*n; j++)
            ASSERT_DOUBLE_EQ(Q_control[j], Q[j]);

        // Clear memory
        delete[] Q;
        delete[] Q_control;
        c_free(Q_sparse);
    }
}

// Testing csc to triangular
TEST(UtilitesTest, CSCtoTriangular) {
    double M_data[4] = { 2.0, 3.0, 3.0, 2.0 };
    int M_i[4] = {0, 1, 0, 1};
    int M_p[3] = {0, 2, 4};

    csc* M = (csc*) malloc(sizeof(csc));

    M->m = 2;
    M->n = 2;
    M->i = M_i;
    M->p = M_p;
    M->x = M_data;
    M->nzmax = 4;
    M->nz = -1;

    csc* M_triag = LCQPow::Utilities::copyCSC(M, true);

    ASSERT_TRUE(M_triag != 0);

    ASSERT_DOUBLE_EQ(M_triag->p[0], 0);
    ASSERT_DOUBLE_EQ(M_triag->p[1], 1);
    ASSERT_DOUBLE_EQ(M_triag->p[2], 3);
    ASSERT_DOUBLE_EQ(M_triag->i[0], 0);
    ASSERT_DOUBLE_EQ(M_triag->i[1], 0);
    ASSERT_DOUBLE_EQ(M_triag->i[2], 1);
    ASSERT_DOUBLE_EQ(M_triag->x[0], 2);
    ASSERT_DOUBLE_EQ(M_triag->x[1], 3);
    ASSERT_DOUBLE_EQ(M_triag->x[2], 2);
    ASSERT_DOUBLE_EQ(M_triag->m, 2);
    ASSERT_DOUBLE_EQ(M_triag->n, 2);
    ASSERT_DOUBLE_EQ(M_triag->nz, -1);
    ASSERT_DOUBLE_EQ(M_triag->nzmax, 3);
}

// Testing solver data storage (sparse to dense and vice versa) 
TEST(LoadDataTest, DenseToSparse) {

    double Q[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double L[1*2] = {1.0, 0.0};
    double R[1*2] = {0.0, 1.0};
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    LCQPow::LCQProblem lcqp( nV, nC, nComp );

	LCQPow::Options options;
    options.setPrintLevel(LCQPow::PrintLevel::NONE);
    lcqp.setOptions( options );
    
    // Load dense data
    LCQPow::ReturnValue retVal = lcqp.loadLCQP( Q, g, L, R );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);
    
    // Switch to dense (should not do anything)
    retVal = lcqp.switchToDenseMode( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

    // Switch to sparse
    retVal = lcqp.switchToSparseMode( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

    // Switch to sparse (should not do anything)
    retVal = lcqp.switchToSparseMode( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

    // Solve in sparse mode
    options.setQPSolver(LCQPow::QPOASES_SPARSE);
    lcqp.setOptions( options );
    retVal = lcqp.runSolver( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

    // Switch back to dense
    retVal = lcqp.switchToDenseMode( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

    // Solve in dense mode
    options.setQPSolver(LCQPow::QPOASES_DENSE);
    lcqp.setOptions( options );
    retVal = lcqp.runSolver( );
    ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);   
}

// Testing LCQPow solver set up
TEST(SolverTest, RunWarmUp) {
    double Q[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double L[1*2] = {1.0, 0.0};
    double R[1*2] = {0.0, 1.0};
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    LCQPow::LCQProblem lcqp( nV, nC, nComp );

	LCQPow::Options options;
    options.setPrintLevel(LCQPow::PrintLevel::NONE);
    lcqp.setOptions( options );

    int numExp = 100;

    // Allocate solution vectors
    double* xOpt = new double[2];
    double* yOpt = new double[nV + nC + 2*nComp];

    for (int i = 0; i < numExp; i++) {

        LCQPow::ReturnValue retVal = lcqp.loadLCQP( Q, g, L, R );
        ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

        retVal = lcqp.runSolver( );
        ASSERT_EQ(retVal, LCQPow::SUCCESSFUL_RETURN);

        lcqp.getPrimalSolution( xOpt );
        lcqp.getDualSolution( yOpt );

        bool sStat1Found = (std::abs(xOpt[0] - 1) <= options.getStationarityTolerance()) && (std::abs(xOpt[1]) <= options.getStationarityTolerance());
        bool sStat2Found = (std::abs(xOpt[1] - 1) <= options.getStationarityTolerance()) && (std::abs(xOpt[0]) <= options.getStationarityTolerance());

        ASSERT_TRUE( sStat1Found || sStat2Found );

        bool stat1 = std::abs(2*xOpt[0] - 2 - yOpt[0] - yOpt[2]) <= options.getStationarityTolerance();
        bool stat2 = std::abs(2*xOpt[1] - 2 - yOpt[1] - yOpt[3]) <= options.getStationarityTolerance();

        ASSERT_TRUE( stat1 );
        ASSERT_TRUE( stat2 );
    }

    // Clear solution vectors
    delete[] xOpt; delete[] yOpt;
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}