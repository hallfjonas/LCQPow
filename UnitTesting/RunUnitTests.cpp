/*
 *	This file is part of LCQPanther.
 *
 *	LCQPanther -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	LCQPanther is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	LCQPanther is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with LCQPanther; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "Utilities.hpp"
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

    LCQPanther::Utilities::MatrixMultiplication(A, B, C, m, n, p);

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

    LCQPanther::Utilities::TransponsedMatrixMultiplication(A, B, C, m, n, p);
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
    LCQPanther::Utilities::MatrixSymmetrizationProduct(A, B, C, m, n);

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

    LCQPanther::Utilities::AffineLinearTransformation(alpha, A, b, c, d, m, n);
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

    LCQPanther::Utilities::WeightedMatrixAdd(alpha, A, beta, B, C, m, n);
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

    LCQPanther::Utilities::WeightedVectorAdd(alpha, a, beta, b, d, m);
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

    double ret = LCQPanther::Utilities::QuadraticFormProduct(Q, p, m);
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

    double ret = LCQPanther::Utilities::DotProduct(a, b, m);
    ASSERT_EQ(ret, 11);
}

// Testing 1 norm
TEST(UtilitiesTest, MaxAbs) {
    int m = 4;

    // a = [0; 1; 2; 3]
    // ret = 3
    double* a = new double[m] { 0, 1, 2, 3 };
    double ret = LCQPanther::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 3);

    // a = [0; -1; 2; 0]
    // ret = 2
    a = new double[m] { 0, -1, 2, 0 };
    ret = LCQPanther::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 2);

    // a = [0; -4; 2; 0]
    // ret = 4
    a = new double[m] { 0, -4, 2, 0 };
    ret = LCQPanther::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 4);
}

// Testing read from file functionality
TEST(UtilitiesTest, ReadFromFile) {
    const char* fpath = "../examples/example_data/one_ivocp_example/lbA.txt";

    // Read first four values from lbA (0, 1, 0, 1)
    double* lbA = new double[4];
    LCQPanther::ReturnValue ret = LCQPanther::Utilities::readFromFile(lbA, 4, fpath);

    // Assert that reading file went ok
    ASSERT_EQ(ret, LCQPanther::ReturnValue::SUCCESSFUL_RETURN);

    // Assert that values have been read correctly
    ASSERT_EQ(lbA[0], 0);
    ASSERT_EQ(lbA[1], 1);
    ASSERT_EQ(lbA[2], 0);
    ASSERT_EQ(lbA[3], 1);
}

// Testing Options constructors, default settings, consistency
TEST(UtilitiesTest, Options) {
    LCQPanther::Options opts;

    // Check changed values
    opts.setInitialPenaltyParameter( 100 );
    opts.setPenaltyUpdateFactor( 100 );
    ASSERT_EQ(opts.getInitialPenaltyParameter(), 100);
    ASSERT_EQ(opts.getPenaltyUpdateFactor(), 100);

    // Check copy constructor
    LCQPanther::Options opts2(opts);
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
    int H_nnx = 3;
    double H_data[3] = { 2.0, 1.0, 2.0 };
    int H_i[3] = {0, 0, 1};
    int H_p[4] = {0, 1, 3, 4};

    csc* H = csc_matrix(m, n, H_nnx, H_data, H_i, H_p);

    double* H_full = LCQPanther::Utilities::csc_to_dns(H);
    ASSERT_TRUE(H_full != 0);

    ASSERT_DOUBLE_EQ(H_full[0], 2);
    ASSERT_DOUBLE_EQ(H_full[1], 1);
    ASSERT_DOUBLE_EQ(H_full[2], 0);
    ASSERT_DOUBLE_EQ(H_full[3], 0);
    ASSERT_DOUBLE_EQ(H_full[4], 2);
    ASSERT_DOUBLE_EQ(H_full[5], 0);
    delete[] H_full;

    // Modify some values: Second test matrix
    // | 2 0 0  |
    // | 1 0 10 |
    H_data[2] = 10;
    H_p[1] = 2;
    H_p[2] = 2;
    H_i[1] = 1;
    H = csc_matrix(m, n, H_nnx, H_data, H_i, H_p);

    H_full = LCQPanther::Utilities::csc_to_dns(H);
    ASSERT_TRUE(H_full != 0);

    ASSERT_EQ(H_full[0], 2);
    ASSERT_EQ(H_full[1], 0);
    ASSERT_EQ(H_full[2], 0);
    ASSERT_EQ(H_full[3], 1);
    ASSERT_EQ(H_full[4], 0);
    ASSERT_EQ(H_full[5], 10);
    delete[] H_full;

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

    double* T_full = LCQPanther::Utilities::csc_to_dns(T);
    ASSERT_TRUE(T_full != 0);

    ASSERT_EQ(T_full[0], 2.0);
    ASSERT_EQ(T_full[1], 1.0);
    ASSERT_EQ(T_full[2], 0.0);
    ASSERT_EQ(T_full[3], 0.0);
    ASSERT_EQ(T_full[4], 10.0);
    ASSERT_EQ(T_full[5], 0.0);
    delete[] T_full;
}

// Testing csc to dns and vice versa
TEST(UtilitiesTest, SparseDenseBackAndForth) {

    int numExp = 100;
    int m = 2;
    int n = 5;

    srand((unsigned int)time(NULL));

    for (int i = 0; i < numExp; i++) {
        double* H = new double[m*n]();

        // Monitor number of nonzeros
        int nnx = 0;

        // Randomly fill values
        for (int j=0; j < m*n; j++)
        {
            // Get random integers
            int rd = std::rand();

            // Only fill even integers (i.e should be rougly 25% filled).
            if (rd % 4 == 0) {
                H[j] = rd % 9;
                nnx++;
            }
        }

        // Convert to sparse
        csc* H_sparse = LCQPanther::Utilities::dns_to_csc(H, m, n);

        double* H_control = LCQPanther::Utilities::csc_to_dns(H_sparse);

        for (int j = 0; j < m*n; j++)
            ASSERT_DOUBLE_EQ(H_control[j], H[j]);

        // Clear memory
        delete[] H;
        delete[] H_control;
        c_free(H_sparse);
    }
}

// Testing LCQPanther solver set up
TEST(SolverTest, RunWarmUp) {
    double H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double S1[1*2] = {1.0, 0.0};
    double S2[1*2] = {0.0, 1.0};
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    LCQPanther::LCQProblem lcqp( nV, nC, nComp );

	LCQPanther::Options options;
    options.setPrintLevel(LCQPanther::PrintLevel::NONE);
    lcqp.setOptions( options );

    int numExp = 100;

    // Allocate solution vectors
    double* xOpt = new double[2];
    double* yOpt = new double[nV + nC + 2*nComp];

    for (int i = 0; i < numExp; i++) {

        LCQPanther::ReturnValue retVal = lcqp.loadLCQP( H, g, S1, S2 );
        ASSERT_EQ(retVal, LCQPanther::SUCCESSFUL_RETURN);

        retVal = lcqp.runSolver( );
        ASSERT_EQ(retVal, LCQPanther::SUCCESSFUL_RETURN);

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