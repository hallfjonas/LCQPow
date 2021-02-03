/*
 *	This file is part of lcqpOASES.
 *
 *	lcqpOASES -- A Solver for Quadratic Programs with Commplementarity Constraints.
 *	Copyright (C) 2020 - 2021 by Jonas Hall et al.
 *
 *	lcqpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	lcqpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with lcqpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


#include "Utilities.hpp"
#include "LCQProblem.hpp"
#include <gtest/gtest.h>

// Testing Options constructors, default settings, consistency
TEST(OptionsTest, DefaultSetting) {
    // Check constructor, default settings, consistency
    lcqpOASES::Options opts;
    ASSERT_EQ(opts.ensureConsistency(), lcqpOASES::returnValue::SUCCESSFUL_RETURN);
    
    // Check changed values
    opts.initialComplementarityPenalty = 100;
    opts.complementarityPenaltyUpdate = 100;
    ASSERT_EQ(opts.initialComplementarityPenalty, 100);
    ASSERT_EQ(opts.complementarityPenaltyUpdate, 100);
 
    // Check copy constructor
    lcqpOASES::Options opts2(opts);
    ASSERT_EQ(opts2.initialComplementarityPenalty, 100);
    ASSERT_EQ(opts2.complementarityPenaltyUpdate, 100);
}

// Testing standard and symmetrization matrix multiplications
TEST(UtilitiesTest, MatrixMultiplicationTest) {
    // 1) Standard matrix multiplication
    // A = [1 0 2; 3 1 1]
    // B = [2 0; 1 0; 0 -1]
    // C = A*B = [2 -2; 7 -1]
    int m = 2;
    int n = 3;

    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1};
    double* B = new double[n*m] { 2, 0, 1, 0, 0, -1};
    double* C = new double[m*m];

    lcqpOASES::Utilities::MatrixMultiplication(A, B, C, m, n, m);
    ASSERT_EQ(C[0], 2);
    ASSERT_EQ(C[1], -2);
    ASSERT_EQ(C[2], 7);
    ASSERT_EQ(C[3], -1);

    // 2) Matrix Symmetrization
    // A = [1 0 2; 3 1 1]
    // B = [2 0 1; 0 0 -1]
    // C = A'*B + B'*A = [4 1 -2; 1 0 0; -2 0 -4]
    double* D = new double[n*n];
    lcqpOASES::Utilities::MatrixSymmetrizationProduct(A, B, D, m, n);

    ASSERT_EQ(D[0], 4);
    ASSERT_EQ(D[1], 0);
    ASSERT_EQ(D[2], 2);
    ASSERT_EQ(D[3], 0);
    ASSERT_EQ(D[4], 0);
    ASSERT_EQ(D[5], -1);
    ASSERT_EQ(D[6], 2);
    ASSERT_EQ(D[7], -1);
    ASSERT_EQ(D[8], 2);
}

// Testing read from file functionality
TEST(ReadFromFileTest, ReadC) {
    const char* fpath = "../examples/example_data/C.txt";

    double* C = new double[4];
    lcqpOASES::Utilities::readFromFile(C, 4, fpath);

    ASSERT_EQ(C[0], 0);
    ASSERT_EQ(C[1], 1);
    ASSERT_EQ(C[2], 1);
    ASSERT_EQ(C[3], 0);
}

// Testing lcqpOASES solver set up
TEST(SolverSetUpTest, SetUpSolver) {
    int nV = 2;
    int nC = 1;
    int nComp = 1;

    double H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double lb[2] = { 0.0, 0.0 };
    double ub[2] = { 10000.0, 10000.0 };
    double S1[1*2] = {1.0, 0.0};
    double S2[1*2] = {0.0, 1.0};
    double A[1*2] = { 1.0, 0.0 };
    double lbA[1] = { -10.0 };
    double ubA[1] = { 100.0};

    // Maximal working set iterations (internally for qpOASES)
    int nwsr = 10000;

    // TODO: Linking qpOASES fails here :'(
    // lcqpOASES::LCQProblem lcqp(nV, nC, nComp);
    // TODO: Test solve
    // lcqp.solve(H, g, A, lb, ub, lbA, ubA, S1, S2, nwsr);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}