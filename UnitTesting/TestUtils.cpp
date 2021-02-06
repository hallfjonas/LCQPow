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

// Testing standard and symmetrization matrix multiplications
TEST(UtilitiesTest, MatrixMultiplicationTest) {
    // A = [1 0 2; 3 1 1]
    // B = [2 0; 1 0; 0 -1]
    // C = A*B = [2 -2; 7 -1]
    int m = 2;
    int n = 3;

    double* A = new double[m*n] { 1, 0, 2, 3, 1, 1 };
    double* B = new double[n*m] { 2, 0, 1, 0, 0, -1 };
    double* C = new double[m*m];

    lcqpOASES::Utilities::MatrixMultiplication(A, B, C, m, n, m);
    ASSERT_EQ(C[0], 2);
    ASSERT_EQ(C[1], -2);
    ASSERT_EQ(C[2], 7);
    ASSERT_EQ(C[3], -1);
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
    lcqpOASES::Utilities::MatrixSymmetrizationProduct(A, B, C, m, n);

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

    lcqpOASES::Utilities::AffineLinearTransformation(alpha, A, b, c, d, m, n);
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

    lcqpOASES::Utilities::WeightedMatrixAdd(alpha, A, beta, B, C, m, n);
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

    lcqpOASES::Utilities::WeightedVectorAdd(alpha, a, beta, b, d, m);
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
    
    double ret = lcqpOASES::Utilities::QuadraticFormProduct(Q, p, m);
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
    
    double ret = lcqpOASES::Utilities::DotProduct(a, b, m);
    ASSERT_EQ(ret, 11);
}

// Testing 1 norm
TEST(UtilitiesTest, MaxAbs) {
    int m = 4;

    // a = [0; 1; 2; 3]
    // ret = 3
    double* a = new double[m] { 0, 1, 2, 3 };    
    double ret = lcqpOASES::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 3);

    // a = [0; -1; 2; 0]
    // ret = 2
    a = new double[m] { 0, -1, 2, 0 };    
    ret = lcqpOASES::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 2);

    // a = [0; -4; 2; 0]
    // ret = 4
    a = new double[m] { 0, -4, 2, 0 };    
    ret = lcqpOASES::Utilities::MaxAbs(a, m);
    ASSERT_EQ(ret, 4);
}

// Testing read from file functionality
TEST(UtilitiesTest, ReadFromFile) {
    const char* fpath = "examples/example_data/C.txt";

    double* C = new double[4];
    lcqpOASES::Utilities::readFromFile(C, 4, fpath);

    ASSERT_EQ(C[0], 0);
    ASSERT_EQ(C[1], 1);
    ASSERT_EQ(C[2], 1);
    ASSERT_EQ(C[3], 0);
}

// Testing Options constructors, default settings, consistency
TEST(UtilitiesTest, Options) {
    lcqpOASES::Options opts;
    ASSERT_EQ(opts.ensureConsistency(), lcqpOASES::SUCCESSFUL_RETURN);
    
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

// Testing lcqpOASES solver set up
TEST(SolverTest, RunWarmUp) {
    double H[2*2] = { 2.0, 0.0, 0.0, 2.0 };
    double g[2] = { -2.0, -2.0 };
    double lb[2] = { 0, 0 };
    double ub[2] = { INFINITY, INFINITY };
    double S1[1*2] = {1.0, 0.0};
    double S2[1*2] = {0.0, 1.0};
    double x0[2] = { 1.0, 1.0 };
    int nV = 2;
    int nC = 0;
    int nComp = 1;

    lcqpOASES::LCQProblem lcqp( nV, nC, nComp );

	lcqpOASES::Options options;
    options.printLvl = lcqpOASES::printLevel::NONE;
    lcqp.setOptions( options );

	lcqpOASES::returnValue retVal = lcqp.solve( H, g, lb, ub, S1, S2, (double*)0, (double*)0, x0);

    ASSERT_EQ(retVal, lcqpOASES::SUCCESSFUL_RETURN);

    // Get solutions
    double* xOpt = new double[2];
	double* yOpt = new double[nV + nC + 2*nComp];

	lcqp.getPrimalSolution( xOpt );
	lcqp.getDualSolution( yOpt );

    bool sStat1Found = (std::abs(xOpt[0] - 1) <= options.stationarityTolerance) && (std::abs(xOpt[1]) <= options.stationarityTolerance);
    bool sStat2Found = (std::abs(xOpt[1] - 1) <= options.stationarityTolerance) && (std::abs(xOpt[0]) <= options.stationarityTolerance);

    ASSERT_TRUE(  sStat1Found || sStat2Found );

    bool stat1 = std::abs(2*xOpt[0] - 2 - yOpt[0] - yOpt[2]) <= options.stationarityTolerance;
    bool stat2 = std::abs(2*xOpt[1] - 2 - yOpt[1] - yOpt[3]) <= options.stationarityTolerance;
    
    ASSERT_TRUE( stat1 );
    ASSERT_TRUE( stat2 );
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}