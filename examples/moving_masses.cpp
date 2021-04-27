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


#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <LCQProblem.hpp>
#include <casadi/casadi.hpp>

using casadi::SX;
using casadi::DM;
using casadi::Function;
using std::vector;
using std::endl;
using std::string;
using lcqpOASES::LCQProblem;

// Settings and flags
int nMasses = 2;                                // Number of masses
int T = 2;                                      // Time span [0, T]
int N = 75;                                     // Number of discretization nodes
double h = (double)T/N;                         // Step size
bool useIPOPT = false;                          // Use qpOASES LCQP or IPOPT NLP
double u_max = 1000;                            // Bound on control
double controlPenalty = 1;                      // Penalty factor on control cost
double constraintTol = 1e-7;                    // Terminal constraint tolerance

// State space
int nu = 1;                                     // Number of control variables per node
int nx = 2*nMasses + 1;                         // Number of state variables per node
int nz = 2*nMasses;                             // Number of algebraic variables per node
int NV = 1 + nx + N*(nu + nx + nz);           // Number of optimization variables
int NComp = 4*N*nMasses;                      // Number of complementarity constraints
int NConstr = N*nMasses + nx*N + nMasses;     // Number of constraints
double* C = new double[NV*NV];                  // S1'*S2

// Initial state of dynamic system
double* x0 = new double[nx];
double* u0 = new double[nu];
double* z0 = new double[nz];
double* lbx = new double[nx];
double* ubx = new double[nx];
double* lbu = new double[nu];
double* ubu = new double[nu];
double* lbz = new double[nz];
double* ubz = new double[nz];

// Complementarity matrices (Store a la data[row*ncol + col] = data)
double* S1 = new double[NComp*NV];
double* S2 = new double[NComp*NV];
int compCounter = 0;

// State variables
SX p = SX::sym("p", nMasses);
SX v = SX::sym("p", nMasses);
SX t = SX::sym("t", 1);
SX x = vertcat(p, v, t);

// Control variables
SX u = SX::sym("u", nu);

// Algebraic variables
SX y = SX::sym("y", nMasses);
SX lambda0 = SX::sym("lambda0", nMasses);
SX z = vertcat(y, lambda0);

// Parameters
SX one = SX::sym("one", 1);

// Trajectory variables, bounds and indices
vector<SX> w_k;
vector<double> lbw;
vector<double> ubw;
vector<int> ind_all;
vector<int> ind_x;
vector<int> ind_u;
vector<int> ind_z;

// ODE constraints
vector <SX> g_k;
vector<double> lbg;
vector<double> ubg;

// Initial guess vector
vector<double> w0;

// Initialize cost function
SX J = 0;

// dx/dt = f(x,u)
SX f(const SX& x, const SX& u, const SX& z) {

    SX dot = SX::sym("f", nx);

    // p_dot = v
    for (int i = 0; i < nMasses; i++) {
        dot(i) = x(nMasses + i);
    }

    // v_dot = ...
    dot(nMasses) = -2*x(0) + x(1) - x(nMasses) - 0.3*(2*z(0) - 1);
    for (int i = 1; i < nMasses - 1; i++) {
        dot(nMasses + i) = x(i - 1) - 2*x(i) + x(i + 1) - x(nMasses + i) - 0.3*(2*z(i) - 1);
    }
    dot(2*nMasses - 1) = x(nMasses - 2) - 2*x(nMasses - 1) - x(2*nMasses - 1) - 0.3*(2*z(nMasses - 1) - 1) + u;

    // t_dot = 1
    dot(2*nMasses) = 1;

    return dot;
}

void FillC(const double* const S1, const double* const S2) {
    for (int i = 0; i < NV; i++) {
        for (int j = 0; j <= i; j++) {
            double res = 0;
            for (int k = 0; k < NComp; k++) {
                res += S1[k*NV + i]*S2[k*NV + j];
            }

            C[i*NV + j] = res;
            C[j*NV + i] = res;
        }
    }
}

void CreateQPfromCasADi(    SX J, SX g, SX w,
                            const double* const lbw, const double* const ubw,
                            const double* const lbg, const double* const ubg,
                            double* Q, double* d, double* A, double* lb, double* ub, double* lbA, double* ubA) {

    // Build zero argument
    vector<double> zero_vec;
    for (int i = 0; i < NV; i++) {
        zero_vec.push_back(0);
    }
    vector<DM> arg_zero = {reshape(DM(zero_vec), NV, 1)};

    // Create QP
    casadi::Dict dict = casadi::Dict();

    // 1) Objective
    SX gradJ = gradient(J, w, dict);
    SX hessJ = jacobian(gradJ, w, dict);
    Function GradJ("GradJ", {w}, {gradJ});
    Function HessJ("HessianJ", {w}, {hessJ});

    // Evaluate constraint linarization at 0
    vector<double> HessJ_eval = HessJ(arg_zero).at(0).get_elements();
    vector<double> GradJ_eval = GradJ(arg_zero).at(0).get_elements();

    // 2) Constraints
    SX jacG = jacobian(g, w, dict);
    Function G("G", {w}, {g});
    Function JacG("JacG", {w}, {jacG});

    // Evaluate constraint linarization at 0
    vector<double> Gk_zero_eval = G(arg_zero).at(0).get_elements();
    vector<double> JacGk_eval = JacG(arg_zero).at(0).get_elements();

    // Store objective
    for (int i = 0; i < NV; i++) {
        d[i] = GradJ_eval[i];
        for (int j = 0; j < NV; j++){
            Q[i*NV + j] = HessJ_eval[i*NV + j];
        }
    }

    // Store constraints
    double* Constr_const = new double[NConstr];
    for (int i = 0; i < NConstr; i++) {
        // Adjust constraint bounds
        Constr_const[i] = Gk_zero_eval[i];
        lbA[i] = lbg[i] - Constr_const[i];
        ubA[i] = ubg[i] - Constr_const[i];

        // Store jacobian in A (Careful: Jacobian is returned transposed for some reason?)
        for (int j = 0; j < NV; j++)
            A[i*NV + j] = JacGk_eval[j*NConstr + i];
    }

    // Fill box constraints
    for (int i = 0; i < NV; i++) {
        lb[i] = lbw[i];
        ub[i] = ubw[i];
    }
}

int main() {
    std::cout << "Building moving masses LCQP" << endl;

    // Bounds and initial guess for states
    for (int i = 0; i < nx; i++) {
        x0[i] = i % 2;
        lbx[i] = -casadi::inf;
        ubx[i] = casadi::inf;
    }
    x0[nx - 1] = 0;

    // Bounds and initial guess for controls
    for (int i = 0; i < nu; i++) {
        u0[i] = 0;
        lbu[i] = -u_max;
        ubu[i] = u_max;
    }

    // Bounds and initial guess for algebraic variables
    for (int j = 0; j < nMasses; j++) {
        // Bounds on y
        lbz[j] = 0;
        ubz[j] = 1;

        // Bounds on lambda0
        lbz[nMasses + j] = 0;
        ubz[nMasses + j] = casadi::inf;

        // Initial guess
        if (x0[nMasses + j] < 0) {
            z0[nMasses + j] = -x0[nMasses + j];
            z0[j] = 1;
        } else {
            z0[nMasses + j] = 0;
            z0[j] = 0;
        }
    }

    // Begin with parameters (one)
    w_k.push_back(one);
    w0.push_back(1);
    lbw.push_back(1);
    ubw.push_back(1);
    ind_all.push_back(1);

    // Lift initial state
    SX xk = SX::sym("x_0", nx);
    w_k.push_back(xk);
    for (int i = 0; i < nx; i++) {
        // Add to initial guess
        w0.push_back(x0[i]);

        // Add equality box constraints
        lbw.push_back(x0[i]);
        ubw.push_back(x0[i]);

        // Add to index vectors
        ind_x.push_back(1 + i);
        ind_all.push_back(1 + i);
    }

    // Node loop
    for (int k = 0; k < N; k++) {
        // New index name suffix
        string nameInd = std::to_string(k);

        SX uk = SX::sym("u_" + nameInd, nu);
        if (nu > 0) {
            w_k.push_back(uk);
            int new_u_start_ind = ind_all.back() + 1;
            for (int i = 0; i < nu; i++) {
                // Add to initial guess
                w0.push_back(u0[i]);

                // Add box constraints
                lbw.push_back(lbu[i]);
                ubw.push_back(ubu[i]);

                // Add to index vectors
                ind_u.push_back(new_u_start_ind + i);
                ind_all.push_back(new_u_start_ind + i);

                // Control penalization
                J = J + controlPenalty*uk(i)*uk(i);
            }
        }

        // New state variables
        SX xkj = SX::sym("x_" + nameInd, nx);
        w_k.push_back(xkj);
        int new_x_start_ind = ind_all.back() + 1;
        for (int i = 0; i < nx; i++) {
            // Add to initial guess
            w0.push_back(x0[i]);

            // Add box constraints
            lbw.push_back(lbx[i]);
            ubw.push_back(ubx[i]);

            // Add to index vectors
            ind_x.push_back(new_x_start_ind + i);
            ind_all.push_back(new_x_start_ind + i);
        }

        // New algebraic variables
        SX zkj = SX::sym("z_" + nameInd, nz);
        w_k.push_back(zkj);
        int new_z_start_ind = ind_all.back() + 1;
        for (int i = 0; i < nz; i++) {
            // Add to initial guess
            w0.push_back(z0[i]);

            // Add box constraints
            lbw.push_back(lbz[i]);
            ubw.push_back(ubz[i]);

            // Add to index vectors
            ind_z.push_back(new_z_start_ind + i);
            ind_all.push_back(new_z_start_ind + i);
        }

        // New complementarity constraints
        for (int i = 0; i < nMasses; i++) {
            // New indices
            int yInd = new_z_start_ind + i;
            int lam0Ind = new_z_start_ind + nMasses + i;
            int vInd = new_x_start_ind + nMasses + i;

            // y*lambda0
            S1[compCounter*NV + yInd] = 1;
            S2[compCounter*NV + lam0Ind] = 1;
            compCounter++;

            // lambda0*y (for symmetry)
            S2[compCounter*NV + yInd] = 1;
            S1[compCounter*NV + lam0Ind] = 1;
            compCounter++;

            // (1-y)*(lambda0 + v)
            S1[compCounter*NV + 0] = 1;
            S1[compCounter*NV + yInd] = -1;
            S2[compCounter*NV + lam0Ind] = 1;
            S2[compCounter*NV + vInd] = 1;
            compCounter++;

            // (lambda0 + v)*(1-y) (for symmetry)
            S2[compCounter*NV + 0] = 1;
            S2[compCounter*NV + yInd] = -1;
            S1[compCounter*NV + lam0Ind] = 1;
            S1[compCounter*NV + vInd] = 1;
            compCounter++;
        }

        // Impose posivity on lambda1 = v + lambda0
        SX gj = SX::sym("g_" + nameInd, nMasses);
        for (int i = 0; i < nMasses; i++) {
            gj(i) = xkj(nMasses + i) + zkj(nMasses + i);
            lbg.push_back(0);
            ubg.push_back(casadi::inf);
        }

        g_k.push_back(gj);

        SX xk_end = SX::sym("xk_end", nx);
        for (int i = 0; i < nx; i++)
            xk_end(i) = -1*xk(i);

        // Loop over collocation points
        SX xp = SX::sym("xp", nx);
        for (int i = 0; i < nx; i++)
            xp(i) = -2*xk(i) + 2*xkj(i);

        SX fj = f(xkj, uk, zkj);
        gj = SX::sym("gj", nx);

        for (int i = 0; i < nx; i++) {
            gj(i) = h*fj(i) - xp(i);
            lbg.push_back(0);
            ubg.push_back(0);
        }

        g_k.push_back(gj);

        for (int i = 0; i < nx; i++) {
            xk_end(i) += 2*xkj(i);
            xk(i) = xk_end(i);
        }
    }

    // Terminal constraint
    SX gj = SX::sym("gj", nMasses);

    for (int i = 0; i < nMasses; i++) {
        gj(i) = xk(nMasses + i);
        lbg.push_back(-constraintTol);
        ubg.push_back(constraintTol);
    }

    g_k.push_back(gj);

    std::cout << "Building moving masses LCQP done ..." << endl;

    FillC(S1, S2);

    // Create optimization variables and constraints
    SX w = vertcat(w_k);
    SX g = vertcat(g_k);

    vector<double> xopt;


    std::cout << "Calling LCQP solver ..." << endl;

    double* b = new double[NComp];
    for (int i = 0; i < NComp; i++)
        b[i] = 1;

    // Allocate LCQP arrays to be filled
    double* Q = new double[NV*NV];
    double* d = new double[NV];
    double* A = new double[NConstr*NV];
    double* lbA = new double[NConstr];
    double* ubA = new double[NConstr];
    double* lb = new double[NV];
    double* ub = new double[NV];

    // Interface to generate LCQP from CasADi functions and variables
    CreateQPfromCasADi(J, g, w, &lbw[0], &ubw[0], &lbg[0], &ubg[0], Q, d, A, lb, ub, lbA, ubA);

    // Initialize LCQP solver
    LCQProblem lcqp( NV, NConstr, NComp);

    // lcqpOASES options and settings
    lcqpOASES::Options options;
    options.initialComplementarityPenalty = 1.0;
    options.complementarityPenaltyUpdate = 10.0;
    lcqp.setOptions( options );

    // Set maximum number of working set changes
    int nWSR = 10000;

    // Set maximum cpu time to unlimited
    double* cputime = 0;

    // Solve LCQP
    lcqp.loadLCQP(Q, d, lb, ub, S1, S2, A, lbA, ubA, &w0[0]);
    lcqp.runSolver();

    double* wk_opt = new double[NV];
    lcqp.getPrimalSolution( wk_opt );
}