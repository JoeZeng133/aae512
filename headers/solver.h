#pragma once

#include<iostream>
#include<vector>
#include<cmath>
#include<utils.h>

class solver {
private:
    int N, M;
    double d1, d2;
    double a, dt;
    double tol = 1e-7;
    int maxItr = 20;
    
    Mat xc, yc;
    Mat left, right, top, bottom, center, corner;
    Mat Told, tmp;
    Mat T;
    Mat g11, g12, g22;
    Mat L1, L2;
    
    double *dd, *bb;
    
    std::fstream fo;
public:
    solver() = default;
	void init(const Mat& xc, const Mat& yc, const double T0, const double T1, const double T2, const double _a, const double _dt);
    void snapshot(std::fstream& ft);
	double advance();
    void ex_advance();
};

void solver::init(const Mat& _xc, const Mat& _yc, const double T0, const double T1, const double T2, const double _a, const double _dt) {
    a = _a;
    dt = _dt;
    _xc.getSize(N, M);
    
    d1 = 1. / (N - 1);
    d2 = 1. / (M - 1);
    
    //padded with one periodic layers
    xc.zeros(N + 1, M);
    yc.zeros(N + 1, M);
    memcpy(&xc(1, 0), &_xc(0, 0), sizeof(double) * _xc.getSize());
    memcpy(&yc(1, 0), &_yc(0, 0), sizeof(double) * _yc.getSize());
    memcpy(&xc(0, 0), &xc(N - 1, 0), sizeof(double) * M);
    memcpy(&yc(0, 0), &yc(N - 1, 0), sizeof(double) * M);
    
    //dirichlet boundary assignment
    T.fill(N + 1, M, T0);
    for (int i = 0; i <= N; ++i) {
        T(i, 0) = T1;
        T(i, M - 1) = T2;
    }
    memcpy(&T(0, 0), &T(N - 1, 0), sizeof(double) * M);
    memcpy(&T(N, 0), &T(1, 0), sizeof(double) * M);
    
    Told.zeros(N + 1, M);
//    memcpy(Told.data(), T.data(), sizeof(double) * Told.getSize());

    // get metrics
    double g_11, g_12, g_22, g, J;
    double J11, J12, J21, J22;
    double sJ11, sJ12, sJ21, sJ22;
    double m1, m2;
    double gamma1_11, gamma1_12, gamma1_22;
    double gamma2_11, gamma2_12, gamma2_22;
    
    g11.zeros(N + 1, M);
    g22.zeros(N + 1, M);
    g12.zeros(N + 1, M);
    L1.zeros(N + 1, M);
    L2.zeros(N + 1, M);
    
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            //jocabian matrix
            J11 = FIR_1(xc, i, j) / (2 * d1);
            J12 = FIR_2(xc, i, j) / (2 * d2);
            J21 = FIR_1(yc, i, j) / (2 * d1);
            J22 = FIR_2(yc, i, j) / (2 * d2);
            
            sJ11 = SEC_1(xc, i, j) / SQR(d1);
            sJ12 = SEC_2(xc, i, j) / SQR(d2);
            sJ21 = SEC_1(yc, i, j) / SQR(d1);
            sJ22 = SEC_2(yc, i, j) / SQR(d2);
            
            m1 = SEC_MIX(xc, i, j) / (4 * d1 * d2);
            m2 = SEC_MIX(yc, i, j) / (4 * d1 * d2);
            
            //covariant
            g_11 = SQR(J11) + SQR(J21);
            g_12 = J11 * J12 + J21 * J22;
            g_22 = SQR(J12) + SQR(J22);
            g = g_11 * g_22 - SQR(g_12);
            J = J11 * J22 - J12 * J21;
            
            //contracovariant
            g11(i, j) = g_22 / g;
            g12(i, j) = -g_12 / g;
            g22(i, j) = g_11 / g;
            
            //Christofel symbols
            gamma1_11 = (J22 * sJ11 - J12 * sJ21) / J;
            gamma1_22 = (J22 * sJ12 - J12 * sJ22) / J;
            gamma1_12 = (J22 * m1 - J12 * m2) / J;
            
            gamma2_11 = (J11 * sJ21 - J21 * sJ11) / J;
            gamma2_22 = (J11 * sJ22 - J21 * sJ12) / J;
            gamma2_12 = (J11 * m2 - J21 * m1) / J;
            
            //laplacian
            L1(i, j) = (2 * g_12 * gamma1_12 - g_11 * gamma1_22 - g_22 * gamma1_11) / g;
            L2(i, j) = (2 * g_12 * gamma2_12 - g_11 * gamma2_22 - g_22 * gamma2_11) / g;
        }
    
    //get stencil constants
    double arg1 = 0.5 * a * dt;
    
    left.zeros(N + 1, M);
    right.zeros(N + 1, M);
    top.zeros(N + 1, M);
    bottom.zeros(N + 1, M);
    center.zeros(N + 1, M);
    corner.zeros(N + 1, M);
    
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            right(i, j) = (g11(i, j) / SQR(d1) + L1(i, j) / (2 * d1)) * arg1;
            left(i, j)  = (g11(i, j) / SQR(d1) - L1(i, j) / (2 * d1)) * arg1;
            top(i, j)   = (g22(i, j) / SQR(d2) + L2(i, j) / (2 * d2)) * arg1;
            bottom(i, j)= (g22(i, j) / SQR(d2) - L2(i, j) / (2 * d2)) * arg1;
            corner(i, j)= (g12(i, j) / (2 * d1 * d2)) * arg1;
            center(i, j)= -(right(i, j) + left(i, j) + top(i, j) + bottom(i, j));
        }
    
    //tridiagnal arrays init
    dd = new double[M];
    bb = new double[M];
}

void solver::ex_advance() {
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            Told(i, j) = right(i, j) * T(i + 1, j) + left(i, j) * T(i - 1, j) + top(i, j) * T(i, j + 1) + bottom(i, j) * T(i, j - 1) + (center(i, j)) * T(i, j) + corner(i, j) * SEC_MIX(T, i, j);
        }
    
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            T(i, j) = 2 *  Told(i, j) + T(i, j);
        }
    
    memcpy(&T(0, 0), &T(N - 1, 0), sizeof(double) * M);
    memcpy(&T(N, 0), &T(1, 0), sizeof(double) * M);
}

double solver::advance() {
    double error = 0;
    
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            Told(i, j) = right(i, j) * T(i + 1, j) + left(i, j) * T(i - 1, j) + top(i, j) * T(i, j + 1) + bottom(i, j) * T(i, j - 1) + center(i, j) * T(i, j) + corner(i, j) * SEC_MIX(T, i, j);
            Told(i, j) += T(i, j);
        }
    
    for(int itr = 0; itr < maxItr; ++itr) {
        tmp = T;

        for (int i = 1; i < N; ++i) {
            for (int j = 1; j < M - 1; ++j) {
                dd[j] = center(i, j) - 1;
                bb[j] = -right(i, j) * T(i + 1, j) - left(i, j) * T(i - 1, j) - corner(i, j) * SEC_MIX(T, i, j) - Told(i, j);
            }

            //dirichlet boundary
            bb[1] -= bottom(i, 1) * T(i, 0);
            bb[M - 2] -= top(i, M - 2) * T(i, M - 1);
            
            tridiag(&bottom(i, 1), dd + 1, &top(i, 1), bb + 1, &T(i, 1), M - 2);
        }
        
        //periodic boundary
        memcpy(&T(0, 0), &T(N - 1, 0), sizeof(double) * M);
        memcpy(&T(N, 0), &T(1, 0), sizeof(double) * M);

        error = calcError(tmp, T);
        std::cout << "\tError at Iteration " << itr << " = " << error << "\n";
        if (error < tol) {
            std::cout << "\tSolution Converged\n";
            break;
        }
    }
    
    return error;
}

void solver::snapshot(std::fstream& ft) {
    for(int j = 0; j < M; ++j)
        for (int i = 1; i <= N; ++i)
            ft.write((char*)&T(i, j), sizeof(double));
}
