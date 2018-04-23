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
    double tol = 1e-5;
    
    int maxItr = 20;
    
    Mat xc, yc;
    Mat left, right, top, bottom, center, corner;
    
    Mat Told, tmp;
    Mat T;
    Mat g11, g12, g22;
    Mat L1, L2;
    
    double *dd, *bb;
public:
    solver() = default;
	void init(const Mat& _xc, const Mat& _yc);
	void advance();
};

void solver::init(const Mat& _xc, const Mat& _yc, const int T1, const int T2) {
    _xc.getSize(N, M);
    d1 = 1. / N;
    d2 = 1. / M;
    
    //padded with one periodic layers
    xc.zeros(N + 1, M);
    yc.zeros(N + 1, M);
    memcpy(&xc(1, 0), &_xc(0, 0), sizeof(double) * _xc.getSize());
    memcpy(&yc(1, 0), &_yc(0, 0), sizeof(double) * _yc.getSize());
    memcpy(&xc(0, 0), &xc(N - 1, 0), sizeof(double) * M);
    memcpy(&yc(0, 0), &yc(N - 1, 0), sizeof(double) * M);
    
    //dirichlet boundary assignment
    T.zeros(N + 1, M);
    for (int i = 0; i <= N; ++i) {
        T(i, 0) = T1;
        T(i, M - 1) = T2;
    }

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
        for (int j = 0; j < M - 1; ++j) {
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
            J = sqrt(g);
            
            //contracovariant
            g11(i, j) = g_22 / g;
            g12(i, j) = -g_12 / g;
            g22(i, j) = g_11 / g;
            
            //Christofel symbols
            gamma1_11 = (J22 * sJ11 - J12 * sJ21) / J;
            gamma1_22 = (J22 * sJ12 - J12 * sJ22) / J;
            gamma1_12 = (J22 * m1 - J12 * m2) / J;
            
            gamma2_11 = (J11 * sJ22 - J21 * sJ12) / J;
            gamma2_22 = (J11 * sJ22 - J21 * sJ12) / J;
            gamma2_12 = (J11 * m2 - J21 * m1) / J;
            
            //laplacian
            L1(i, j) = (2 * g_12 * gamma1_12 - g_11 * gamma1_22 - g_22 * gamma1_11) / g;
            L2(i, j) = (2 * g_12 * gamma2_12 - g_11 * gamma2_22 - g_22 * gamma2_11) / g;
        }
    
    //get stencil constants
    left.zeros(N + 1, M);
    right.zeros(N + 1, M);
    top.zeros(N + 1, M);
    bottom.zeros(N + 1, M);
    center.zeros(N + 1, M);
    
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            right(i, j) = g11(i, j) / SQR(d1) + L1(i, j) / (2 * d1);
            left(i, j)  = g11(i, j) / SQR(d1) - L1(i, j) / (2 * d1);
            top(i, j)   = g22(i, j) / SQR(d2) + L2(i, j) / (2 * d2);
            bottom(i, j)= g22(i, j) / SQR(d2) - L2(i, j) / (2 * d2);
            corner(i, j)= g12(i, j) / (2 * d1 * d2);
            center(i, j)= -(right(i, j) + left(i, j) + top(i, j) + bottom(i, j));
        }
    
    //tridiagnal arrays init
    dd = new double[M];
    bb = new double[M];
}

void solver::advance() {
    double arg1 = 2 / (a * dt);
    
    memcpy(&T(0, 0), &T(N - 1, 0), sizeof(double) * M);
    for (int i = 1; i < N; ++i)
        for (int j = 1; j < M - 1; ++j) {
            Told(i, j) = right(i, j) * T(i + 1, j) + left(i, j) * T(i - 1, j) + top(i, j) * T(i, j + 1) + bottom(i, j) * T(i, j - 1) + (arg1 + center(i, j)) * T(i, j) + corner(i, j) * SEC_MIX(T, i, j);
        }
    
    for(int itr = 0; itr < maxItr; ++itr) {
        tmp = T;
        
        memcpy(&T(0, 0), &T(N - 1, 0), sizeof(double) * M);
        for (int i = 1; i < N; ++i) {
            if (i == N - 1)
                memcpy(&T(N, 0), &T(1, 0), sizeof(double) * M);
            
            for (int j = 1; j < M - 1; ++j) {
                dd[j] = center(i, j) - arg1;
                bb[j] = right(i, j) * T(i + 1, j) + left(i, j) * T(i - 1, j) + corner(i, j) * SEC_MIX(T, i, j) - T(i, j);
            }
            tridiag(&bottom(i, 1), dd, &top(i, 1), bb, &T(i, 1), M - 2);
        }
        
        double error = calcError(tmp, T);
        cout << "Error at Iteration " << itr << " = " << error << endl;
        if (error < tol) {
            cout << "Solution Converged" << endl;
        }
    }
}
