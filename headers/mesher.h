#pragma once

#include<cmath>
#include<iostream>
#include<fstream>
#include<utils.h>

#define FIR_DER_X(x, i, j) (FIR_1(x, i, j) / (2 * d1))
#define FIR_DER_Y(x, i, j) (FIR_2(x, i, j) / (2 * d2))
#define SEC_DER_X(x, i, j) (SEC_1(x, i, j) / (d1 * d1))
#define SEC_DER_Y(x, i, j) (SEC_2(x, i, j) / (d2 * d2))


class mesher {
private:
    int clustering;
    int N, M;
    int maxItr = 200;
    double d2, d1;
    double tol = 1e-5;
    
    double omega;
    double* aa, *dd, *cc, *bbx, *bby, *tmp, *xtmp, *ytmp;
    
    Mat xc, yc, alpha, gamma, beta, psi, phi, xprev, yprev;
    
    std::fstream fo;
    
public:
    
    void init(std::string filename);
    void getMesh();
    void getPQ();
    void updateJacobian();
    void slor(int i);
    void output(std::string filename);
    Mat getxc();
    Mat getyc();
};

Mat mesher::getxc() {
    return xc;
}

Mat mesher::getyc() {
    return yc;
}

void mesher::slor(int i) {
    double arg1, arg2, arg3, arg4;
    double tau = d2 / d1;
    
    //build equations
    for (int j = 1; j < M - 1; ++j) {
        arg1 = 0.5 * gamma(i, j) * psi(i, j) * d2;
        arg2 = alpha(i, j) * tau * tau;
        arg3 = 0.5 * alpha(i, j) * phi(i, j) * tau * d2;
        arg4 = 0.5 * beta(i, j) * tau;
        
        aa[j] = gamma(i, j) * (1 - arg1);
        dd[j] = -2 * (gamma(i, j) + arg2);
        cc[j] = gamma(i, j) * (1 + arg1);
        
        //cout << MIXED(yc, j, M) << " "  << MIXED(xc, j, M) << "\n";
        
        bbx[j] = -arg2 * AVE_1(xc, i, j) - arg3 * FIR_1(xc, i, j) + arg4 * SEC_MIX(xc, i, j);
        bby[j] = -arg2 * AVE_1(yc, i, j) - arg3 * FIR_1(yc, i, j) + arg4 * SEC_MIX(yc, i, j);
    }
    
    //assign boundary conditions
    cc[0] = 0;
    aa[M - 1] = 0;
    dd[0] = dd[M - 1] = 1;
    bbx[0] = xc(i, 0);
    bbx[M - 1] = xc(i, M - 1);
    bby[0] = yc(i, 0);
    bby[M - 1] = yc(i, M - 1);
    
    /*    for (int j = 0; j < M; ++j) {
     printf("%5.3e %5.3e %5.3e %5.3e %5.3e\n", aa[j], dd[j], cc[j], bbx[j], bby[j]);
     }*/
    
    //line gaussian-seidel
    memcpy(tmp, dd, sizeof(double) * M);
    tridiag(aa, dd, cc, bbx, xtmp, M);
    
    memcpy(dd, tmp, sizeof(double) * M);
    tridiag(aa, dd, cc, bby, ytmp, M);
    
    //over relaxation
    for (int j = 1; j < M - 1; ++j) {
        xc(i, j) = omega * xtmp[j] + (1 - omega) * xc(i, j);
        yc(i, j) = omega * ytmp[j] + (1 - omega) * yc(i, j);
    }
}

void mesher::updateJacobian() {
    double x1, x2, y1, y2;
    
    for (int i = 1; i < N - 1; ++i)
        for (int j = 1; j < M - 1; ++j) {
            x1 = FIR_DER_X(xc, i, j);
            x2 = FIR_DER_Y(xc, i, j);
            y1 = FIR_DER_X(yc, i, j);
            y2 = FIR_DER_Y(yc, i, j);
            
            //printf("%5.3e %5.3e %5.3e %5.3e\n", x1, x2, y1, y2);
            alpha(i, j) = x2 * x2 + y2 * y2;
            beta(i, j) = x1 * x2 + y1 * y2;
            gamma(i, j) = x1 * x1 + y1 * y1;
            
            //printf("%5.3e %5.3e %5.3e\n", alpha(i, j), beta(i, j), gamma(i, j));
        }
}


void mesher::getPQ() {
    double xp, xpp, yp, ypp;
    
    for (int i = 1; i < N - 1; ++i) {
        xp = FIR_DER_X(xc, i, M - 1);
        xpp = SEC_DER_X(xc, i, M - 1);
        yp = FIR_DER_X(yc, i, M - 1);
        ypp = SEC_DER_X(yc, i, M - 1);
        phi(i, M - 1) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);
        
        xp = FIR_DER_X(xc, i, 0);
        xpp = SEC_DER_X(xc, i, 0);
        yp = FIR_DER_X(yc, i, 0);
        ypp = SEC_DER_X(yc, i, 0);
        phi(i, 0) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);
        
    }
    
    for (int j = 1; j < M - 1; ++j) {
        xp = FIR_DER_Y(xc, N - 1, j);
        xpp = SEC_DER_Y(xc, N - 1, j);
        yp = FIR_DER_Y(yc, N - 1, j);
        ypp = SEC_DER_Y(yc, N - 1, j);
        psi(N - 1, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);
        
        xp = FIR_DER_Y(xc, 0, j);
        xpp = SEC_DER_Y(xc, 0, j);
        yp = FIR_DER_Y(yc, 0, j);
        ypp = SEC_DER_Y(yc, 0, j);
        psi(0, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);
        
        //printf("%9.3e %9.3e %9.3e %9.3e\n", xp, xpp, yp, ypp);
    }
    
    //linear interpolation
    for (int i = 1; i < N - 1; ++i)
        for (int j = 1; j < M - 1; ++j) {
            psi(i, j) = (i * psi(N - 1, j) + (N - 1 - i) * psi(0, j)) / (N - 1);
            phi(i, j) = (j * phi(i, M - 1) + (M - 1 - j) * phi(i, 0)) / (M - 1);
        }
}



void mesher::init(std::string filename) {
    fo.open(filename, std::ios::binary | std::ios::in);
    if (!fo.is_open()) {
        std::cout << "Fail to open " << filename << "\n";
        exit(1);
    }
    fo.read((char*)&clustering, sizeof(int));
    fo.read((char*)&N, sizeof(int));
    fo.read((char*)&M, sizeof(int));
    xc.rand(N, M);
    yc.rand(N, M);
    
    //read boundary point coordinates
    for (int i = 0; i < N - 1; ++i) {
        fo.read((char*)&xc(i, 0), sizeof(double));
        fo.read((char*)&yc(i, 0), sizeof(double));
    }
    
    for (int j = 0; j < M - 1; ++j) {
        fo.read((char*)&xc(N - 1, j), sizeof(double));
        fo.read((char*)&yc(N - 1, j), sizeof(double));
    }
    
    for (int i = N - 1; i > 0; --i) {
        fo.read((char*)&xc(i, M - 1), sizeof(double));
        fo.read((char*)&yc(i, M - 1), sizeof(double));
    }
    
    for (int j = M - 1; j > 0; --j) {
        fo.read((char*)&xc(0, j), sizeof(double));
        fo.read((char*)&yc(0, j), sizeof(double));
    }
    
    fo.close();
    
    d2 = 1. / (M - 1);
    d1 = 1. / (N - 1);
    
    //governing equations coefficients
    psi.zeros(N, M);
    phi.zeros(N, M);
    alpha.ones(N, M);
    gamma.ones(N, M);
    beta.ones(N, M);
    
    //slor buffer arrays
    omega = 1.7;
    aa = new double[M];
    dd = new double[M];
    cc = new double[M];
    bbx = new double[M];
    bby = new double[M];
    tmp = new double[M];
    xtmp= new double[M];
    ytmp = new double[M];
}

void mesher::getMesh() {
    for(int itr = 0; itr < maxItr; ++itr) {
        updateJacobian();
        xprev = xc;
        yprev = yc;
        //Successive Line Over Relaxation
        for (int i = 1; i < N - 1; ++i) {
            slor(i);
        }
        
        double xerror = calcError(xprev, xc);
        double yerror = calcError(yprev, yc);
        std::cout << "At iteration" << itr << ", the error is : xerror = " << xerror << ", yerror = " << yerror << "\n";
        if (xerror < tol && yerror < tol) {
            std::cout << "Solution Converged" << "\n";
            break;
        }
    }
    
    if (clustering) {
        getPQ();
        for (int itr = 0; itr < maxItr; ++itr) {
            updateJacobian();
            xprev = xc;
            yprev = yc;
            //Successive Line Over Relaxation
            for (int i = 1; i < N - 1; ++i) {
                slor(i);
            }
            
            double xerror = calcError(xprev, xc);
            double yerror = calcError(yprev, yc);
            std::cout << "At iteration" << itr << ", the error is : xerror = " << xerror << ", yerror = " << yerror << "\n";
            if (xerror < tol && yerror < tol) {
                std::cout << "Solution Converged" << "\n";
                break;
            }
        }
    }
    fo.close();
}

void mesher::output(std::string filename) {
    fo.open(filename, std::ios::out | std::ios::binary);
    if (!fo.is_open()) {
        std::cout << "Fail to open" << filename << "\n";
        exit(1);
    }
    
    std::cout << xc.getSize() << " " << yc.getSize() << "\n";
    fo.write((char*)xc.data(), sizeof(double) * xc.getSize());
    fo.write((char*)yc.data(), sizeof(double) * yc.getSize());
    fo.close();
}
