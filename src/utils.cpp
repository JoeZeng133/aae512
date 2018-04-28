#include <utils.h>
#include <stdio.h>

using namespace std;

void tridiag(double *aa, double *dd, double *cc, double *bb, double *x, int imax) {
    //
    // Thomas's Tridiagonal Algorithm
    //
    // Description:
    //
    //   Solve a tridiagonal system of the form:
    //
    //   d[0] x[0] + c[0] x[1] = b[0]
    //   a[1] x[0] + d[1] x[1] + c[1] x[2] = b[1]
    //   a[2] x[1] + d[2] x[2] + c[2] x[3] = b[2]
    //   ...
    //   a[N-2] x[N-3] + d[N-2] x[N-2] + c[N-2] x[N-1] = b[N-2]
    //   a[N-1] x[N-2] + d[N-1] x[N-1] = b[N-1]
    //
    // Reference:
    //    Cheney and Kincaid (1994), Numerical Mathematics and Computing,
    //    3rd ed., Brooks/Cole Publishing Co., Pacific Grove CA, 1994,
    //    sec. 6.3, pp. 249-253
    //
    int i;
    double xmult;
    
    for (i = 1; i < imax; i++) {
        xmult = aa[i] / dd[i - 1];
        dd[i] = dd[i] - xmult*cc[i - 1];
        bb[i] = bb[i] - xmult*bb[i - 1];
    }
    
    x[imax - 1] = bb[imax - 1] / dd[imax - 1];
    
    for (i = imax - 2; i >= 0; i--) {
        x[i] = (bb[i] - cc[i] * x[i + 1]) / dd[i];
    }
}

double calcError(const Mat& A, const Mat& B) {
    //p2 norm
    int N, M;
    double arg, diff = 0, norm = 0;
    double res = 1e-16;
    A.getSize(N, M);
    
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) {
            arg = A(i, j) - B(i, j);
            if (A(i, j) != 0)
                res = max(res, abs(arg / A(i, j)));
        }
    return res;
}
