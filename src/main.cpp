#include<cmath>
#include<iostream>
#include<fstream>
#include<utils.h>

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
#define FIR_DER_X(x, i, j) ((x(i + 1, j) - x(i - 1, j)) / 2 / ksid)
#define FIR_DER_Y(x, i, j) ((x(i, j + 1) - x(i, j - 1)) / 2 / etad)
#define SEC_DER_X(x, i, j) ((x(i + 1, j) + x(i - 1, j) - 2 * x(i, j)) / (ksid * ksid))
#define SEC_DER_Y(x, i, j) ((x(i, j + 1) + x(i, j - 1) - 2 * x(i, j)) / (etad * etad))


namespace slor {
#define MIXED(x, j, N) (x[j + N + 1] + x[j - N - 1] - x[j + N - 1] - x[j - N + 1])
#define AVE(x, j, N) (x[j + N] + x[j - N])
#define DIF(x, j, N) (x[j + N] - x[j - N])

	int N, M;
	double omega;
	double* aa, *dd, *cc, *bbx, *bby, *tmp, *xtmp, *ytmp;

	void init(int _N, int _M) {
		N = _N;
		M = _M;
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

	void destructor() {
		delete[] aa;
		delete[] dd;
		delete[] cc;
		delete[] bbx;
		delete[] bby;
		delete[] tmp;
		delete[] xtmp;
		delete[] ytmp;
	}

	void slor(double etad, double ksid, double* xc, double *yc, double* alpha, double* beta, double* gamma, double* phi, double* psi) {
		double arg1, arg2, arg3, arg4;
		double tau = etad / ksid;

		//build equations
		for (int j = 1; j < M - 1; ++j) {
			arg1 = 0.5 * gamma[j] * psi[j] * etad;
			arg2 = alpha[j] * tau * tau;
			arg3 = 0.5 * alpha[j] * phi[j] * tau * etad;
			arg4 = 0.5 * beta[j] * tau;

			aa[j] = gamma[j] * (1 - arg1);
			dd[j] = -2 * (gamma[j] + arg2);
			cc[j] = gamma[j] * (1 + arg1);

			//cout << MIXED(yc, j, M) << " "  << MIXED(xc, j, M) << endl;

			bbx[j] = -arg2 * AVE(xc, j, M) - arg3 * DIF(xc, j, M) + arg4 * MIXED(xc, j, M);
			bby[j] = -arg2 * AVE(yc, j, M) - arg3 * DIF(yc, j, M) + arg4 * MIXED(yc, j, M);
		}

		//assign boundary conditions
		cc[0] = 0;
		aa[M - 1] = 0;
		dd[0] = dd[M - 1] = 1;
		bbx[0] = xc[0];
		bbx[M - 1] = xc[M - 1];
		bby[0] = yc[0];
		bby[M - 1] = yc[M - 1];

	/*	for (int j = 0; j < M; ++j) {
			printf("%5.3e %5.3e %5.3e %5.3e %5.3e\n", aa[j], dd[j], cc[j], bbx[j], bby[j]);
		}*/

		//line gaussian-seidel
		memcpy(tmp, dd, sizeof(double) * M);
		tridiag(aa, dd, cc, bbx, xtmp, M);

		memcpy(dd, tmp, sizeof(double) * M);
		tridiag(aa, dd, cc, bby, ytmp, M);

		//over relaxation
		for (int i = 1; i < M - 1; ++i) {
			xc[i] = omega * xtmp[i] + (1 - omega) * xc[i];
			yc[i] = omega * ytmp[i] + (1 - omega) * yc[i];
		}
	}
}

void updateJacobian(const Mat& xc, const Mat& yc, Mat& alpha, Mat& beta, Mat& gamma, double etad, double ksid) {
	int N, M;
	double x1, x2, y1, y2;
	xc.getSize(N, M);

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


void getPQ(const Mat& xc, const Mat& yc, Mat& psi, Mat& phi, double etad, double ksid) {
	int N, M;
	double xp, xpp, yp, ypp;

	xc.getSize(N, M);

	//get boundary values
	/*for (int i = 1; i < N - 1; ++i)
		for (int j = 1; j < M - 1; ++j) {
			xp = FIR_DER_X(xc, i, j);
			xpp = SEC_DER_X(xc, i, j);
			yp = FIR_DER_X(yc, i, j);
			ypp = SEC_DER_X(yc, i, j);
			phi(i, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);

			xp = FIR_DER_X(xc, i, j);
			xpp = SEC_DER_X(xc, i, j);
			yp = FIR_DER_X(yc, i, j);
			ypp = SEC_DER_X(yc, i, j);
			phi(i, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);

			xp = FIR_DER_Y(xc, i, j);
			xpp = SEC_DER_Y(xc, i, j);
			yp = FIR_DER_Y(yc, i, j);
			ypp = SEC_DER_Y(yc, i, j);
			psi(i, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);

			xp = FIR_DER_Y(xc, i, j);
			xpp = SEC_DER_Y(xc, i, j);
			yp = FIR_DER_Y(yc, i, j);
			ypp = SEC_DER_Y(yc, i, j);
			psi(i, j) = -(xp * xpp + yp * ypp) / (xp * xp + yp * yp);
		}*/

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

double calcError(const Mat& A, const Mat& B) {
	//p2 norm
	int N, M;
	double arg, diff = 0, norm = 0;
	A.getSize(N, M);


	for (int i = 0; i < N; ++i)
		for (int j = 0; j < M; ++j) {
			arg = A(i, j) - B(i, j);
			diff += arg * arg;
			norm += A(i, j) * A(i, j);
		}

	return sqrt(diff / norm);
}

void read(fstream& fo, Mat& xc, Mat& yc, int& N, int& M, int& clustering) {
	fo.read((char*)&clustering, sizeof(int));
	fo.read((char*)&N, sizeof(int));
	fo.read((char*)&M, sizeof(int));
	xc.rand(N, M);
	yc.rand(N, M);

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
}



int main() {
	int clustering;
	int N, M;
	int maxItr = 200;
	double etad, ksid;
	double tol = 1e-5;
	fstream fo;
	Mat xc, yc, alpha, gamma, beta, psi, phi, xprev, yprev;

	fo.open("../output/bnd.bin", ios::in | ios::binary);
	if (!fo.is_open()) {
		std::cout << "Fail to Open File" << endl;
		exit(1);
	}
	read(fo, xc, yc, N, M, clustering);
	etad = 1. / (M - 1);
	ksid = 1. / (N - 1);

	psi.zeros(N, M);
	phi.zeros(N, M);
	alpha.ones(N, M);
	gamma.ones(N, M);
	beta.ones(N, M);

	slor::init(N, M);
	for(int itr = 0; itr < maxItr; ++itr) {
		updateJacobian(xc, yc, alpha, beta, gamma, etad, ksid);
		xprev = xc;
		yprev = yc;
		//Successive Line Over Relaxation
		for (int i = 1; i < N - 1; ++i) {
			slor::slor(etad, ksid, &xc(i, 0), &yc(i, 0), &alpha(i, 0), &beta(i, 0), &gamma(i, 0), &phi(i, 0), &psi(i, 0));
		}
		
		double xerror = calcError(xprev, xc);
		double yerror = calcError(yprev, yc);
		std::cout << "At iteration" << itr << ", the error is : xerror = " << xerror << ", yerror = " << yerror << endl;
		if (xerror < tol && yerror < tol) {
			std::cout << "Solution Converged" << endl;
			break;
		}
	}

	if (clustering) {
		getPQ(xc, yc, psi, phi, etad, ksid);
		for (int itr = 0; itr < maxItr; ++itr) {
			updateJacobian(xc, yc, alpha, beta, gamma, etad, ksid);
			xprev = xc;
			yprev = yc;
			//Successive Line Over Relaxation
			for (int i = 1; i < N - 1; ++i) {
				slor::slor(etad, ksid, &xc(i, 0), &yc(i, 0), &alpha(i, 0), &beta(i, 0), &gamma(i, 0), &phi(i, 0), &psi(i, 0));
			}

			double xerror = calcError(xprev, xc);
			double yerror = calcError(yprev, yc);
			std::cout << "At iteration" << itr << ", the error is : xerror = " << xerror << ", yerror = " << yerror << endl;
			if (xerror < tol && yerror < tol) {
				std::cout << "Solution Converged" << endl;
				break;
			}
		}
	}
	

	fo.close();


	fo.open("../output/mesh.bin", ios::out | ios::binary);
	std::cout << xc.getSize() << " " << yc.getSize() << endl;
	fo.write((char*)xc.data(), sizeof(double) * xc.getSize());
	fo.write((char*)yc.data(), sizeof(double) * yc.getSize());
	fo.close();

	system("PAUSE");
	return 0;
}