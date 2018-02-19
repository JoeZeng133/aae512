#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#define spaceSteps 20
#define dx 0.05
#define SQR(x) ((x) * (x))

void tridiag(double *aa, double *dd, double *cc, double *bb, double *x, int imax) {
    int i;
    double xmult;

    for (i = 1; i < imax; ++i) {
        xmult = aa[i] / dd[i - 1];
        dd[i] = dd[i] - xmult * cc[i - 1];
        bb[i] = bb[i] - xmult * bb[i - 1];
    }

    x[imax - 1] = bb[imax - 1] / dd[imax - 1];

    for (i = imax - 2; i >= 0; --i) {
        x[i] = (bb[i] - cc[i] * x[i + 1]) / dd[i];
    }
}

void ftcs(string filename, double r) {
    fstream fo;
    fo.open(filename, ios :: out | ios :: binary);
    double endTime = 1.0 / 16 / 1e-6;
    double dt = r * SQR(dx) / 1e-6;

    int timeSteps = endTime / dt;
    int now = 0,prev = 1;

    vector<double> u[2];
    u[0].resize(spaceSteps + 1);
    u[1].resize(spaceSteps + 1);

    fill(u[prev].begin(), u[prev].end(), 0);
    for(int t = 1; t <= timeSteps; ++t) {
        prev = now ^ 1;
        //force boundary conditions
        u[now].front() = 1;
        u[now].back() = 0;
        //explicit update
        for(int j = 1; j < spaceSteps; ++j) {
            u[now][j] = u[prev][j] + r * (u[prev][j + 1] - 2 * u[prev][j] + u[prev][j - 1]);
        }
        //alternate arrays
        now = now ^ 1;
    }

    fo.write((char*)&u[now ^ 1][0], sizeof(double) * (spaceSteps + 1));
    fo.close();
}

void btcs(string filename, double r) {
    fstream fo;
    fo.open(filename, ios :: out | ios :: binary);
    double endTime = 1.0 / 16 / 1e-6;
    double dt = r * SQR(dx) / 1e-6;

    int timeSteps = endTime / dt;
    int now = 0,prev = 1;

    vector<double> u[2];
    vector<double> a(spaceSteps - 1);
    vector<double> d(spaceSteps - 1);
    vector<double> c(spaceSteps - 1);

    u[0].resize(spaceSteps + 1);
    u[1].resize(spaceSteps + 1);

    fill(u[prev].begin(), u[prev].end(), 0);
    for(int t = 1; t <= timeSteps; ++t) {
        prev = now ^ 1;
        u[now].front() = 1;
        u[now].back() = 0;

        fill(a.begin(), a.end(), -r);
        fill(d.begin(), d.end(), 2 * r + 1);
        fill(c.begin(), c.end(), -r);

        //add boundary term and solve the implicit equation
        u[prev][1] += r * 1;
        tridiag(&a[0], &d[0], &c[0], &u[prev][1], &u[now][1], spaceSteps - 1);
        //alternate arrays
        now = now ^ 1;
    }

    fo.write((char*)&u[now ^ 1][0], sizeof(double) * (spaceSteps + 1));
    fo.close();
}

int main (int argv, char* args[]) {
    ftcs("explicit_r=0.25.bin", 0.25);
    ftcs("explicit_r=0.52.bin", 0.52);

    btcs("implicit_r=0.25.bin", 0.25);
    btcs("implicit_r=0.52.bin", 0.52);
}
