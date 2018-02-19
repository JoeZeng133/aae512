#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <vector>

#define spaceSteps 70
#define dx 1
#define SQR(x) ((x) * (x))

using namespace std;

typedef vector<double> vec;

double initCondition(double x) {
    if (x < 0 || x > 70) {
        cout << "Outside range [0, 70]" << endl;
        exit(1);
    }

    if (x <= 5)
        return 0;
    else if (x <= 15)
        return 2 * x - 10;
    else if (x <= 25)
        return -2 * x + 50;
    else
        return 0;
}

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

void laxWendroff (string filename, double dt) {
    fstream fo;
    fo.open(filename, ios :: out | ios :: binary);

    const int timeSteps = 0.15 / dt;
    const double c1 = 200 * dt / 2 / dx;
    const double c2 = SQR(c1) * 2;
    int now = 0, prev = 1;

    vector<double> u[2];
    u[0].resize(spaceSteps + 1);
    u[1].resize(spaceSteps + 1);

    //initializing
    for (int i = 0; i <= spaceSteps; ++i)
        u[prev][i] = initCondition(i * dx);

    //stepping
    for (int t = 1; t <= timeSteps; ++t) {
        prev = now ^ 1;
        //force boundary conditions
        u[now][0] = u[now][spaceSteps] = 0;
        //explicit one-step
        for(int j = 1; j < spaceSteps; ++j)
            u[now][j] = u[prev][j] - c1 * (u[prev][j + 1] - u[prev][j - 1])
            + c2 * (u[prev][j + 1] + u[prev][j - 1] - 2 * u[prev][j]);
        //alternate arrays
        now ^= 1;
    }
    //writing to file
    fo.write((char*)&u[now ^ 1][0], sizeof(double) * (spaceSteps + 1));
    fo.close();
}

void implicitEuler(string filename, double dt) {
    fstream fo;
    fo.open(filename, ios :: out | ios :: binary);

    const int timeSteps = 0.15 / dt;
    const double c1 = 200 * dt / dx;
    int now = 0, prev = 1;

    vec u[2];
    u[0].resize(spaceSteps + 1);
    u[1].resize(spaceSteps + 1);
    vec a(spaceSteps - 1);
    vec d(spaceSteps - 1);
    vec c(spaceSteps - 1);

    for (int i = 0; i <= spaceSteps; ++i)
        u[prev][i] = initCondition(i * dx);

    for (int t = 1; t <= timeSteps; ++t) {
        prev = now ^ 1;
        u[now][0] = u[now][spaceSteps] = 0;
        fill(a.begin(), a.end(), - c1 / 2);
        fill(d.begin(), d.end(), 1);
        fill(c.begin(), c.end(), c1 / 2);

        //tridiagnal solve
        tridiag(&a[0], &d[0], &c[0], &u[prev][1], &u[now][1], spaceSteps - 1);
        //force boundary conditions
        u[now].front() = u[now].back() = 0;
        //alternate array
        now ^= 1;
    }
    //write to file
    fo.write((char*)&u[now ^ 1][0], sizeof(double) * (spaceSteps + 1));
    fo.close();
}

int main(int argc, char* args[]) {
    laxWendroff("Lax_dt=2.5e-3.bin", 2.5e-3);
    laxWendroff("Lax_dt=5e-3.bin", 5e-3);

    implicitEuler("Implicit_dt=2.5e-3.bin", 2.5e-3);
    implicitEuler("Implicit_dt=5e-3.bin", 5e-3);
}
