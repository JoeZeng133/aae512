#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

fstream fo;

void derivative(double *res, double *u) {
    res[0] = - sin(u[1]);
    res[1] = u[0];
}

int main()  {
    fo.open("output.bin", fstream :: out | fstream :: binary);
    
    double delT = 1;
    int num = 10;
    
    double u[2], du[2];
    
    for (int i = 0; i < 9; ++i) {
        u[0] = 0;
        u[1] = 1;
        
        for(int j = 0; j < num; ++j) {
            derivative(du, u);
            u[0] += du[0] * delT;
            u[1] += du[1] * delT;
        }
        
        fo.write((char*)u, sizeof(double) * 2);
        num *= 10;
        delT /= 10;
        
    }
    
    fo.close();
}
