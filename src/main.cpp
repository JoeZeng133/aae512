//
//  main.cpp
//  mesh
//
//  Created by 曾舟 on 4/21/18.
//

#include <stdio.h>
#include <mesher.h>
#include <solver.h>

using namespace std;

int main() {
    solver s;
    Mat xc, yc;
    int N, M;
    
    fstream fo;
    fo.open("../output/meshtest.bin", ios::in | ios::binary);
    fo.read((char*)&N, sizeof(int));
    fo.read((char*)&M, sizeof(int));
    cout << "Domain size is " << N << " by " << M << "\n";
    xc.zeros(N, M);
    yc.zeros(N, M);
    fo.read((char*)xc.data(), sizeof(double) * xc.getSize());
    fo.read((char*)yc.data(), sizeof(double) * yc.getSize());
    fo.close();
    
    s.init(xc, yc, 0, 100, 20, 6.4241e-5, 0.1);
    
    for (int t = 0; t < 1000; ++t) {
        cout << "Time step " << t + 1 << "\n";
        s.advance();
    }
    
    s.snapshot();
    return 0;
}
