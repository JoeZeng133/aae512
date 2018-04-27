//
//  main.cpp
//  mesh
//
//  Created by 曾舟 on 4/21/18.
//

#include <stdio.h>
#include <mesher.h>
#include <solver.h>
#include <string>

using namespace std;

void generateMesh(string bndfile, string meshfile) {
    mesher msh;
    msh.init(bndfile);
    msh.getMesh();
    msh.output(meshfile);
}

void solve(string meshfile, string configfile, string outputfile) {
    solver s;
    Mat xc, yc;
    int N, M;
    fstream fo;
    
    fo.open(meshfile, ios::in | ios::binary);
    if (!fo.is_open())
        throw "Unable to open file " + meshfile;
    
    fo.read((char*)&N, sizeof(int));
    fo.read((char*)&M, sizeof(int));
    xc.zeros(N, M);
    yc.zeros(N, M);
    fo.read((char*)xc.data(), sizeof(double) * xc.getSize());
    fo.read((char*)yc.data(), sizeof(double) * yc.getSize());
    fo.close();

    double T0, T1, T2, a, dt;
    int tstep = 0;
    fo.open(configfile, ios::in);
    if (!fo.is_open())
        throw runtime_error("Unable to open file " + configfile);
    
    fo >> T0 >> T1 >> T2 >> a >> dt >> tstep;
    fo.close();

    s.init(xc, yc, T0, T1, T2, a, dt);

    fo.open(outputfile, ios::out | ios::binary);
    if (!fo.is_open())
        throw runtime_error("Unable to open file " + outputfile);
    
    for (int t = 0; t < tstep; ++t) {
        cout << "Time step " << t + 1 << "\n";
        s.advance();
        s.snapshot(fo);
    }
    
    cout << "Total time steps " << tstep << "\n";
    cout << "Domain size is " << N << " by " << M << "\n";
    fo.close();
}

int main(int argv, char* args[]) {
    if (argv == 3)
        generateMesh(args[1], args[2]);
    else if (argv == 4)
        solve(args[1], args[2], args[3]);
    else
        throw runtime_error("Error, invalid arguments\n");
    
    return 0;
}
