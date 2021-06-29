// excercise 01.3: Buffon's experiment

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include "random.h"
#include "func.h"

using namespace std;

int main () {
    
    // initializing random generator
    Random rnd;
    set_rnd(rnd);
    
//============================= 01.3 ================================
    
    // Choose two horizontal lines at y=0 and y=d (the other lines are reachable via translation). Generate a random number for the first extreme of the needle, generate an angle theta and calculate the other extreme. Find if a line is between these two extremes.
    
    int M = 100000;
    int N = 100;
    
    double L = 0.8;     // length of the needle
    double d = 1.;      // distance between horizontal lines
    double x;
    vector<double> pi(N, 0.0);
    vector<double> pi2(N, 0.0);
    vector<double> err(N, 0.0);

    for (int i = 0; i < N; i++) {
        
        // esteem of pi with M/N throws
        x = pi_esteem(rnd, L, d, M/N);
        
        // progressive sums of pi and pi^2
        if (i == 0) {
            pi[i] = x;
            pi2[i] = pow(x, 2);
        }
        else {
            pi[i] = (pi[i-1] + x);
            pi2[i] = (pi2[i-1] + pow(x, 2));
        }
    }
    
    // cumulative average, square average and error
    for (int i = 0; i < N; i++) {
        pi[i] = pi[i] / (i+1);
        pi2[i] = pi2[i] / (i+1);
        err[i] = error(pi, pi2, i);
    }
    
    // writing computed data on file
    ofstream out("pi.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
        out << pi[i] << "," << err[i] << endl;
    out.close();
    
    
//===================================================================
    
    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}


