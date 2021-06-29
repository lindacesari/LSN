// excercise 01.2

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
    
//============================= 01.2 ================================
    
    vector<int> N{1, 2, 10, 100};       // values used to compute average
    double s_unif = 0;
    double s_exp = 0;
    double s_cl = 0;
    
    ofstream out("dice.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    
    for (int j = 0; j < N.size(); j++) {
        for (int h = 0; h < 10000; h++) {
            s_unif = 0;
            s_exp = 0;
            s_cl = 0;
            for (int i = 0; i < N[j]; i++) {
                s_unif += rnd.Rannyu();
                s_exp += rnd.Exponential(1.);
                s_cl += rnd.Cauchy_Lorentz(1., 0.);
            }
            out << s_unif/N[j] << "," << s_exp/N[j] << "," << s_cl/N[j] << endl;
        }
    }
    
    out.close();

//===================================================================
    
    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}


