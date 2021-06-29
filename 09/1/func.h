#include <iomanip>
#include <iostream>
#include <fstream>
#include <armadillo>

#include "popolazione.h"

using namespace std;
using namespace arma;

mat make_positions(bool circ, int n_cities) {
    
    mat positions(n_cities, 2);
    
    if (circ == true) {     // on a circumference, r = 1
        double theta;
        for (int i = 0; i < n_cities; i++) {
            theta = 2*M_PI*randu();
            positions.row(i).col(0) = cos(theta);
            positions.row(i).col(1) = sin(theta);
        }
    } else {                // in a square, [-1, 1]x[-1, 1]
        double lato = 2.;
        for (int i = 0; i < n_cities; i++) {
            positions.row(i).col(0) = randu()*lato - 1.;
            positions.row(i).col(1) = randu()*lato - 1.;
        }
    }
    
    return positions;
}
          
double error(vector<double> ave, vector<double> ave2, int n) {
    
    if (n == 0) {
        return 0.;
    }
    else {
        return sqrt( (ave2[n] - pow(ave[n], 2)) / n );
    }
}

