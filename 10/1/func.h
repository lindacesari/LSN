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
          
rowvec move(Popolazione& popu, double T) {
    
    mat pop = popu.get_population();
    
    rowvec cromo = pop.row(0);
    double f = popu.fitness(cromo);
    
    rowvec new_cromo = popu.mutate(cromo);
    double new_f = popu.fitness(new_cromo);
    
    if (new_f <= f)
        cromo = new_cromo;
    else if ( randu() < exp(-1./T*(new_f-f)) )
        cromo = new_cromo;
    
    pop.row(0) = cromo;
    popu.set_population(pop);
    
    return cromo;
}
