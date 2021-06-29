#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include "random.h"

using namespace std;

double error(vector<double> ave, vector<double> ave2, int n) {
    if (n == 0) {
        return 0.;
    }
    else {
        return sqrt( (ave2[n] - pow(ave[n], 2)) / n );
    }
}

bool Intersect(Random &rnd, double L, double d) {
    
    double y1 = rnd.Rannyu();       // first random extreme of the needle
    
    double u = rnd.Rannyu(-1., 1.); // values to sample sin(theta)
    double v = rnd.Rannyu(0., 1.);
    while(pow(u,2)+pow(v,2) > 1) {
        u = rnd.Rannyu(-1., 1.);
        v = rnd.Rannyu(0., 1.);
    }
    
    double y2 = y1 + L*(2*u*v/(pow(u,2)+pow(v,2))); // y2=y1+Lsin(theta)
    
    // check if the needle intersected a line at y=d or y=0
    if ((y1>=d && y2<=d) || (y1<=d && y2>=d) || (y1>=0 && y2<=0) || (y1>=0 && y2<=0))
        return true;
    else
        return false;
}

double pi_esteem (Random &rnd, double L, double d, int throws) {
    
    // count how many of the N_throws are successful
    int n = 0;
    for (int j = 0; j < throws; j++) {
        n += int(Intersect(rnd, L, d));
    }
    
    // compute pi as (2*L*N_throws)/(N_hit*d)
    return 2*L*throws/(double(n)*d);
}

void set_rnd(Random& rnd) {
    
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    return;
}
