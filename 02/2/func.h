#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <vector>

#include "random.h"

using namespace std;

// step of length = 1 on a cubic lattice (discrete)
vector<int> step1(vector<int> pos, Random &rnd) {
    
    int n = int(rnd.Rannyu(0, 3));      // choose direction: either 0=x, 1=y or 2=z
    pos[n] += 2*int(rnd.Rannyu(0, 2)) - 1;    // choose step: either +1=forward or -1=backward
    /*
    if (x >= 0)
        step[n] += 1;
    else
        step[n] -= 1;*/
    
    return pos;
}

// step of length = 1 in a random direction in continuum
vector<double> step2(vector<double> pos, Random &rnd) {
    
    double theta = rnd.Rannyu(0, M_PI);     // choose direction (r=1, theta, phi)
    double phi = rnd.Rannyu(0, 2*M_PI);
    
    pos[0] += cos(phi)*sin(theta);          // x = cos(phi)sin(theta)
    pos[1] += sin(phi)*sin(theta);          // y = sin(phi)sin(theta)
    pos[2] += cos(theta);                   // z = cos(theta)
   
    return pos;
}

// compute squared distance from origin (0,0,0)
double dist(vector<int> pos) {
    return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
}

double dist(vector<double> pos) {
    return pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
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
