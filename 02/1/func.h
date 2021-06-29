#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include "random.h"
#include "funzionebase.h"

using namespace std;

// compute integral with method of the average from a uniform distribution
double IntegraleAVE(double xmin, double xmax, FunzioneBase &f, int npunti, Random& rnd) {

    double sum = 0;
    double x = 0;
    for (int i = 0; i < npunti; i++) {
        x = rnd.Rannyu(xmin, xmax);
        sum = sum + f.Eval(x);
    }

    return (sum * (xmax-xmin)/(double)npunti);
}

// compute integral with method of the average from importance sampling
double IntegraleAVE_IS(double xmin, double xmax, double dmax, FunzioneBase &f, FunzioneBase &d, int npunti, Random &rnd) {

    double sum = 0;
    double x = 0;
    for (int i = 0; i < npunti; i++) {
        x = rnd.Reject(xmin, xmax, dmax, d);
        sum = sum + f.Eval(x);
    }

    return (sum * (xmax-xmin)/(double)npunti);
}


double error(vector<double> ave, vector<double> ave2, int n) {
    
    if (n == 0) {
        return 0.;
    }
    else {
        return sqrt( (ave2[n] - pow(ave[n], 2)) / n );
    }
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
