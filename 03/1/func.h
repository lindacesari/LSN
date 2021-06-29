#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>

#include "random.h"

using namespace std;

double C_i(double S0, double r, double sigma, double T, double K, Random& rnd) {
    
    double z = rnd.Gauss(0., 1.);
    double S = S0 * exp( (r-sigma*sigma/2.)*T + sigma*z*pow(T, 0.5) );
    return exp(-r*T)*max(0., S-K);
}

double P_i(double S0, double r, double sigma, double T, double K, Random& rnd) {
    
    double z = rnd.Gauss(0., 1.);
    double S = S0 * exp( (r-sigma*sigma/2.)*T + sigma*z*pow(T, 0.5) );
    return exp(-r*T)*max(0., K-S);
}

double C_i_GBM(double S0, double r, double sigma, double T, double K, Random& rnd) {
    
    double z;
    double S = S0;
    for (int i = 0; i < 100; i++) {
        z = rnd.Gauss(0., 1.);
        S = S * exp((r-0.5*sigma*sigma)*0.01 + sigma*z*pow(0.01, 0.5));
    }
    return exp(-r*T)*max(0., S-K);
}

double P_i_GBM(double S0, double r, double sigma, double T, double K, Random& rnd) {
    
    double z;
    double S = S0;
    for (int i = 0; i < 100; i++) {
        z = rnd.Gauss(0., 1.);
        S = S * exp((r-0.5*sigma*sigma)*0.01 + sigma*z*pow(0.01, 0.5));
    }
    return exp(-r*T)*max(0., K-S);
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
