#include <iomanip>
#include <iostream>
#include <fstream>

#include "random.h"
#include "funzionebase.h"

using namespace std;

double Metropolis_step (double pos, double mu, double sigma, FunzioneBase& psi, Random& rnd) {
    
    double step = 2.6;

    double pos_new = rnd.Rannyu(pos-step, pos+step);
        
    double A = min( 1., pow(psi.Eval(pos_new, mu, sigma), 2.)/pow(psi.Eval(pos, mu, sigma), 2.) );
    
    if (rnd.Rannyu() <= A)
        pos = pos_new;
    
    return pos;
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
