#include <iomanip>
#include <iostream>
#include <fstream>

#include "random.h"
#include "funzionebase.h"

using namespace std;

double Metropolis (int L, vector<double> pos, FunzioneBase& p, Random& rnd, bool uni, bool gs) {
    
    vector<double> pos_1(3, 0.0);
    double step = 1.2;
    double var = 0.75;
    int count = 0;
    double s = 0;

    for (int j = 0; j < L; j++) {

        fill(pos_1.begin(), pos_1.end(), 0.);
    
        // proposed step
        if (gs == true) {
            if (uni == true) {
                step = 1.2;
                pos_1[0] = rnd.Rannyu(pos[0]-step, pos[0]+step);
                pos_1[1] = rnd.Rannyu(pos[1]-step, pos[1]+step);
                pos_1[2] = rnd.Rannyu(pos[2]-step, pos[2]+step);
            }
            else {  //gauss
                var = 0.75;
                pos_1[0] = rnd.Gauss(pos[0], var);
                pos_1[1] = rnd.Gauss(pos[1], var);
                pos_1[2] = rnd.Gauss(pos[2], var);
            }
        }
        else {
            if (uni == true) {
                step = 3;
                pos_1[0] = rnd.Rannyu(pos[0]-step, pos[0]+step);
                pos_1[1] = rnd.Rannyu(pos[1]-step, pos[1]+step);
                pos_1[2] = rnd.Rannyu(pos[2]-step, pos[2]+step);
            }
            else {      // gauss
                var = 1.9;
                pos_1[0] = rnd.Gauss(pos[0], var);
                pos_1[1] = rnd.Gauss(pos[1], var);
                pos_1[2] = rnd.Gauss(pos[2], var);
            }
        }
        
        // acceptance
        double A = min(1., p.Eval(pos_1)/p.Eval(pos));
        double r = rnd.Rannyu();
        if (r <= A) {
            pos = pos_1;
            count++;
        }
        s += sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

    }
    
    // cout << double(count)/L << endl;         // to check the acceptance ratio
    
    return s/L;
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
