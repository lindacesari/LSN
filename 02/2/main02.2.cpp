// excercise 02.2

#include "func.h"

using namespace std;

int main () {
    
    // initializing random generator
    Random rnd;
    set_rnd(rnd);
    
//============================= 02.2.1 ================================
    
    int M = 10000;      // total
    int N = 100;        // number of blocks
    int L = int(M/N);   // length of each block
    int nstep = 100;    // length of each random walk
    
    vector<double> sum(nstep, 0.0);     // initializing vectors of zeros
    vector<double> sum2(nstep, 0.0);
    vector<double> err(nstep, 0.0);
    
    vector<double> distances(nstep, 0.0);   // store distance after _ steps
    
    vector<int> pos;      // stores the current position of the walker on the RW
    
    for (int i = 0; i < N; i++) {
        
        fill(distances.begin(), distances.end(), 0);

        // accumulate |r_N|^2 in distances for each of the L RWs
        for (int j = 0; j < L; j++) {
            
            pos = {0, 0, 0};
            for(int k = 1; k < nstep; k++) {
                pos = step1(pos, rnd);          // new step
                distances[k] += dist(pos);      // accumulate distance after k steps
            }
        }
        
        // accumulate block averages beacuse only the last evaluation is to be plotted
        for (int k = 0; k < nstep; k++) {
            sum[k] += distances[k]/double(L);          // average distances
            sum2[k] += pow(distances[k]/double(L), 2.);
        }
        
    }
    
    // cumulative average, square average and error
    for (int i = 1; i < nstep; i++) {
        sum[i] = sum[i] / double(N);
        sum2[i] = sum2[i] / double(N);
        err[i] = sqrt( (sum2[i] - pow(sum[i], 2)) / double(N) );
    }
    
    ofstream out("RW_discrete.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << M << "," << nstep << endl;
    for (int i = 0; i < nstep; i++)
        out << sqrt(sum[i]) << "," << 0.5*pow(sum[i], -0.5)*err[i] << endl;
    out.close();

//============================= 02.2.2 ================================
    
    // set every element to zero to re-use the same vectors
    fill(sum.begin(), sum.end(), 0);
    fill(sum2.begin(), sum2.end(), 0);
    fill(err.begin(), err.end(), 0);
    
    vector<double> poss;      // stores the current position of the walker on the RW
    
    for (int i = 0; i < N; i++) {
        
        fill(distances.begin(), distances.end(), 0);

        // accumulate |r_N|^2 in distances for each of the L RWs
        for (int j = 0; j < L; j++) {
            
            poss = {0., 0., 0.};
            for(int k = 1; k < nstep; k++) {
                poss = step2(poss, rnd);          // new step
                distances[k] += dist(poss);      // accumulate distance after k steps
            }
        }
        
        // accumulate block averages beacuse only the last evaluation is to be plotted
        for (int k = 0; k < nstep; k++) {
            sum[k] += distances[k]/double(L);          // average distances
            sum2[k] += pow(distances[k]/double(L), 2.);
        }
        
    }
    
    // cumulative average, square average and error
    for (int i = 1; i < nstep; i++) {
        sum[i] = sum[i] / double(N);
        sum2[i] = sum2[i] / double(N);
        err[i] = sqrt( (sum2[i] - pow(sum[i], 2)) / double(N) );
    }
    
    out.open("RW_continuum.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << M << "," << nstep << endl;
    for (int i = 0; i < nstep; i++)
        out << sqrt(sum[i]) << "," << 0.5*pow(sum[i], -0.5)*err[i] << endl;
    out.close();
    
    
//===================================================================
    
    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}
