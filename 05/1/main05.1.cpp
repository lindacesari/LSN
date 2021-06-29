// excercise 03.1

#include "func.h"

using namespace std;

int main () {
    
    // initializing random generator
    Random rnd;
    set_rnd(rnd);
    
//============================= 03.1 ================================
    
    probabilita1 p1;
    probabilita2 p2;
    
    int M = 1.0e6;     // total
    int N = 100;        // number of blocks
    int L = int(M/N);   // length of each block
    
    vector<double> sum(N, 0.0); // initializing vectors of zeros
    vector<double> sum2(N, 0.0);
    vector<double> err(N, 0.0);
    
    vector<double> pos = {0.5, 0., 0.};
    
    // compute average for each block and store it in the progressive sums vectors
    double s;

// Uniform p1si
    
    for (int i = 0; i < N; i++) {
        
        s = 0.;
        s = Metropolis(L, pos, p1, rnd, true, true);
        
        if (i == 0) {
            sum[i] = s;
            sum2[i] = pow(s, 2);
        }
        else {
            sum[i] = (sum[i-1] + s);
            sum2[i] = (sum2[i-1] + pow(s, 2));
        }
    }
    
    // cumulative average, square average and error
    for (int i = 0; i < N; i++) {
        sum[i] = sum[i] / (i+1);
        sum2[i] = sum2[i] / (i+1);
        err[i] = error(sum, sum2, i);
    }
    
    // writing computed data on file
    ofstream out("unif_p1.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();

// Uniform psi2
    
    fill(sum.begin(), sum.end(), 0);    // set every element to zero to re-use the same vectors as above
    fill(sum2.begin(), sum2.end(), 0);
    fill(err.begin(), err.end(), 0);

    for (int i = 0; i < N; i++) {
        
        s = 0.;
        s = Metropolis(L, pos, p2, rnd, false, true);
        
        if (i == 0) {
            sum[i] = s;
            sum2[i] = pow(s, 2);
        }
        else {
            sum[i] = (sum[i-1] + s);
            sum2[i] = (sum2[i-1] + pow(s, 2));
        }
    }
    
    // cumulative average, square average and error
    for (int i = 0; i < N; i++) {
        sum[i] = sum[i] / (i+1);
        sum2[i] = sum2[i] / (i+1);
        err[i] = error(sum, sum2, i);
    }
    
    // writing computed data on file
    out.open("unif_p2.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();

// Gaussian psi1

    fill(sum.begin(), sum.end(), 0);    // set every element to zero to re-use the same vectors as above
    fill(sum2.begin(), sum2.end(), 0);
    fill(err.begin(), err.end(), 0);

    for (int i = 0; i < N; i++) {
        
        s = 0.;
        s = Metropolis(L, pos, p1, rnd, true, false);
        if (i == 0) {
            sum[i] = s;
            sum2[i] = pow(s, 2);
        }
        else {
            sum[i] = (sum[i-1] + s);
            sum2[i] = (sum2[i-1] + pow(s, 2));
        }
    }
    
    // cumulative average, square average and error
    for (int i = 0; i < N; i++) {
        sum[i] = sum[i] / (i+1);
        sum2[i] = sum2[i] / (i+1);
        err[i] = error(sum, sum2, i);
    }
    
    // writing computed data on file
    out.open("gauss_p1.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();
    
// Gaussian psi2
    
    fill(sum.begin(), sum.end(), 0);    // set every element to zero to re-use the same vectors as above
    fill(sum2.begin(), sum2.end(), 0);
    fill(err.begin(), err.end(), 0);

    for (int i = 0; i < N; i++) {
        
        s = 0.;
        s = Metropolis(L, pos, p2, rnd, false, false);
        
        if (i == 0) {
            sum[i] = s;
            sum2[i] = pow(s, 2);
        }
        else {
            sum[i] = (sum[i-1] + s);
            sum2[i] = (sum2[i-1] + pow(s, 2));
        }
    }
    
    // cumulative average, square average and error
    for (int i = 0; i < N; i++) {
        sum[i] = sum[i] / (i+1);
        sum2[i] = sum2[i] / (i+1);
        err[i] = error(sum, sum2, i);
    }
    
    // writing computed data on file
    out.open("gauss_p2.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();

//===================================================================
    
    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}
