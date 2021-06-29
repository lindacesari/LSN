// excercise 01.1

#include "random.h"
#include "func.h"

using namespace std;

int main () {
    
    // initializing random generator
    Random rnd;
    set_rnd(rnd);
    
//============================= 01.1.1 ================================
    
    int M = 1000000;     // total
    int N = 100;        // number of blocks
    int L = int(M/N);   // length of each block
    
    vector<double> sum(N, 0.0);     // vectors of progressive results
    vector<double> sum2(N, 0.0);
    vector<double> err(N, 0.0);
    
    double s = 0;       // accumulator
    
    // compute average for each block and store it in the progressive sums vectors
    for (int i = 0; i < N; i++) {
          
        // compute block average
        s = 0;
        for (int j = 0; j < L; j++) {       // accumulate
            s += rnd.Rannyu();
        }
        s = s/double(L);                    // final value
        
        // store progressive results
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
    
    // write computed data on file
    ofstream out("data_r.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << M << "," << N << endl;
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();
    
//============================= 01.1.2 ================================
    
    // set every element to zero to re-use the same vectors
    fill(sum.begin(), sum.end(), 0);
    fill(sum2.begin(), sum2.end(), 0);
    fill(err.begin(), err.end(), 0);
    
    // compute average for each block and store it in the progressive sums vectors
    for (int i = 0; i < N; i++) {
          
        // compute block average
        s = 0;
        for (int j = 0; j < L; j++) {       // accumulate
            s += pow((rnd.Rannyu() - 0.5), 2);
        }
        s = s/double(L);                    // final value
        
        // store progressive results
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
    
    // write computed data on file
    out.open("data_sigma.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << M << "," << N << endl;
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();
    
//============================= 01.1.3 ================================
    
    // initialize variables
    M = 100;
    int n = 10000;
    vector<int> n_i(M, 0);              // count occupation of every interval
    vector<double> chi2(100, 0.0);      // chi square of every bin
    
    for (int h = 0; h < 100; h++) {
        
        fill(n_i.begin(), n_i.end(), 0);
        
        // count how many (of the 10000) x fall in each interval
        for (int j = 0; j < n; j++) {
            n_i[floor(rnd.Rannyu()/(1./double(M)))]++;
        }
        
        // compute chi2
        for (int i = 0; i < M; i++) {
            chi2[h] += pow(n_i[i] - double(n)/double(M), 2);
        }
        chi2[h] = chi2[h] / (double(n)/double(M));
    }
        
    // writing computed data on file
    out.open("data_chi2.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << M << endl << n << endl;
    for (int i = 0; i < N; i++)
        out << chi2[i] << endl;
    out.close();

//===================================================================
    
    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}


