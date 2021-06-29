// excercise 08.1

#include "func.h"

using namespace std;

double get_H ( double mu, double sigma ) ;

int main () {
    
    // initializing random generator
    Random rnd;
    set_rnd(rnd);
    
//============================= 08.1 ================================
    
    Psi psi;
    Integranda f;
    
    int M = 1e6;     // total
    int N = 100;        // number of blocks
    int L = int(M/N);   // length of each block
    
    vector<double> sum(N, 0.0); // initializing vectors of zeros
    vector<double> sum2(N, 0.0);
    vector<double> err(N, 0.0);
    
    // approximate starting values
    double pos = 0.8;
    double sigma = 0.7;
    double mu = 0.9;
    
    // compute average for each block and store it in the progressive sums vectors
    double s, pos_new;
    
    int accepted, attempted;

    for (int i = 0; i < N; i++) {
        
        s = 0.;
        accepted = 0;
        attempted = 0;
        
        for (int j = 0; j < L; j++) {
            pos_new = Metropolis_step (pos, mu, sigma, psi, rnd);
            if (pos_new != pos)
                accepted++;
            attempted++;
            pos = pos_new;
            s = s + f.Eval(pos, mu, sigma);
        }
        s = s/double(L);
        //cout << "Acceptance rate: " << double(accepted)/double(attempted) << endl;
        
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
    
//============================== 08.2 =================================
    
    // find the best parameters
    
    int ntries = 10000;
    
    double H = sum[N-1];    // taking the energy from the starting (mu, sigma)
    
    double new_mu, new_sigma, new_H;
    
    for (int i = 0; i < ntries; i++) {
        
        // taking new values around the previous best one
        new_mu = rnd.Rannyu(mu-0.3, mu+0.3);
        new_sigma = rnd.Rannyu(sigma-0.3, sigma+0.3);
        new_H = get_H(new_mu, new_sigma);
        
        if ( new_H < H ) {
            H = new_H;
            mu = new_mu;
            sigma = new_sigma;
        }
    }
    
    cout << endl << "Best parameters:\tmu = " << mu << "\tsigma = " << sigma << "\tenergy = " << H << endl;
    
    // print the best parameters energy
    
    for (int i = 0; i < N; i++) {
        
        s = 0.;
        
        for (int j = 0; j < L; j++) {
            pos = Metropolis_step (pos, mu, sigma, psi, rnd);
            s = s + f.Eval(pos, mu, sigma);
        }
        s = s/double(L);
        
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
    ofstream out("energy.csv");
    if (!out) {
        cout << "Unable to open output file." << endl;
        exit(-1);
    }
    out << mu << "," << sigma << endl;
    for (int i = 0; i < N; i++)
        out << sum[i] << "," << err[i] << endl;
    out.close();
    
//=====================================================================
    
    // making the histogram
    
    int sampled = 10000;
    
    out.open("histo.csv");
    if (!out) {
        cout << "Unable to open histo file." << endl;
        exit(-1);
    }
    for (int i = 0; i < sampled; i++) {
        pos = Metropolis_step (pos, mu, sigma, psi, rnd);
        out << pos << endl;
    }
    out.close();
    

//======================================================================

    cout << "\nSUCCESS!" << endl << endl;

    return 0;
}


double get_H (double mu, double sigma) {
    
    double pos = 0.;
    Psi psi;
    Integranda f;
    Random rnd;
    set_rnd(rnd);
    
    int M = 1e5;     // total
    int N = 100;        // number of blocks
    int L = int(M/N);   // length of each block
    
    vector<double> sum(N, 0.0); // initializing vectors of zeros
    vector<double> sum2(N, 0.0);
    vector<double> err(N, 0.0);

    for (int i = 0; i < N; i++) {
        
        double s = 0.;
        
        for (int j = 0; j < L; j++) {
            pos = Metropolis_step (pos, mu, sigma, psi, rnd);
            s = s + f.Eval(pos, mu, sigma);
        }
        s = s/double(L);
        
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
    
    // always check that the energy makes sense :)
    if (sum[N-1] < -0.46) {
        cout << "\n\nENERGIA ROTTA\n\n";
        return 9999999;
    }
    
    return sum[N-1];

};
