// Excercise 09.1:
//      find shortest path amongst 32 cities (on circumference and in square) with simulated annealing genetic algorithm.

#include "func.h"
#include "popolazione.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

    if (argc != 2) {
        cout << "\nUsage: main10.1.exe <seed(int)>\n" << endl;
        return -5;
    }
    
    string shape;
    cout << "Circle(0) or square(1)?" << endl;
    cin >> shape;
    
    // setting params
    int dim = 1;            // how many cromos in the population
    int cities = 32;        // length of each cromo
    int steps = 10000;
    int print = 10;
    
    int n_temp = 14;
    double scale_t = 1.5;
    double T;
    //vector<double> temperatures = {4., 3.5, 3., 2.5, 2., 1.5, 1., 0.75, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05};
    
    ofstream out, t_out;
    rowvec cromo;
    
    // set seed for reproducibility
    int seed = atoi(argv[1]);
    arma_rng::set_seed(seed);

// ==============================================================
    
    string file_pos;
    string file_cost;
    string file_cromo;
    mat pos;
    
    if (shape == "0") {         // circle
        cout << "\n\nCIRCUMFERENCE: searching shortest path..." << endl;
        file_pos = "pos_circle.csv";
        file_cost = "cost_circle.csv";
        file_cromo = "best_cromo_circle.csv";
        pos = make_positions(true, cities);
    } else if (shape == "1") {      // square
        cout << "\n\nSQUARE: searching shortest path..." << endl;
        file_pos = "pos_square.csv";
        file_cost = "cost_square.csv";
        file_cromo = "best_cromo_square.csv";
        pos = make_positions(false, cities);
    } else {
        cout << "\nOnly accepted \"0\" or \"1\" " << endl << endl;
    }
    
    // cities
    out.open(file_pos);
    out << "seed = " << seed << endl;
    for (int i = 0; i < cities; i++)
        out << pos(i,0) << "," << pos(i,1) << endl;
    out.close();
    
    // create population
    Popolazione popu(dim, cities, pos);
    
    // print starting point
    cromo = popu.get_population().row(0);
    out.open(file_cost);
    out << popu.fitness(cromo);
    out.close();
    out.open(file_cromo);
    for (int k = 0; k < cities-1; k++)
        out << cromo(k) << ",";
    out << cromo(cities-1) << endl;
    out.close();
    
    T = 3.*scale_t;     // first temperature is T=1
    t_out.open("temperatures.csv");
    out.open(file_cost, ios::app);
    for (int t = 0; t < n_temp; t++) {
        T = T/scale_t;
        t_out << T << endl;
        cout << "Simulating at temperature " << T << endl;
        for (int j = 0; j < steps; j++) {
            cromo = move(popu, T);
            if (j%print == 0) {
                out << "," << popu.fitness(cromo);
            }
        }
    }
    out.close();
    t_out.close();
    
    cout << "Best path = " << cromo << endl;

    out.open(file_cromo, ios::app);
    for (int k = 0; k < cities-1; k++)
        out << cromo(k) << ",";
    out << cromo(cities-1) << endl;
    out.close();
    
    cout << "\nSUCCESS!" << endl << endl;
    
    return 0;
}
