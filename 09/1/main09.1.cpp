// Excercise 09.1:
//      find shortest path amongst 32 cities (on circumference and in square) with genetic algorithm.

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
        cout << "\nUsage: main09.1.exe <seed(int)>\n" << endl;
        return -5;
    }
    
    string shape;
    cout << "Circle(0) or square(1)?" << endl;
    cin >> shape;
    
    // setting params
    int n_gen = 500;        // #evolutions
    int dim = 300;          // how many cromos in the population
    int cities = 32;        // length of each cromo
    
    ofstream out;
    double L = 0.;
    
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
    mat config(dim, cities);

    // print starting point
    config = popu.get_population();
    out.open(file_cromo);
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();
    
    out.open(file_cost);
    out << popu.fitness(config.row(0)) << ",";
    for (int k = 0; k < dim/2.; k++)
        L += popu.fitness(config.row(k));
    L /= int(dim/2.);
    out << L << endl;
    out.close();
    
    out.open(file_cost, ios::app);
    for (int i = 0; i < n_gen; i++) {
        
        // evolve
        popu.evolve();
        
        // save current best path length
        config = popu.get_population();
        out << popu.fitness(config.row(0)) << ",";
        
        // best half of the population <L>
        L = 0.;
        for (int k = 0; k < dim/2.; k++)
            L += popu.fitness(config.row(k));
        L /= int(dim/2.);
        out << L << endl;
        
    }
    out.close();
    
    cout << "Best path = " << config.row(0) << endl;
    
    // print ending point
    config = popu.get_population();
    out.open(file_cromo, ios::app);
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();
    
// ============================ CIRC ============================
/*
    cout << "\n\nCIRCUMFERENCE: searching shortest path..." << endl;
    
    // cities on circumference r = 1
    mat pos_circle = make_positions(true, cities);
    out.open("pos_circle.csv");
    out << "seed = " << seed << endl;
    for (int i = 0; i < cities; i++)
        out << pos_circle(i,0) << "," << pos_circle(i,1) << endl;
    out.close();
    
    // create population
    Popolazione popu_circ(dim, cities, pos_circle);
    
    mat config(dim, cities);

    // print starting point
    config = popu_circ.get_population();
    out.open("cost_circle.csv");
    out << popu_circ.fitness(config.row(0));
    out.close();
    double L = 0.;
    for (int k = 0; k < dim/2.; k++)
        L += popu_circ.fitness(config.row(k));
    L /= int(dim/2.);
    out.open("ave_cost_circle.csv");
    out << L;
    out.close();
    out.open("best_cromo_circle.csv");
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();

    for (int i = 0; i < n_gen; i++) {
        
        // evolve
        popu_circ.evolve();
        
        // save current best path length
        config = popu_circ.get_population();
        out.open("cost_circle.csv", ios::app);
        out << "," << popu_circ.fitness(config.row(0));
        out.close();
        
        // best half of the population <L>
        L = 0.;
        for (int k = 0; k < dim/2.; k++)
            L += popu_circ.fitness(config.row(k));
        L /= int(dim/2.);
        out.open("ave_cost_circle.csv", ios::app);
        out << "," << L;
        out.close();
        
    }
    
    cout << "Best path = " << config.row(0) << endl;
    
    // print ending point
    config = popu_circ.get_population();
    out.open("best_cromo_circle.csv", ios::app);
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();
    
// ============================ SQUARE ============================
        
    // set seed for reproducibility
    arma_rng::set_seed(atoi(argv[1]));
    
    cout << "\n\nSQUARE: searching shortest path..." << endl;
        
    // cities in square [-1, 1]x[-1, 1]
    mat pos_square = make_positions(false, cities);
    out.open("pos_square.csv");
    out << "seed = " << seed << endl;
    for (int i = 0; i < cities; i++)
        out << pos_square(i,0) << "," << pos_square(i,1) << endl;
    out.close();
        
    // create population
    Popolazione popu_square(dim, cities, pos_square);

    // print starting point
    config = popu_square.get_population();
    out.open("cost_square.csv");
    out << popu_square.fitness(config.row(0));
    out.close();
    L = 0.;
    for (int k = 0; k < dim/2.; k++)
        L += popu_square.fitness(config.row(k));
    L /= int(dim/2.);
    out.open("ave_cost_square.csv");
    out << L;
    out.close();
    out.open("best_cromo_square.csv");
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();

    for (int i = 0; i < n_gen; i++) {
        
        // evolve
        popu_square.evolve();
        
        // save current best path length
        config = popu_square.get_population();
        out.open("cost_square.csv", ios::app);
        out << "," << popu_square.fitness(config.row(0));
        out.close();
        
        // best half of the population <L>
        L = 0.;
        for (int k = 0; k < dim/2.; k++)
            L += popu_square.fitness(config.row(k));
        L /= int(dim/2.);
        out.open("ave_cost_square.csv", ios::app);
        out << "," << L;
        out.close();
        
    }
    cout << "Best path = " << config.row(0) << endl;
    
    // print ending point
    config = popu_square.get_population();
    out.open("best_cromo_square.csv", ios::app);
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();
  */
// =========================================================
        
    cout << "\nSUCCESS!" << endl << endl;
    
    return 0;
}
