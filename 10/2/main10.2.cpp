// Exercise 10.2:
//      find shortest path amongst 32 cities placed in a unitary square with genetic parallelized algorithm.

#include "func.h"
#include <mpi.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc , char* argv[]) {

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (argc != 2) {
        cout << "\nUsage: main101.exe <seed(int)>\n" << endl;
        return -5;
    }
    
    cout << "\nProcess n° " << rank << ": \tsimulating on a square..." << endl;
    
    // set seed for reproducibility (comparison with one core)
    int seed;
    if (rank == 0)
        seed = atoi(argv[1]);
    else
        seed = rank;
    arma_rng::set_seed(seed);

    // setting params
    int n_gen = 1000;       // #evolutions
    int dim = 300;          // how many cromos in the population
    int n_migr = 25;        // how often to migrate
    int n_send = 10;        // best cromos to migrate
    int cities = 32;        // length of each cromo
    
    // params for results
    int n_print = 1;        // print config every n_print generations

    // generate positions on square
    mat pos_square = make_positions(false, cities);
    
    // convert from arma::colvec to double[]
    double pos_x[cities];
    double pos_y[cities];
    for (int i = 0; i < cities; i++){
        pos_x[i] = pos_square(i,0);
        pos_y[i] = pos_square(i,1);
    }
    
    // broadcast cities positions from core 0 to all
    if (size > 1) {
        MPI_Bcast(pos_x, cities, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
        MPI_Bcast(pos_y, cities, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
    }

    // make sure every core uses the same city configuration
    for (int i = 0; i < cities; i++){
        pos_square(i,0) = pos_x[i];
        pos_square(i,1) = pos_y[i];
    }

    // print it out
    ofstream out;
    if (rank == 0 && size > 1) {
        out.open("pos_square.csv");
        out << "seed = " << seed << endl;
        for (int i = 0; i < cities; i++)
            out << pos_x[i] << "," << pos_y[i] << endl;
        out.close();
    }
    
    // useful things
    int best_cromos[n_send*cities];
    int migrated_cromos[n_send*cities*size];
    mat config(dim, cities);
    rowvec cromo(cities);
    
    // create the population and sort it
    Popolazione popu(dim, cities, pos_square);
    popu.sort();
    
    // running on one core -> files named with 4
    if (size == 1)
        rank = 4;
    
    // print starting point
    config = popu.get_population();
    out.open("cost"+to_string(rank)+".csv");
    out << popu.fitness(config.row(0)) << ",";
    out.close();
    out.open("best_cromo"+to_string(rank)+".csv");
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();

    // evolve the population and exchange the best cromos
    
    for (int i = 0; i < n_gen; i++) {
        
        popu.sort();
        
        // show evolution of best cromo length
        if (i%n_print == 0 && i != 0) {
            config = popu.get_population();
            out.open("cost"+to_string(rank)+".csv", ios::app);
            out << popu.fitness(config.row(0)) << ",";
            out.close();
        }
        // perform one evolution
        popu.evolve();

        // migrate best cromos every n_migr generations
        if (i%n_migr == 0 && size > 1) {
            
            popu.sort();
            config = popu.get_population();
            
            // prepare vector to send
            for (int j = 0; j < n_send; j++)
                for (int h = 0; h < cities; h++)
                   best_cromos[h+j*cities] = config(j, h);

            // migration
            MPI_Allgather(best_cromos, n_send*cities, MPI_INTEGER, migrated_cromos, n_send*cities, MPI_INTEGER, MPI_COMM_WORLD);
            
            // update population with migrated best cromos:
            //      keep the first dim-(n_send*size) old ones
            //      substitute the last n_send*size rows with the new cromos
            for (int i = 0; i < n_send*size; i++) {
                for (int k = 0; k < cities; k++)
                    cromo(k) = migrated_cromos[k+i*cities];
                config.row(dim-i-1) = cromo;
            }
            popu.set_population(config);
        }
    }
    
    // final sorting to extract the best cromo
    popu.sort();
    
    // print best cromo length
    config = popu.get_population();
    out.open("cost"+to_string(rank)+".csv", ios::app);
    out << popu.fitness(config.row(0)) << ",";
    out.close();
    // print best cromo
    out.open("best_cromo"+to_string(rank)+".csv", ios::app);
    for (int k = 0; k < cities-1; k++)
        out << config(0,k) << ",";
    out << config(0,cities-1) << endl;
    out.close();
    
    cout << "SUCCESS!\n";
    
    MPI_Finalize();

    return 0;
    
}
