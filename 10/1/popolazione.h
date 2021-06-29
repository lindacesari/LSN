// genetic algorithm

#ifndef __Popolazione_h__
#define __Popolazione_h__

#include <string>
#include <armadillo>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace arma;

class Popolazione {

  public:

    // constructor
    Popolazione();
    Popolazione(int dimension, int cities, mat positions);
    
    // destructor
    ~Popolazione();

    // methods
    double fitness(rowvec cromo) const;
    void sort();
    mat select();
    mat crossover(mat);
    rowvec mutate(rowvec offspring);
    void evolve();
    bool check() const;
    bool check(rowvec) const;
    
    int get_dim() const { return m_dimension; };
    int get_cities() const { return m_cities; };
    mat get_population() const { return m_config; };
    mat get_positions() const { return m_positions; };
    void set_population(mat p) { m_config = p; };

  protected:

    mat m_config;       // current population
    mat m_positions;    // cities layout
    int m_dimension;    // number of cromos (size of population)
    int m_cities;       // length of cromo
    
};


#endif // __Popolazione_h__

