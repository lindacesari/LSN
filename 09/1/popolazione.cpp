#include "popolazione.h"

using namespace std;
using namespace arma;

Popolazione::Popolazione() {
    
    m_dimension = 100;
    m_cities = 32;
    m_config.set_size(m_dimension, m_cities);
    rowvec cromo(m_cities);
    
    for (int i = 0; i < m_dimension; i++) {
        cromo = regspace(1, m_cities);
        cromo.subvec(1, m_cities-1) = shuffle(cromo.subvec(1, m_cities-1));
        m_config.row(i) = cromo;
    }
      
    check();
    
    // if no positions were provided, fill m_positions with zeros
    m_positions = zeros<mat>(m_cities, 2);
}

Popolazione::Popolazione (int dimension, int cities, mat positions) {
    
    m_dimension = dimension ;
    m_cities = cities ;
    m_config.set_size(m_dimension, m_cities);
    rowvec cromo;
    
    for (int i = 0; i < m_dimension; i++) {
        cromo = regspace<rowvec>(1, 1, m_cities);
        cromo.subvec(1, m_cities-1) = shuffle(cromo.subvec(1, m_cities-1));
        m_config.row(i) = cromo;
    }
    
    check();
    
    m_positions = positions;
}

Popolazione::~Popolazione() {}

bool Popolazione::check() const {

    bool ok = true;
    colvec a(m_dimension, fill::ones);
    
    if ( arma::size( find(m_config.col(0) == a) ) != arma::size(a) ) {
        ok = false;
        exit(-1);
    }
        
    for (int i = 0; i < m_dimension; i++)
        if (size(unique(m_config.row(i))) != size(m_config.row(i))) {
            ok = false;
            exit(-2);
        }
    
    return ok;
}

bool Popolazione::check(rowvec cromo) const {

    bool ok = true;
    colvec a(m_dimension, fill::ones);
    
    if ( cromo(0) != 1 ) {
        ok = false;
        exit(-1);
    }
    
    if (size(unique(cromo)) != size(cromo)) {
        ok = false;
        exit(-2);
    }
    
    return ok;
}


rowvec Popolazione::mutate(rowvec offspring) {
    
    int index, m, r1, r2, n;
    double prob = 0.5;
    
    // mutation 1
    //      swap two random cities
    if (randu() < prob){
        r1 = randi<double>( distr_param(1, m_cities - 1) );
        r2 = randi<double>( distr_param(1, m_cities - 1) );
        offspring.swap_cols(r1, r2);
    }
    
    // mutation 1.1
    //      swap two random adjacent cities
    if (randu() < prob){
        r1 = randi<double>( distr_param(1, m_cities - 2) );
        double f = fitness(offspring);
        offspring.swap_cols(r1, r1+1);
        if (fitness(offspring) > f)
            offspring.swap_cols(r1, r1+1);
    }
           
    // mutation 2
    //      shift part of the offspring
    if (randu() < prob) {
        index = int(randi<double>(distr_param(1, m_cities-2)));
        n = int(randi<double>(distr_param(1, m_cities-index-1)));     // length of shift push
        offspring.subvec(index, m_cities-1) = shift(offspring.subvec(index, m_cities-1), n);
    }
    
    // mutation 3
    //      permutation among m contiguous cities with other m contiguous cities
    if (randu() < prob) {
        m = int(randi<double>(distr_param( 1, int((m_cities-1)/2.) )));
        r1 = int(randi<double>(distr_param( 1, m_cities-2 )));
        if ( ((r1 + m + m_cities)%m_cities) < ((r1 - m + 2*m_cities - 1)%m_cities) )
            r2 = randi<double>( distr_param( (r1 + m + m_cities)%m_cities, (r1 - m + 2*m_cities - 1)%m_cities ));
        else r2 = randi<double>( distr_param( (r1 - m + 2*m_cities - 1)%m_cities, (r1 + m + m_cities)%m_cities ));
        
        rowvec helper = offspring.tail(m_cities-1); // helper vector without 1
        for (int i = 0; i < m; i++)             // swap col by col the two subvectors starting at r1 and r2, respectively
            helper.swap_cols( (i+r1)%(m_cities-1), (r2+i)%(m_cities-1));
        offspring = join_rows(rowvec{1}, helper);        // mutation permformed
    }
    
    // mutation 4
    //      invert m contiguous cities
    if (randu() < prob) {
        index = int(randi<double>(distr_param( 1, m_cities-2 )));
        m = int(randi<double>(distr_param( 1, m_cities-index )));
        offspring.subvec(index, index+m-1) = reverse(offspring.subvec(index, index+m-1));
    }
    
    check(offspring);
    return offspring;
    
}

mat Popolazione::crossover(mat parents) {
    
    // read the two parents (they are rows in the matrix parents)
    rowvec mother = parents.row(0);
    rowvec father = parents.row(1);
    
    // random position for the cut
    int index = int(randi<double>(distr_param( 1, m_cities-1 )));
    //cout << "Index: \t" << index <<endl;
    
    // cut vectors
    rowvec head_mother = mother.head(index+1);
    rowvec head_father = father.head(index+1);
    rowvec tail_mother = mother.tail(m_cities-1-index);
    rowvec tail_father = father.tail(m_cities-1-index);
        
    int helper_index = 0;
    rowvec temp;

    // offspring 1: mother's head
    rowvec tail_offspring1(m_cities-index-1);
    for (int i = 0; i < m_cities; i++) {
        temp = intersect(rowvec{father(i)}, tail_mother);
        if (temp.n_cols > 0) {
            tail_offspring1(helper_index) = temp(0);
            helper_index++;
        }
    }
    
    // offspring 2: father's head
    helper_index = 0;
    rowvec tail_offspring2(m_cities-index-1);
    for (int i = 0; i < m_cities; i++) {
        temp = intersect(rowvec{mother(i)}, tail_father);
        if (temp.n_cols > 0) {
            tail_offspring2(helper_index) = temp(0);
            helper_index++;
        }
    }

    // join heads and tails for the final offsprings
    rowvec offspring1 = join_rows(head_mother, tail_offspring1);
    rowvec offspring2 = join_rows(head_father, tail_offspring2);
    
    // matrix to be returned (each row is an offspring)
    check(offspring1);
    check(offspring2);
    mat offsprings = join_cols(offspring1, offspring2);
    
    //offsprings.print("offsprings");
    
    return offsprings;
}

mat Popolazione::select() {
    
    mat parents(2, m_cities);
    
    double p = 4.;       // p large <-> small index <-> best fitness cromo are favoured
    double r = randu();
    double s = randu();
    while ( r == 1 )        // 1 should not be included
        r = randu();
    while( s == 1 )
        s = randu();
    
    int r1 = int( m_dimension*pow(r,p) ) + 1;
    if (r1 == m_dimension) r1--;
    int r2 = int( m_dimension*pow(s,p) ) + 1;
    if (r2 == m_dimension) r2--;
    
    parents.row(0) = m_config.row(r1);
    parents.row(1) = m_config.row(r2);
    
    return parents;
}

double Popolazione::fitness(rowvec cromo) const {
    
    // cost function
    double L = 0.;
    
    rowvec pos1, pos2;
    int city;
    
    // coordinates are cartesian
    pos1 = m_positions.row(0);
    
    for (int i = 1; i < m_cities; i++) {
        city = cromo(i);
        pos2 = m_positions.row(city-1);
        L += sqrt( pow(pos2(0)-pos1(0), 2) + pow(pos2(1)-pos1(1), 2) );
        pos1 = pos2;
    }
    
    // distance from last city to first
    pos2 = m_positions.row(0);
    L += sqrt( pow(pos2(0)-pos1(0), 2) + pow(pos2(1)-pos1(1), 2) );
    
    return L;
}


void Popolazione::sort() {
    
    vec distances(m_dimension);
    
    for (int i = 0; i < m_dimension; i++) {
        distances(i) = fitness(m_config.row(i));
    }
    
    // sorting in ascending order
    uvec indices = sort_index(distances, "ascend");
    
    mat copy = m_config;
    for (int i = 0; i < m_dimension; i++) {
        m_config.row(i) = copy.row(indices(i));
    }
    
}

void Popolazione::evolve() {
    
    int keep = int(m_dimension/2.);
    mat new_generation(m_dimension-keep, m_cities);
    
    // sort by goodness of fitness
    sort();
    
    // one new generation changes m_dimension/2 elements
    for (int i = 0; i < m_dimension-keep; i += 2) {
        
        // select parents, preferably the first ones (the best ones)
        mat parents = select();
        
        // crossover with given probability
        mat offsprings(2, m_cities);
        if (randu() < 0.6)
            offsprings = crossover(parents);
        else
            offsprings = parents;
            
        // mutate
        rowvec offspring1 = mutate(offsprings.row(0));
        rowvec offspring2 = mutate(offsprings.row(1));
        
        // store new cromosomes
        new_generation.row(i) = offspring1;
        new_generation.row(i+1) = offspring2;
    }
    
    // update the population
    m_config.rows( keep, m_dimension-1 ) = new_generation.rows( 0, m_dimension-keep-1 );
    
    assert(size(m_config.col(0)) == size(colvec(m_dimension)));
    check();
    
    // final sorting
    sort();
}
