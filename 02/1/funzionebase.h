#ifndef __funzionebase_h__
#define __funzionebase_h__

#include <vector>

using namespace std;

class FunzioneBase {

    public:
    virtual double Eval (double x) const =0;

};

// initial integrand
class Integranda : public FunzioneBase {

    public:
    Integranda() {};
    ~Integranda() {};
    double Eval(const double x) const override { return ( M_PI/2. * cos(M_PI*x/2.) ); };

};

// integrand, normalized to the importance sampling distribution
class Integranda_IS : public FunzioneBase {
    
    public:
    Integranda_IS() {};
    ~Integranda_IS() {};
    double Eval(const double x) const override { return ( M_PI/3. * cos(M_PI*x/2.) / (1-pow(x,2)) ); };

};

// probability density for the importance sampling
class Densita : public FunzioneBase {
    
    public:
    Densita() {};
    ~Densita() {};
    double Eval(const double x) const override { return ( 3./2. * (1-pow(x,2)) ); };

};

#endif
