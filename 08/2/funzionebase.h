#ifndef __funzionebase_h__
#define __funzionebase_h__

#include <vector>
#include <cmath>

using namespace std;

class FunzioneBase {

    public:
    virtual double Eval (const double x, double mu, double sigma) const = 0;

};

class Psi : public FunzioneBase {

    public:
    Psi() {};
    ~Psi() {};
    double Eval(const double x, double mu, double sigma) const override { return exp( -pow((x-mu), 2.)/(2*sigma*sigma) ) + exp( -pow((x+mu), 2.)/(2*sigma*sigma) ); };

};

class Integranda : public FunzioneBase {
    
    public:
    Integranda() {};
    ~Integranda() {};
    double Eval(const double x, double mu, double sigma) const override {
        
        double psi = exp( -pow((x-mu), 2.)/(2.*sigma*sigma) ) + exp( -pow((x+mu), 2.)/(2.*sigma*sigma) );
        double kin = -( -psi/(sigma*sigma) + pow(x-mu, 2.) * exp(-1.*pow((x-mu), 2.)/(2*pow(sigma, 2.))) / pow(sigma, 4.) + pow(x+mu, 2.)* exp(-1.*pow((x+mu), 2.)/(2*pow(sigma, 2.))) / pow(sigma, 4.) )/2.;
        double pot = psi*pow(x, 4.) - 5./2.*x*x*psi;
        
        return (pot+kin)/psi;
        
    };

};


#endif
