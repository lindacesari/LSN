#ifndef __funzionebase_h__
#define __funzionebase_h__

#include <vector>
#include <cmath>

using namespace std;

class FunzioneBase {

    public:
    virtual double Eval (vector<double> x) const =0;

};

class probabilita1 : public FunzioneBase {

    public:
    probabilita1() {};
    ~probabilita1() {};
    double Eval(const vector<double> x) const override { return 1./M_PI * exp(-2*sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])); };

};

class probabilita2 : public FunzioneBase {
    
    public:
    probabilita2() {};
    ~probabilita2() {};
    double Eval(const vector<double> x) const override { return 1./(32*M_PI) * exp(-sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])) * x[2]*x[2]; };

};

#endif
