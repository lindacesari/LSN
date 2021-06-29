/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
//parameters, observables
const int m_props=4;
int n_props;
const int nbins = 100;
double walker[ nbins ];
int iv,ik,it,ie,igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
int nblocks , N;
double bin_size;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(std::string);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(std::string);
double Force(int, int);
double Pbc(double);
void data_blocking(std::ifstream&, std::ofstream&);
void data_blocking_bid(std::ifstream&, std::ofstream&);
double error(std::vector<double>, std::vector<double>, int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
