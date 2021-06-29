/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <sstream>
#include <string>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(int argc, char**argv) {
    
    if (argc != 2) {
        cout << "\nUsage: MolDyn_NVE.cpp <state(solid, liquid, gas)>\n\n" ;
    }
    string state = argv[1];
    
    //Inizialization
    Input(state);
    nblocks = 20;
    N = (nstep/10)/nblocks;
    
    // move and measure
    for (int istep = 1; istep <= nstep; ++istep) {
        
        Move();           //Move particles with Verlet algorithm
        if (istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if (istep%10 == 0) {
            Measure(state);     //Properties measurement
        }
    }
    
    // block averages reading from files
    cout << endl << "Computing block averages..." << endl;

    // to write
    ifstream inT, inEK, inEP, inET, inP;
    inT.open(state+"/output_temp.dat");
    inEK.open(state+"/output_ekin.dat");
    inEP.open(state+"/output_epot.dat");
    inET.open(state+"/output_etot.dat");
    inP.open(state+"/output_press.dat");
    
    // to read
    ofstream outEP, outEK, outET, outT, outP;
    outT.open(state+"/ave_temp.dat", ios::app);
    outEK.open(state+"/ave_ekin.dat", ios::app);
    outEP.open(state+"/ave_epot.dat", ios::app);
    outET.open(state+"/ave_etot.dat", ios::app);
    outP.open(state+"/ave_press.dat", ios::app);

    data_blocking(inT, outT);
    data_blocking(inEK, outEK);
    data_blocking(inEP, outEP);
    data_blocking(inET, outET);
    data_blocking(inP, outP);
    
    inT.close();
    inEK.close();
    inEP.close();
    inET.close();
    inP.close();
    outT.close();
    outEK.close();
    outEP.close();
    outET.close();
    outP.close();
    
    // g(r)
    ofstream outG;
    ifstream inG;

    outG.open(state+"/ave_gofr.dat");
    inG.open(state+"/output_gofr.dat");

    data_blocking_bid(inG, outG);
    
    // write final configuration to restart
    ConfFinal();

    return 0;
}


//Prepare all stuff for the simulation
void Input(string state) {
    
    ifstream ReadInput,ReadConf,ReadPreConf;
    string equil, start;

    cout << "First run? (y/n)" << endl;
    cin >> start;
    cout << "Equilibrate? (y/n)" << endl;
    cin >> equil;
    
    if (equil != "y" && equil != "n") {
        cout << endl << "Only accepted \"y\" or \"n\" " << endl << endl;
        exit(-1);
    }
    if (start != "y" && start != "n") {
        cout << endl << "Only accepted \"y\" or \"n\" " << endl << endl;
        exit(-1);
    }

    cout << "\n\nClassic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator
  
    string input_file = "input."+state;
    ReadInput.open(input_file); //Read input

    ReadInput >> temp;
    cout << "Temperature = " << temp << endl;
    
    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();

    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables

    //Measure of g(r)
    igofr = 4;
    n_props = n_props + nbins;
    bin_size = (box/2.0)/(double)nbins;

    //Read initial configuration
    if (start == "y") {
        cout << "Read initial configuration from file config.0" << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
    } else {
        cout << "Read initial configuration from file config.final" << endl << endl;
        ReadConf.open("config.final");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
    }
    
    //Prepare initial velocities
    if ( start == "y" & equil =="y" ){
       cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
       double sumv[3] = {0.0, 0.0, 0.0};
       for (int i=0; i<npart; ++i){
         vx[i] = rand()/double(RAND_MAX) - 0.5;
         vy[i] = rand()/double(RAND_MAX) - 0.5;
         vz[i] = rand()/double(RAND_MAX) - 0.5;

         sumv[0] += vx[i];
         sumv[1] += vy[i];
         sumv[2] += vz[i];
       }

       for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
       double sumv2 = 0.0, fs;
       for (int i=0; i<npart; ++i){
         vx[i] = vx[i] - sumv[0];
         vy[i] = vy[i] - sumv[1];
         vz[i] = vz[i] - sumv[2];

         sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
       }
       sumv2 /= (double)npart;

       fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
       for (int i=0; i<npart; ++i){
         vx[i] *= fs;
         vy[i] *= fs;
         vz[i] *= fs;

         xold[i] = Pbc(x[i] - vx[i] * delta);
         yold[i] = Pbc(y[i] - vy[i] * delta);
         zold[i] = Pbc(z[i] - vz[i] * delta);
       }
    }
    else if (start != "y" & equil== "y") {
        double sumv[3] = {0.0, 0.0, 0.0};

        ReadPreConf.open("config.almostfinal");

        for (int i=0; i<npart; i++){
            ReadPreConf >> xold[i] >> yold[i] >> zold[i] ;
            xold[i] *= box ;
            yold[i] *= box ;
            zold[i] *= box ;
        }

        ReadPreConf.close();

        for (int i=0; i<npart; i++){
            vx[i] = Pbc(x[i] - xold[i])/(2.0 * delta);
            vy[i] = Pbc(y[i]- yold[i])/(2.0 * delta);
            vz[i] = Pbc(z[i] - zold[i])/(2.0 * delta);

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }

        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        double sumv2 = 0., fs ;
        for (int i=0; i<npart; ++i){
          vx[i] = vx[i] - sumv[0];
          vy[i] = vy[i] - sumv[1];
          vz[i] = vz[i] - sumv[2];

          sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
          vx[i] *= fs;
          vy[i] *= fs;
          vz[i] *= fs;

          xold[i] = Pbc(x[i] - vx[i] * delta);
          yold[i] = Pbc(y[i] - vy[i] * delta);
          zold[i] = Pbc(z[i] - vz[i] * delta);
      }
    }

    else if ( equil != "y" ){
        ReadPreConf.open("config.almostfinal");

        for (int i=0; i<npart; i++){
            ReadPreConf >> xold[i] >> yold[i] >> zold[i] ;
            xold[i] *= box ;
            yold[i] *= box ;
            zold[i] *= box ;
        }

        ReadPreConf.close();
    }

    return;
}

//Move particles with Verlet algorithm
void Move(void) {
    
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

    for (int i=0; i<npart; ++i) { //Force acting on particle i
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
    }

    for (int i=0; i<npart; ++i) { //Verlet integration scheme

        xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

        vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
        vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
        vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

//Compute forces as -Grad_ip V(r)
double Force(int ip, int idir) {
    
    double f=0.0;
    double dvec[3], dr;

    for (int i=0; i<npart; ++i){
        if(i != ip){
            dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
            dvec[1] = Pbc( y[ip] - y[i] );
            dvec[2] = Pbc( z[ip] - z[i] );

            dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
            dr = sqrt(dr);
            
            if(dr < rcut){
            f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
            }
        }
    }
  
  return f;
}

//Properties measurement
void Measure(string state) {
    
    int bin;
    double v, w, t, vij, wij;
    double dx, dy, dz, dr, dv;
    ofstream Epot, Ekin, Etot, Temp, Press, Gofr;

    Temp.open(state+"/output_temp.dat",ios::app);
    Press.open(state+"/output_press.dat",ios::app);
    Ekin.open(state+"/output_ekin.dat",ios::app);
    Epot.open(state+"/output_epot.dat",ios::app);
    Etot.open(state+"/output_etot.dat",ios::app);
    Gofr.open(state+"/output_gofr.dat",ios::app);

    v = 0.0; //reset observables
    t = 0.0;
    w = 0.0;

    //reset the hystogram of g(r)
    for (int k=0; k < nbins ; ++k)
        walker[k]=0.0;

    //cycle over pairs of particles
    for (int i=0; i<npart-1; ++i) {
        for (int j=i+1; j<npart; ++j) {

            dx = Pbc(x[i] - x[j]);
            dy = Pbc(y[i] - y[j]);
            dz = Pbc(z[i] - z[j]);

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);

            //update of the histogram of g(r) (cycling over pairs -> 2 increments)
            bin = int(dr/bin_size);
            if (bin < 100)
                walker[bin] += 2.;
            
            if(dr < rcut){
                vij = 4.0/pow(dr,12.) - 4.0/pow(dr,6.);
                wij = 1.0/pow(dr,12.) - 0.5/pow(dr,6.);
                //Potential energy
                v += vij;
                w += wij;
            }
        }
    }
    w *= 48./3.;

    //Kinetic energy
    for (int i=0; i<npart; ++i)
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
    // g(r)
    for(int i=0; i < nbins ; i++){
          dv = 4.*M_PI*( pow(bin_size*(double)(i+1.), 3. ) - pow( bin_size*(double)(i), 3.))/3.;
          walker[i] /= rho * dv * (double)npart;
      }
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_pres = rho*temp + w/vol;    // Pressure
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << setprecision(10) << stima_pot  << endl;
    Ekin << setprecision(10) << stima_kin  << endl;
    Temp << setprecision(10) << stima_temp << endl;
    Press << setprecision(10) << stima_pres << endl;
    Etot << setprecision(10) << stima_etot << endl;

    Temp.close();
    Epot.close();
    Ekin.close();
    Press.close();
    Etot.close();
    
    // g(r)
    for (int i = 0 ; i < nbins ; i++ ){
        Gofr << walker[i];
        if (i != nbins-1) Gofr << ",";
        else Gofr << endl;
    }
    Gofr.close();

    return;
}

//Write final configuration
void ConfFinal(void) {
    
    ofstream WriteConf, WritePreConf;

    cout << "Print final configuration to file config.final" << endl << endl;
    
    WriteConf.open("config.final");
    WritePreConf.open("config.almostfinal");

    for (int i=0; i<npart; ++i) {
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
        WritePreConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    
    WriteConf.close();
    WritePreConf.close();
    
    return;
}

//Write configuration in .xyz format
void ConfXYZ(int nconf) {
    
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

//Algorithm for periodic boundary conditions with side L=box
double Pbc(double r){
    return r - box * rint(r/box);
}

// data blocking
void data_blocking(ifstream& in, ofstream& out) {
    
    // vectors of progressive results
    vector<double> sum(nblocks, 0.0);
    vector<double> sum2(nblocks, 0.0);
    vector<double> err(nblocks, 0.0);

    double input, s;

    for (int i = 0; i < nblocks; i++) {
        
        input = 0.;
        s = 0.;
        
        // accumulate
        for (int j = 0; j < N; j++) {
            in >> input;
            s += input;
        }
        s /= double(N);

        // store progressive results
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
    for (int i = 0; i < nblocks; i++) {
        sum[i] = sum[i] / (i+1);
        sum2[i] = sum2[i] / (i+1);
        err[i] = error(sum, sum2, i);
    }


    for (int i = 0; i < nblocks; i++)
        out << setprecision(10) << sum[i] << "," << setprecision(10) << err[i] << endl;
    
    return;
}

void data_blocking_bid(ifstream& in, ofstream& out) {
    
    //Data blocking for g(r)
    vector<vector<double>> sum(nbins, vector<double>(nblocks, 0.));
    vector<vector<double>> sum2(nbins, vector<double>(nblocks, 0.));
    vector<vector<double>> err(nbins, vector<double>(nblocks, 0.));

    vector<double> input;
    vector<double> s(nbins, 0.);
    string aux;
    
    for (int i = 0; i < nblocks; i++) {
        
        // clean
        fill(s.begin(), s.end(), 0.);
        
        for (int j = 0; j < N; j++) {
            
            // get the string and turn it into vector<double>
            input.clear();
            getline(in, aux);
            stringstream ss(aux);
            while(ss.good()) {
                string substr;
                getline(ss, substr, ',');
                input.push_back(stod(substr));
            }
            for (int h = 0; h < nbins; h++)
                s[h] += input[h];
        }
        for (int h = 0; h < nbins; h++)
            s[h] /= double(N);

        // average
        if (i==0) {
          for (int k = 0; k < nbins; k++) {
              sum[k][i] = s[k];
              sum2[k][i] = s[k]*s[k];
          }
      } else {
          for (int k = 0; k < nbins; k++) {
              sum[k][i] = sum[k][i-1] + s[k];
              sum2[k][i] = sum2[k][i-1] + s[k]*s[k];
          }
        }

    }


  for (int i = 0; i < nblocks; i++) {
      for (int k = 0; k < nbins; k++) {
          sum[k][i] /= double(i+1);
          sum2[k][i] /= double(i+1);
      }
  }

  for (int k = 0; k < nbins; k++) {
      for (int i = 0; i < nblocks; i++) {
          if (i==0) err[k][i] = 0.;
          else {
              err[k][i] = error(sum[k], sum2[k], i);
          }
      }
  }

  for (int k = 0; k < nbins; k++)
    out << sum[k][nblocks-1] << " , " << err[k][nblocks-1] << endl;

}

// error for datablocking
double error(vector<double> ave, vector<double> ave2, int n) {
    if (n == 0) {
        return 0.;
    }
    else {
        return sqrt( (ave2[n] - pow(ave[n], 2)) / n );
    }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
