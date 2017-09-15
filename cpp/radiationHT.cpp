#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <thread>
#include <vector>
#include <cmath>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <boost/filesystem.hpp>

#include "radiationHT.h"

const int width = 15; // output width
namespace fs = boost::filesystem;

RadiationHT::RadiationHT(int flavor, std::string fileName, bool plain, bool refresh)
:   HQflavor_(flavor),
    Ntime_(2), Ntemp_(11), NE_(10),
    timeL_(2), timeH_(4), tempL_(0.15), tempH_(0.45), EL_(10), EH_(100),
    dtime_((timeH_-timeL_)/(Ntime_-1.)), dtemp_((tempH_ - tempL_)/(Ntemp_-1.)), dE_((EH_-EL_)/(NE_-1.)),
    RadTab(boost::extents[Ntime_][Ntemp_][NE_]), Max_dNg(boost::extents[Ntime_][Ntemp_][NE_])
{
    if (flavor == 4)    HQmass_ = cMass;
    else if (flavor == 5)   HQmass_ = bMass;
    else
    {
        std::cerr << "Flavor " << flavor << " is currently not implemented!" << std::endl;
        exit (EXIT_FAILURE);
    }

    bool fileExist = fs::exists(fileName);
    if (!fileExist || (fileExist && refresh)) 
    {
        std::cout << "Generating radiation table by HT formula" << std::endl;
        std::cout << "threads " << std::thread::hardware_concurrency() << std::endl;
        std::vector<std::thread> threads;

        size_t Ncores = std::min(size_t(std::thread::hardware_concurrency()), NE_);
        size_t call_per_core = size_t(NE_ * 1./Ncores);
        size_t call_last_core = NE_ - call_per_core * (Ncores -1);

        auto code = [this](size_t Nstart_, size_t dN_) { this->tabulate(Nstart_, dN_);};

        for (size_t i=0; i<Ncores; i++)
        {
            size_t Nstart = i*call_per_core;
            size_t dN = (i==Ncores-1) ? call_last_core : call_per_core;
            threads.push_back(std::thread(code, Nstart, dN));
        }
        
        for (std::thread& t: threads) t.join();

        save_to_file(fileName, plain);
    } else {
        //std::cout << "Reading in radiation table:" << std::endl;
        read_from_file(fileName, plain);
    }
}



void RadiationHT::save_to_file(std::string fileName, bool plain)
// save radiation table, check if plain .dat or h5 
{
    if (plain) {
        fs::ofstream fRad{fileName};
        fRad.precision(6);
        fRad << "#   time     temp     HQenery     Ng     Max_dNgdxdydt" << std::endl;
        double time, temp, energy;
        for (size_t i=0; i<Ntime_; i++){
          time = timeL_ + i*dtime_;
          for (size_t j=0; j<Ntemp_; j++){
            temp = tempL_ + j*dtemp_;
            for (size_t k=0; k<NE_; k++) {
              energy = EL_ + k*dE_;
              fRad << std::setw(8) << time 
                   << std::setw(8) << temp
                   << std::setw(8) << energy
                   << std::setw(width) << RadTab[i][j][k]
                   << std::setw(width) << Max_dNg[i][j][k] << std::endl;
            }
          }
        }
        std::cout << "Writing to file successfully :)" << std::endl;
    }
}



void RadiationHT::read_from_file(std::string fileName, bool plain)
// read from radiation table, check if .dat or .h5
{
    if (plain) {
        fs::ifstream fRad(fileName);
        double time, temp, energy;
        std::string header;
        std::getline(fRad, header);
        std::cout << header << std::endl;
        for(size_t i=0; i<Ntime_; i++) 
          for (size_t j=0; j<Ntemp_; j++)
            for (size_t k=0; k<NE_; k++){ 
              fRad >> time >> temp >> energy >> RadTab[i][j][k] >> Max_dNg[i][j][k];
            }
        //std::cout << "Read in radiative table successfully :)" << std::endl;

    }
}


double RadiationHT::AlphaS(double kT, double temp)
// strong coupling constant, kT is the gluon transverse momentum
{
    if (kT < M_PI*temp)    kT = M_PI * temp;

    int nflavor;
    double lambdas;
    if (kT<cMass) {
        nflavor = 3;
        lambdas = 0.2;
    } else if (kT<bMass) {
        nflavor = 4;
        lambdas = 0.172508;
    } else {
        nflavor = 5;
        lambdas = 0.130719;
    }

    double result=4.*M_PI/((11.-2.*nflavor/3.)*(2.*log(kT/lambdas)));
    return result;
}

double RadiationHT::SplittingP(double x) 
{
    return (2.-2.*x + x*x)*CF/x;
}



double RadiationHT::TauF(double x, double y, double HQenergy)
// gluon formation time
{
    double k0_gluon = x*HQenergy;
    double kT_gluon = y*k0_gluon;
    return 2.*k0_gluon * (1.-x)/(pow(kT_gluon, 2) + pow(x*HQmass_,2));
}



double RadiationHT::dNg_over_dxdydt(double time, double temp, double HQenergy, double x, double y, double &maxdNg)
// HT formula (as a function of seperated time, temperature, HQenergy, x, y), 
{
    double qhat = 4.*ALPHA*pow(temp,3)/CF;

    double scale = M_PI*temp;
    double k0_gluon = x*HQenergy;
    double kT_gluon = y*k0_gluon;
    double tauf = TauF(x, y, HQenergy);

    double result;
    if (k0_gluon < scale)
        return 1e-12;
    if (tauf < 1./scale)
        return 1e-12;

    result = 4./M_PI*Nc*AlphaS(kT_gluon, temp) * SplittingP(x) * qhat * pow(y,5) 
            * pow(sin((time)/2./tauf/HBARC),2)
            * pow(HQenergy*HQenergy/(y*y*HQenergy*HQenergy + HQmass_*HQmass_), 4)
            / (k0_gluon * k0_gluon * HBARC);
    if (result > maxdNg)
        maxdNg = result;
    //std::cout << x << " " << y << "  " << result <<"  " << maxdNgdxdydt << std::endl;
    return result;
}


double wrapper_dNgdxdydt(double*args, size_t dim, void *params)
{
    (void) (dim);
    struct gsl_NgIntegral_params *p = (struct gsl_NgIntegral_params*) params;
    double x = args[0];
    double y = args[1];

    double result = p->ptr->dNg_over_dxdydt(p->time, p->temp, p->HQenergy, x, y, p->maxdNgdxdydt);
    return result;
}

void RadiationHT::calculate(double time, double temp, double HQenergy, double& result, double &maxdNg)
{
    double xl[2] = {0,0};
    double xu[2] = {1,1};

    struct gsl_NgIntegral_params ps;
    ps.time = time;
    ps.temp = temp;
    ps.HQenergy = HQenergy;
    ps.maxdNgdxdydt = 0; // initialize the maxdNg before each integration
    ps.ptr = this;

    double res,err;
    const gsl_rng_type *T;
    gsl_rng* r;
    gsl_monte_function G;
    
    G.f = &wrapper_dNgdxdydt;
    G.dim = 2;
    G.params = &ps;

    size_t calls = 500000;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    {
        gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
        gsl_monte_vegas_integrate(&G, xl, xu, 2, 100000, r, s, &res, &err);
        int ctl = 0;

        do {
            gsl_monte_vegas_integrate(&G, xl, xu, 2, calls/5, r, s, &res, &err);
            ctl ++;
        } while (fabs(gsl_monte_vegas_chisq(s) -1.0) > 0.1 && ctl <=20);

        gsl_monte_vegas_free(s);
    }

    gsl_rng_free(r);

    result = res;
    maxdNg = ps.maxdNgdxdydt;
}


void RadiationHT::tabulate(size_t NEstart, size_t dNE)
{
    double time, temp, energy;
    double result, maxdNg;
    for (size_t i=0; i < Ntime_; ++i) {
      time = timeL_ + i*dtime_;
      for (size_t j=0; j<Ntemp_; ++j) {
        temp = tempL_ + j*dtemp_;
        for (size_t k=NEstart; k < (NEstart + dNE); ++ k) {
          energy = EL_ + k * dE_;
          calculate(time, temp, energy, result, maxdNg);
          RadTab[i][j][k] = result;
          Max_dNg[i][j][k] = maxdNg;
//   std::cout << time << " " << temp << " " << energy << " " << result << " " << maxdNg<< std::endl;
      }
    }
  } 
}


double RadiationHT::interpR(double time, double temp, double HQenergy)
{
    if (time > timeH_) {
        time = timeH_;
        std::cerr << "tabulate in time is not big enough! " << std::endl;
    }

    if (temp > tempH_) {
        temp = tempH_;
        std::cerr << "tabulate in temperature is not big enough!" << std::endl;
    }

    if (HQenergy > EH_) {
        HQenergy = EH_;
        std::cerr << "tabulate in energy is not big enough! " << std::endl;
    }

    double xTime, xTemp, xE, rTime, rTemp, rE;
    size_t iTime, iTemp, iE;

    xTime = (time - timeL_) / dtime_;  iTime = floor(xTime+1e-6); rTime = xTime - iTime;
    xTemp = (temp - tempL_) / dtemp_;  iTemp = floor(xTemp+1e-6); rTemp = xTemp - iTemp;
    xE = (HQenergy - EL_) / dE_;              iE = floor(xE+1e-6);       rE = xE - iE;
  
    return interpolate3d(&RadTab, iTime, rTime, iTemp, rTemp, iE, rE);
}




int main()
{
    RadiationHT rad = RadiationHT(4, "dNg_over_dxdydt.dat", true, true);

    double time=2.0;
    double temp=0.18;
    double HQenergy = 10;

    double result, maxdNg;
    rad.calculate(time, temp, HQenergy, result, maxdNg);
    std::cout << result << " " << maxdNg << std::endl;
    //rad.tabulate(0,1);
    for (int i=1; i<=200; i++)
      for (int j=1; j<=200; j++)
      {
        double x = i*(1./200);
        double y = j*(1./200);
        //double result = rad.dNg_over_dxdydt(time, temp, HQenergy, x, y);
        //std::cout << x << "  " << y << "  " <<rad.dNg_over_dxdydt(time, temp, HQenergy, x, y) << std::endl;
      }

}
