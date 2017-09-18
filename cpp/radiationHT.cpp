#include <stdlib.h>
#include <time.h>
#include <random>
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
        fRad << "#   time     temp     HQenery     Ng/qhat     Max_dNgdxdydt/qhat" << std::endl;
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
        //std::cout << header << std::endl;
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
// arxiv: 1605.06447 eqn.14
// returns actually dNg_over_dxdydt/qhat, I figured this is the best way without all those hastles..
{
    double qhat = 4. * ALPHA * pow(temp, 3) / CF;
    double scale = M_PI*temp;
    double k0_gluon = x*HQenergy;
    double kT_gluon = y*k0_gluon;
    double tauf = TauF(x, y, HQenergy);

    double result;
    if (k0_gluon < scale)
        return 1e-12;
    if (tauf < 1./scale)
        return 1e-12;

    /*
    result = 4./M_PI*Nc*AlphaS(kT_gluon, temp) * SplittingP(x) * qhat * pow(y,5) 
            * pow(sin((time)/2./tauf/HBARC),2)
            * pow(HQenergy*HQenergy/(y*y*HQenergy*HQenergy + HQmass_*HQmass_), 4)
            / (k0_gluon * k0_gluon * HBARC);
    */

    result = 4./M_PI*Nc*AlphaS(kT_gluon, temp) * SplittingP(x) * pow(y,5)
            *pow(sin((time)/2./tauf/HBARC),2)
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

double RadiationHT::interpMax(double time, double temp, double HQenergy)
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
  
    return interpolate3d(&Max_dNg, iTime, rTime, iTemp, rTemp, iE, rE);
}


double RadiationHT::getNg(double time, double temp, double HQenergy, double qhat)
// get the average number of radiated gluons from a hard parton with energy E, temp T, require qhat input
{
    double delta_ng = interpR(time, temp, HQenergy);
    return delta_ng * qhat;
}



double RadiationHT::get_MaxdNg(double time, double temp, double HQenergy)
// get the largest integrand incountered during integration, usedful when use rejection sampling
{
    double maxNg = interpMax(time, temp, HQenergy);
    return maxNg;
}

double RadiationHT::get_dNg(double time, double temp, double HQenergy, double x, double y)
{
    double dum = 0;
    double result = dNg_over_dxdydt(time, temp, HQenergy, x, y, dum);
    return result;
}

bool RadiationHT::emitGluon(double time, double temp, double HQenergy, double qhat, double deltat)
// check whether a gluon can be emiited in this timestep
{ 
    double delta_ng = getNg(time, temp, HQenergy, qhat) * deltat;
    if (delta_ng > 1)
    {
        std::cerr << "Error! Gluon emission probability exceeds 1! " << std::endl;
    } 
    
    //else
   // {
   //     std::cout << "Gluon emission probability: " << delta_ng << std::endl;
   // }


    double xLow = M_PI * temp / HQenergy;
    double xHigh = 1.0;
    double xRange = xHigh - xLow;
    if ((xRange > eSMALL) && (randUniform() < delta_ng - eSMALL) && (2.*HQenergy*(HQenergy - M_PI*temp) > HQmass_*HQmass_))
        return true;
    else
        return false;
}

bool RadiationHT::sampleGluon(double time, double temp, double HQenergy, double qhat, double deltat, std::vector<double> & gluon)
// given time-interval, temp, HQenergy, qhat, deltat, sample emitted gluon momentum
// the gluon momentum should be retated into 
{
    gluon.resize(3);
    gluon[0] = 0.; gluon[1] = 0.; gluon[2] = 0.;


    double omega0 = M_PI * temp;
    double xLow = omega0 / HQenergy;
    double xHigh = 1.0;
    double xRange = xHigh - xLow;
    double maxf = get_MaxdNg(time, temp, HQenergy);

    // gluon energy Eg = x * HQenergy, gluon perp energy E_perp = y * x * HQenergy;
    double x, y, tauf, dNgdxdydt;
    int count = 0;

    do {
        x = xLow + randUniform() * xRange;
        y = randUniform();
        tauf = TauF(x, y, HQenergy);
        dNgdxdydt = get_dNg(time, temp, HQenergy, x, y);
        count ++;

        if (count > 1e5)
        {
            std::cout << "Give up loop at point 1 ..." << std::endl;
            return false;
        }
    } while((tauf < 1./omega0) || (maxf*randUniform() > dNgdxdydt));

//<<<<< successfully sampled Gluon <<<<<
    double gluon_kperp = x * y * HQenergy;
    double theta = 2. * M_PI * randUniform();
    gluon[0] = gluon_kperp * cos(theta);
    gluon[1] = gluon_kperp * sin(theta);
    gluon[2] = x * HQenergy * sqrt(1 - y*y);

// remove the momentum broadening from gluons
// in SCao's Langevin, using senario narrow 4, to narrow the emitted gluon
    int count2 = 0;
    double widthGluon = sqrt(Nc * qhat * tauf/2.);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0, widthGluon);
    std::vector<double> gluon_dk = {0.,0.,0.};
    double gluon_dkT;
    do {
        gluon_dk[0] = distribution(generator);
        gluon_dk[1] = distribution(generator);
        gluon_dkT = sqrt(gluon_dk[0] * gluon_dk[0] + gluon_dk[1] * gluon_dk[1]);
        count2 ++;

        if (count2 > 1e6)
        {
            std::cout << "Give up loop at point 2 ... " << gluon_dk[0] << " " << gluon_dk[1] << std::endl;
            gluon[0] = 0.; gluon[1] = 0.; gluon[2] = 0.;
            return false;
        }
    } while ((gluon_dkT > gluon_kperp) || (gluon_dkT*gluon_dkT - 2*gluon_dk[0]*gluon[0] - 2*gluon_dk[1]*gluon[1])>0);
            
    gluon[0] = gluon[0] - gluon_dk[0];
    gluon[1] = gluon[1] - gluon_dk[1];

    return true;
}

int main()
{
    RadiationHT rad = RadiationHT(4, "dNg_over_dxdydt.dat", true, false);

    double time=2.0;
    double temp=0.18;
    double HQenergy = 50;
    double deltat_lrf = 0.2;
    double qhat = M_PI * pow(temp, 3);
 
    double result, maxdNg;
    rad.calculate(time, temp, HQenergy, result, maxdNg);
    std::cout << "gluon emission probability: " << result*qhat*deltat_lrf << " " << maxdNg << std::endl;

    bool emit, success;
    int count1 = 0, count2=0, trial = 1e5;
    std::vector<double> gluon;

    for (int i=0; i < trial; i++)
    {
        emit = rad.emitGluon(time, temp, HQenergy, qhat, deltat_lrf);
        if (emit)
        {
            count1 ++;
            success = rad.sampleGluon(time, temp, HQenergy, qhat, deltat_lrf, gluon);
            if (success)   count2 ++;
        }
    }

    std::cout << "reality emitted gluon: " << 1.*count1/trial << " " <<  1.*count2/trial << std::endl;

}
