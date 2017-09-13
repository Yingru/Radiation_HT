#ifndef Rad_HT
#define Rad_HT

#include <iostream>
#include <iomanip>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <boost/multi_array.hpp>

#include "utility.h"

class RadiationHT
{
private:
    int HQflavor_;
    double HQmass_;
    int Ntime_, Ntemp_, NE_;
    double timeH_, tempL_, tempH_, EL_, EH_;
    double dtime_, dtemp_, dE_;
    
    boost::multi_array<double, 3> RadTab, Max_dNg;

    void save_to_file(std::string fileName, bool plain);
    void read_from_file(std::string fileName, bool plain);

    double AlphaS(double kT, double temp); // a fancy running coupling (depend on flavor)
    double SplittingP(double x);  // splitting function for quarks/gluons
    double TauF(double x, double y, double HQenergy); // formation time for emitted gluon

public:
    RadiationHT(int flavor, std::string fileName, bool plain, bool refresh);
    double dNg_over_dxdydt(double time, double temp, double HQnergy, double x, double y); // gluon emission distribution
    double Ngluon(double time, double temp, double HQenergy);//gluon emission prob
};

struct gsl_NgIntegral_params
{
    double time;
    double temp;
    double HQenergy;
    RadiationHT *ptr;
};


#endif
