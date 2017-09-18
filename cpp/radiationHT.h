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
    size_t Ntime_, Ntemp_, NE_;
    double timeL_, timeH_, tempL_, tempH_, EL_, EH_;
    double dtime_, dtemp_, dE_;
    
    boost::multi_array<double, 3> RadTab, Max_dNg;

    void save_to_file(std::string fileName, bool plain);
    void read_from_file(std::string fileName, bool plain);

    double AlphaS(double kT, double temp); // a fancy running coupling (depend on flavor)
    double SplittingP(double x);  // splitting function for quarks/gluons
    double TauF(double x, double y, double HQenergy); // formation time for emitted gluon

public:
    RadiationHT(int flavor, std::string fileName, bool plain, bool refresh);
    double dNg_over_dxdydt(double time, double temp, double HQnergy, double x, double y, double& maxdNg); // gluon emission distribution
    void calculate(double time, double temp, double HQenergy, double& result, double& maxdNg);//gluon emission prob, also the maximum integrand it encountered
    void tabulate(size_t NEstart, size_t dNE);
    double interpR(double time, double temp, double HQenergy); // interpolate emission prob
    double interpMax(double time, double temp, double HQenergy); // interpolate max_dNg
    double getNg(double time, double temp, double HQenergy, double qhat);
    double get_MaxdNg(double time, double temp, double HQenegry);
    double get_dNg(double time, double temp, double HQenergy, double x, double y);
    bool emitGluon(double time, double temp, double HQenergy, double qhat, double deltat);
    bool sampleGluon(double time, double temp, double HQenergy, double qhat, double deltat, std::vector<double> & gluon);
};

struct gsl_NgIntegral_params
{
    double time;
    double temp;
    double HQenergy;
    double maxdNgdxdydt;
    RadiationHT *ptr;
};


#endif
