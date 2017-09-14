#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <iostream>

// >>>>>>>>>>>> constants >>>>>>>>>>>>>>>>>>>
const double cMass = 1.27;
const double bMass = 4.19;

const int Nc=3;
const double CF = 4.0/3;
const double ALPHA = 1.047;
const double HBARC = 0.1973;
const double fixAlphas = 0.3;



const double TINIT = 0.0;
const int HQener_gn = 200;
const int t_gn = 10;
const int temp_gn = 120;

const double HQener_max = 100.0;
const double t_min = 12.0;
const double t_max = 14.0;
const double temp_max = 0.75;
const double temp_min = 0.15;


// >>>>>>>>>>>> interpolate function >>>>>>>>>>>>>>
double interpolate3d(boost::multi_array<double, 3> *A, const int ni, const double ri,
                                                      const int nj, const double rj,
                                                      const int nk, const double rk)
{
    double wi[2] = {1.-ri, ri}, wj[2] = {1.-rj, rj}, wk[2] = {1.-rk, rk};
    double result = 0;
    for (int i=0; i<2; ++i) 
      for (int j=0; j<2; ++j)
        for (int k=0; k<2; ++k)
        { 
          result += (*A)[ni+i][nj+j][nk+k] * wi[i] * wj[j] * wk[k];
        }
    return result;
}
#endif
