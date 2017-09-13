#include <iostream>
#include <iomanip> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include "utility.h"

using namespace std;

// a function of the strong coupling constant, 
// temp is the medium temperature, alphas = alphas(pi * temp) if kT < pi*temp due to the screening effect
double AlphaS(double kT, double temp)
{
    if (kT < M_PI*temp) 
        kT = M_PI * temp;
   
    int nflavor;
    double lambdas;
    if (kT < cMass) {
        nflavor = 3;
        lambdas = 0.2;
    } else if (kT < bMass) {
        nflavor = 4;
        lambdas = 0.172508;
    } else {
        nflavor = 5;
        lambdas = 0.130719;
    }

    double result = 4.*M_PI/((11.0 - 2.0*nflavor/3.0)*(2.0*log(kT/lambdas)));
    //std::cout << kT << " " << lambdas << " " << nflavor <<" " << log(kT/lambdas) << " " << result << std::endl;
    return result;
}



// the splitting function
double SplittingP(double x)
{
    double result = (2. - 2.*x + x*x) * CF / x;
    return result;
}


// gluon formation time
double TauF(double x, double y, double HQenergy, double HQmass)
{
    double k0_gluon = x * HQenergy;
    double kT_gluon = y * k0_gluon;
    return 2.0 * k0_gluon * (1.-x) / (pow(kT_gluon, 2) + pow(x*HQmass, 2));
}


struct Ng_params
{
    double HQenergy;
    double temp;
    double HQmass;
    double time;
};


// single gluon radiation formular from Wang's paper
//double dNg_over_dxdydt(double* args, double temp, double HQenergy, double HQmass)
// *args=[x,y], *params=Ng_params{HQenergy, temp, HQmass, time}

double dNg_over_dxdydt(double* args, size_t dim, void *params)
{
    (void) (dim);

    struct Ng_params *p = (struct Ng_params*) params;
    double HQenergy = p->HQenergy;
    double temp = p->temp;
    double HQmass = p->HQmass;
    double time = p->time;

    double x = args[0]; // x is the fraction of energy carried by emitted gluon
    double y = args[1]; // y is the ratio between kT / k0 => kT = x*y*HQenergy

    double qhat = 4.0*alpha*pow(temp,3)/CF;

    double scale = M_PI * temp;
    double k0_gluon = x * HQenergy;
    double kT_gluon = y * k0_gluon;
    double tauf = TauF(x, y, HQenergy, HQmass);

    double result;
    if (k0_gluon < scale)
        return 1e-12;
    
    if (tauf < 1./scale)
        return 1e-12;

    result = 4./M_PI*Nc*AlphaS(kT_gluon, temp) * SplittingP(x) * qhat * pow(y,5)
            * pow(sin((time-TINIT)/2./tauf/HBARC),2) 
            * pow(HQenergy*HQenergy/(y*y*HQenergy*HQenergy + HQmass*HQmass), 4) 
            / (k0_gluon*k0_gluon * HBARC);
   
    //cout << time << " " << temp << " " << HQenergy << " " << x << " " << y << " " << result << endl;
    return result;
}

double Ng_gluon(double HQenergy, double HQmass, double temp, double time)
{
    double res, err;
    double xl[2] = {0,0};
    double xu[2] = {1,1};

    const gsl_rng_type *T;
    gsl_rng * r;
    gsl_monte_function G;

    struct Ng_params ps;
    ps.HQenergy = HQenergy;
    ps.HQmass = HQmass;
    ps.temp = temp;
    ps.time = time;

    G.f = &dNg_over_dxdydt;
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
        } while (fabs(gsl_monte_vegas_chisq(s) -1.0) > 0.2 && ctl <=20);

        gsl_monte_vegas_free(s);
    }

    gsl_rng_free(r);
    return res;

}


int main()
{
    double time = 2.0;
    double temp0 = 0.15;
    double Energy0 = 10;
    double HQmass = 1.27;
    double dE = 10;
    double result;
    double temp, Energy;

    for (int i=0; i<10; i++)
    {
        temp = temp0 + i * 0.03;
        Energy = Energy0;
        for (int j=0; j<10; j++)
        {
            result = Ng_gluon(Energy, HQmass, temp, time);
            cout << setw(10) << time << setw(10) << temp << setw(10) << Energy << setw(10) << result << endl;
            Energy += dE;
        }
    }

}
