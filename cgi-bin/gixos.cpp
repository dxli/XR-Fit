#include "gixosnph.h"

double multilayer::cT(double t)
//Transimission factoi T(sin\theta)
{
//    double cosk0=t*t;
 //   cosk0 =  cosk0*(1./2.-cosk0*(1./24.-1./720.*cosk0)); // 1 - cos(t)
    //double cosk0=1-a;
    //c=sqrt(a+a);
    //double sink0=sqrt(cosk0*(2.0-cosk0));
    complex<double> sink0 =n_bulk0+t*t;
    complex<double> sink1 =sink0 - n_bulk1;
    sink0=sqrt(sink0);
    sink1=sqrt(sink1);
    //cosk0=1-cosk0;
    //sink1=(sink1 - n_bulk1)*(1.+n_bulk1); // (1 -a0) /(1-nk)=(1-a0)(1+nk+nk^2)=1 -(a0 +a0nk-nk-nk^2)
    //complex<double> cosk1=(1. - sink1)/cosk0;
    //a0 = sqrt(b0+b0);
    //sink1= sqrt(sink1*(2.0-sink1))/sink0;
    sink0=2.*sink0/(sink0+(complex<double>(1.,0.)-n_bulk1)*sink1);
    return(sink0.real()*sink0.real()+sink0.imag()*sink0.imag());
}

double multilayer::yet_qz(double t)
// yet=k_B*T/(2 Pi \gamma)*qz^2
{
double qz=k0*(sin(alpha0)+t);
//cout<<"k0="<<k0<<" "<<qz<<" yeta="<<(Boltzmann_kB*temperature_T*1.e20/(2*M_PI*gamma)*qz*qz)<<endl;
return(Boltzmann_kB*temperature_T/(2*M_PI*gamma)*qz*qz);
}


double multilayer::sip_sink(double sink)
// structure independent parameter, all
{
double yeta0=yet_qz(sink);
double beta0=asin(sink);
double qr=k0*sqrt(cos(beta0)*cos(beta0)+cos(alpha0)*cos(alpha0)-2*cos(alpha0)*cos(beta0)*cos(dth));
beta0=n_bulk1.real()*M_PI/(lambda*lambda)*gsl_sf_gamma(1-0.5*yeta0)/qr; // constants
double ans=pow(qr*sigmad,yeta0);
yeta0 *= M_PI;
ans*=2.*sin(yeta0)/yeta0;
//cout<<"sink="<<sink<<" qr="<<qr<<" "<<cT0*ans/(gamma)<<endl;
ans *= cT(sink)*Boltzmann_kB*temperature_T/sin(alpha0)*beta0*beta0;
return(ans/(gamma*fresnel(sink)));
//return(cT0*ans/(gamma));
}

double multilayer::bulk(double sink)
// bulk scattering of water
{
double a0=sin(alpha0);
complex<double> qt=sqrt(a0*a0-n_bulk1)*k0;
double beta0=n_bulk1.real()*M_PI/(lambda*lambda); // constants
double ans= cT(sink)*Boltzmann_kB*temperature_T*kappa;
return(ans * beta0*beta0*0.5/qt.imag());
}

