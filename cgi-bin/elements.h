#include<cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
int check_ele(string);
double E_to_l(double);
double my_spline(vector<vector<double> >, double);
class element {
	public:
	double f1,mabs,rho_el,en,amass,density,lambda,ratio;
	string name,symbol;
	element(){};
	element(string ele0,double en0);
	double mu(){
		return( mabs*density);
	};
	double beta(){
		return(mu()*lambda/(4.*M_PI));
	};
	/*
	double mabs(){ return(mabs); };
	double f1(){return(f1);};
	double amass(){return(amass);};
	*/
	int error(){return(! (symbol.size()));};
	private:
		string ele;
	double formFactor_f1();
	
};


class compound {
public:
double rho,f1,beta,mu,mabs,amass,rho_el,lambda;
string formula;
vector<element> els;
 compound(string c0,double en0, double rho0); // (formula, x-ray energy, density (g/cm^3))
 /*
 double f1(){return f1;};
 double beta(){return beta;};
 */
 };
