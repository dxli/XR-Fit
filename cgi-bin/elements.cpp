// extract absorption coefficient for a given energy with given element
// using data files mass_abs_coeffs.dat, *.nff

#include "nph.h"
#include "elements.h"

using namespace std;

element::element (string ele0,double en0)
//get real part of an element (ele0) at a given photon energy en0
{
    if(check_ele(ele0)) return;
    ele=ele0;
    en=en0;
    lambda=E_to_l(en);
    string err1("Incorrect data file format! Exiting now!\n"),
    fn("mass_abs_coeffs.dat"),
    err2("Run with parameters: element_Symbol energy(KeV)\n");
    ifstream in0(fn.c_str());
    if (! in0.is_open())
    {
        cerr<<"Error in opening the source file: "<<fn<<endl;
        exit (1);
    }
    /*
    string fout0("tmpoutput.txt");
    ofstream out0(fout0.c_str());
    if (! out0.is_open())
      {
        cerr<< "Error in creating the target file: "<<fout0<<endl;
        exit (1);
      }
      */
    if( en<1. || en >100.) {
        cerr<<"Energy = "<<en<<" KeV? Maybe you specified in eV, while it should be in KeV\n";
    }
    //Find the element_symbol
    string linebuf;
    unsigned i=0;
    while(i<ele.size() && isalpha(ele.at(i))){
        if(i) ele.at(i)=tolower(ele.at(i)); else ele.at(i)=toupper(ele.at(i));
        if(i++>=1) break;
    }
    double el_z,el_a,rho0;
    string initial("begin");
    while(getline(in0,linebuf)){
        //Find element
        istringstream iss (linebuf);
        vector < string >vs;
        std::copy (istream_iterator < string >(iss),
                   istream_iterator < string >(), back_inserter (vs));
        if(vs.size()<3) continue;
        if(vs.at(0).compare(initial)) continue;
        if(vs.at(2).compare(ele)) continue;
        if(vs.size()<6) {
            cerr<<"Incorrect structure:\n";
            cout<<linebuf<<endl;
            return;
        }
        //cout<<"Element: "<<vs.at(1)<<"\tEnerge= "<<en<<" KeV";
        name=vs.at(1); //element name
        symbol=vs.at(2);
        if ( ! (stringstream(vs.at(3))>>el_z)) {
            cerr<<"Error reading Z\n";
            return;
        }
        if ( ! (stringstream(vs.at(4))>>el_a)) {
            cerr<<"Error reading A\n";
            return;
        }
        if ( ! (stringstream(vs.at(5))>>rho0)) {
            cerr<<"Error reading rho\n";
            return;
        }
        rho0 *= 1e-24;  // in g/angstrom^3 now
        break;
    }
    vector< vector<double> > abs0;
    vector<vector< vector<double> >::iterator> edges;
    i=0;
    while(getline(in0,linebuf)){
        if(!linebuf.size()) continue;
        string::size_type p0=linebuf.find("end",0);
        if(p0 != string::npos ) break;
        vector < double >vs;
        istringstream iss (linebuf);
        std::copy (istream_iterator < double >(iss),
                   istream_iterator < double >(), back_inserter (vs));
        if (vs.size()<2) continue;
        if(vs.at(1)>0.) vs.at(1) = log(vs.at(1));//we do interpolation of log
        if(abs0.size()) {
            if(abs0.at(abs0.size()-1).at(1) < vs.at(1)) {//edge
                //cout<<vs.at(0)<<'\t'<<vs.at(1)<<endl;
                if( vs.at(0) >= en) break;
                i=abs0.size();
            }
        }
        abs0.push_back(vs);
    }
    in0.close();
    if(i) {
        abs0.erase(abs0.begin(),abs0.begin()+i);
    }
    if(!abs0.size()){
        cerr<<"No absorption data read\n";
        return;
    }
    //cout<<"abs0.size()="<<abs0.size()<<endl;
    //for(i=0;i<abs0.size();i++) cout<<abs0.at(i).at(0)<<' '<<abs0.at(i).at(1)<<endl;
    //cout<<'\t'<<el_z<<' '<<el_a<<' '<<rho0<<" g/cm^3\n";
    mabs=exp(my_spline(abs0, en))*1e16; // in angstrom^2/g
    //cout<<"mabs1="<<abs1<<endl;
    f1=formFactor_f1();
    //beta=mu*lambda/(4.*M_PI);
    amass=el_a;
    density=rho0;
   // cout<<"name="<<name<<"\tdensity="<<density<<endl;
    rho_el=f1*N_Av*density*1e-24/amass;
    /*
    cout<<"x-ray energy= "<<en<<" KeV\tlambda= "<<1.2398424334e1/en<<" angstrom\n";
    cout<<"Mass absorption= "<<abs1<<" cm^2/g\t"<<"Linear absorption= "<<b<<" cm^-1\n";

    cout<<"Absorption Length= "<<1.e4/b<<" micro-m\tbeta= "<<b*1.2398424334e-7/en/(4.*M_PI)<<endl;
    cout<<ele<<" atomic form factor at q=0, Energy= "<<en<<"KeV: "<<formFactor_f1(ele,en)<<endl;
    */
}

double element::formFactor_f1()
// get the real part of form factor of the element for the given energy, using the *.nff files
{
    string s0=symbol;
    double en1=en*1000.;
    //cout<<en1<<endl;
    //std::transform(s0.begin(),s0.end(),s0.begin(),(int(*)(int)) tolower);
    s0.at(0) = tolower(s0.at(0));
    s0 += string(".nff");
    //cout<<s0<<endl;
    ifstream in0(s0.c_str());
    if(!in0.is_open()) {
        cerr<<"Can not open "<<s0<<endl;
        return(0.);
    }
    string linebuf;
    vector<vector<double> > nff;
    vector<vector<double> >::iterator pnff;
//    int il=0;
    double x0=0.,x1=0.,x2=0.;
    while(getline(in0,linebuf)){
        string::size_type loc=linebuf.find("E",0);
        if(loc != string::npos) {
	 if(loc<4) continue;
	}
        vector < double >vs;
        istringstream iss (linebuf);
        std::copy (istream_iterator < double >(iss),
                   istream_iterator < double >(), back_inserter (vs));
        if (vs.size()<3) continue;
        if (vs.at(1)< -10.) continue;
        x0=x1;x1=x2;x2=vs.at(1);
        if (nff.size()>=2) {//find the absorption edges
            pnff=nff.end();
            pnff--;
            if(  x1 <= x0-1e-3 && x1 <= x2- 1e-3) {
                if (vs.at(0) >= en1) break;
                nff.erase(nff.begin(),pnff);//new edge
            }
        }
        nff.push_back(vs);
    }
    in0.close();
   // cout<<"Size of nff "<<nff.size()<<endl;
   // for(int i=0;i<nff.size();i++)cout<<nff.at(i).at(0)<<' '<<nff.at(i).at(1)<<endl;
    return( my_spline(nff,en1));
}

compound::compound(string c0,double en0,double rho0)
//compound, get the real part of form factor and absorption, compound chemical formula in c0, energy in en0, bulk density in rho0 (g/cm^3)
{
rho_el=0.;
beta=0.;
    vector<string> el0;
    vector<int> n0;
    string c1=c0;
    //read formula
    unsigned j=0;
    while(j<c1.size()) { //parse the formula
        string ele;
	//cout<<c1.at(j)<<endl;
        if(! isalpha(c1.at(j))) {
            j++;
            continue;
        }
        ele.push_back(toupper(c1.at(j++)));
	if(j<c1.size()) {
	if( islower(c1.at(j))) ele.push_back(c1.at(j++));
	}
        if(check_ele(ele)) {
            continue; //ignore, if not a known element
        }
	el0.push_back(ele);
        string num;
        while(j<c1.size()) {
            if(!isdigit(c1.at(j))) break;
            num.push_back(c1.at(j++));
        }
        if(num.size()) n0.push_back(atoi(num.c_str()));else n0.push_back(1);
    }
    
    if(!n0.size()) {
    	mu=0.;
	beta=0.;
	rho_el=0.;
    return;
    }
    //found
    ostringstream f0;
    amass=0.;
    double mabs0=0.;
//    double mu0=0.;
    f1=0.;
    for(unsigned i=0;i<el0.size();i++){
        f0<<el0.at(i)<<n0.at(i);
        element ele1(el0.at(i),en0);
        els.push_back(ele1);
        amass +=ele1.amass*n0.at(i);
        f1 +=ele1.f1*n0.at(i);
    //cout<<ele1.f1<<endl;
        mabs0 += ele1.mabs*ele1.amass*n0.at(i);
    }
    formula=f0.str();
if(els.size()==1 && isalpha(c0.at(c0.size()-1))) {
// Element symbol, use density from data_base
rho=els.at(0).density;
} else rho=rho0*1e-24; //in angstrom^-3
    mabs=mabs0/amass;//in angstrom^2/g
    //cout<<mabs<<endl;
    mu = mabs*rho; // in angstrom^-1
    lambda=E_to_l(en0);
    beta = mu*lambda/(4.*M_PI);
    //cout<<f1<<endl;
    rho_el = f1 * N_Av * rho/amass;
}


double my_spline(vector<vector<double> > t0,double x)
// do the general interpolation
{
 vector<double> xi,yi;
 vector<vector<double> >::iterator pnff=t0.begin();
  while( pnff != t0.end()){
          xi.push_back(pnff->at(0));
          yi.push_back(pnff->at(1));
          pnff++;
  } 
  gsl_interp_accel *acc = gsl_interp_accel_alloc (); 
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, xi.size()); 
  gsl_spline_init(spline, &(xi[0]), &(yi[0]), xi.size());
          return(gsl_spline_eval(spline, x, acc));
}

double E_to_l(double en0)
//Energy to wavelength, KeV to angstrom
{
 double  e=1.602176e-19, h=6.626069e-34,  c=2.997925e8;
 return h*c/e*1.e7/en0;
}

int check_ele(string el0)
//if el0 is an element return 1
{
	string etable(" ac ag al ar as at au ba be bi b br ca cd ce cl c co cr cs cu dy er eu fe f fr ga gd ge he hf hg h ho i in ir k kr la li lu mg mn mo na nb nd ne ni n o os pa pb pd pm p po pr pt ra rb re rh rn ru sb sc se si sm s sn sr ta tb tc te th ti tl tm u v w xe yb y zn zr ");
	string el1=" "+el0+" ";
	el1.at(1)=tolower(el1.at(1));
	if( etable.find(el1.c_str(),0) == string::npos) return(1);
		else return(0);
}
