#pragma once
#include <string>
using namespace std;

class Task
{
    public:
        Task();
        ~Task();
        
        int id;
        int progress;
        int maxisteps;
	int start;
        unsigned int nslab;
        string date;
        string email;
        string data;

	string formula0;
	string formula1;
	string formula12;
        int started;
        int running;
        int error;
        bool done;
	double rho0,rho1,rho12;//density of bulk0, bulk1, and maximum density allowed
	double rhoin0,rhoin1,rhoin12;//density from input
	double beta0,beta1,beta12;//beta of bulk0, bulk1, and maximum density allowed
	double slab,energy;
	double qmin,qmax; //qz range
	double sigma0;
	double dth,alpha0,gamma,sigmad;//gixos
};

