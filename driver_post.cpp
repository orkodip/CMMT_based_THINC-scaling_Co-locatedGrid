//POST PROCESSING OF THE RESULT
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#define TRUNC_L 1e-7	//volume fraction lower truncation criteria
#define TRUNC_U 1e-7	//volume fraction upper truncation criteria
#define EPS 1e-6
#define SMALL 1e-8
#define TOL 1e-6	//convergence criteria for BICGSTAB solver
using namespace std;
const int I=128,J=32;	//no of grids in i and j directions
const double WIDTH=0.4,HEIGHT=0.1;	//domain size
#include "func_head.h"
#include "mbase.cpp"
#include "BICGSTAB_SSOR.cpp"
#include "ivf_rec.cpp"
#include "clsvof.cpp"
#include "NS_g.cpp"
#include "post.cpp"
int main()
{
	int CNT=20;
	POST ms;
	ms.MBASE::ini(CNT,1.2,1000.0,1.77634e-4,0.011387,0.0,1.0e-4);	//count,rho_0,rho_1,mu_0,mu_1,sigma,dt
	ms.CLSVOF::ini(0.05715,0.05715,0.05715,0.05715);
	string fname;
	for(int COUNT=CNT;COUNT<=3200;COUNT+=20)
	{
		try	//solution initialization
		{
			fname="inter_"+to_string(COUNT);
			ms.read_bin(fname);
		}
		catch(int a)
		{
			cout<<"FATAL ERROR! INPUT FILE DOES NOT EXIST"<<endl;
			return 1;
		}
		ms.calc_xloc(COUNT);
		ms.calc_yloc(COUNT);
	}
	ms.write();
	return 0;
}
