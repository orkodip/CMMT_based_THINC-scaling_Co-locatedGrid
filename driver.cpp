//SOLUTION OF NAVIER STOKES EQUATION IN A COLOCATED GRID USING THE BALANCED FORCE METHOD.
//BICGSTAB IS USED TO SOLVE THE LINEAR SYSTEMS.
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<csignal>
#define TRUNC_L 1e-7	//volume fraction lower truncation criteria
#define TRUNC_U 1e-7	//volume fraction upper truncation criteria
#define EPS 1e-6
#define SMALL 1e-8
#define TOL 1e-6	//convergence criteria for BICGSTAB solver
using namespace std;
int LC=1;	//loop condition
void control(int n)	//signal control function
{
	cout<<"PROGRAM INTERUPTED!"<<endl;
	cout<<"Enter 0 to stop: ";
	cin>>LC;
	if(LC==0) cout<<"STOP SIGNAL RECIEVED. PROGRAM WILL STOP AFTER SIMULATING CURRENT TIME-STEP."<<endl;
	else {LC=1; cout<<"SIMULATION CONTINUES."<<endl;}
}
const int I=128,J=32;	//no of grids in i and j directions
const double WIDTH=0.4,HEIGHT=0.1;	//domain size
#include "func_head.h"
#include "mbase.cpp"
#include "BICGSTAB_SSOR.cpp"
#include "ivf_rec.cpp"
#include "clsvof.cpp"
#include "NS_g.cpp"
int main()
{
	signal(SIGINT,control);	//define signal and its function (FOR MANUALLY CONTROLLING THE PROGRAM DURING RUNTIME)
	int CNT=0;
	NS ms;
	ms.MBASE::ini(CNT,1.2,1000.0,1.77634e-4,0.011387,0.0,1.0e-4);	//count,rho_0,rho_1,mu_0,mu_1,sigma,dt
	ms.CLSVOF::ini(0.05715,0.05715,0.05715,0.05715);
	ms.grid_write();
	ms.lsvf_write(0);
	ms.write_den_vis(0);
	for(int COUNT=CNT;(LC && (COUNT<3500));COUNT++)	//manually controlled loop
	{
		ms.CLSVOF::solve();
		ms.NS::solve();
		if((((COUNT+1)%20)==0)||(COUNT==0))
		{
			ms.write(COUNT+1);
			ms.lsvf_write(COUNT+1);
			ms.write_bin(COUNT+1);
		}
	}
	return 0;
}
