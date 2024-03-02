//CONSISTENT TRANSPORT CLSVOF ALGORITHM USING THINC-LS SCHEME FOR VOLUME FRACTIONS AND SL SCHEME FOR LEVEL SETS.
//BOUNDARY CONDITIONS ARE AS FOLLOWS.
//LEFT, RIGHT, AND BOTTOM BOUNDARY - NO SLIP AND NO PENETRATION
//TOP BOUNDARY - OUTFLOW
class CLSVOF:public MBASE
{
	int **tag;	//interface cell tag
	int **LS_tag;	//cells having fixed LS
	VECTOR **XT;	//departure points
	double **F1,**F2;	//intermediate Volume Fractions
	double **Phi1,**Phi2;	//intermediate level set functions
	double **A_x1,**A_x2,**A_y1,**A_y2;	//intermediate advection terms for X and Y momentum equations
	class poly2	//stores the 2nd order interface polynomial
	{
		public:
				int tag;
				double a_00,a_10,a_01,a_11,a_20,a_02,LSc;	//coefficients of interface polynomial
	};
	poly2 **psi;	//polynomial surface of the interface
	double om[3],eta[3];	//Gaussian weights and points
	double beta;	//parameter to control slope and thickness of interface jump
	double mass_act;	//actual mass

	void updt_ghost(double **Phi);	//update the ghost cells of cell centered field
	void cell_tag();	//interface cell tagging
	void calc_dp();	//departure point calculation
	void poly_2nd(double **Fa,double **Phia);	//2nd order psi calculation
	double LS_corr(double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02);	//enforce mass conservation to the level set function
	double vol_frac_flux(int flag_xy,int flag_f,int i,int j);	//volume fraction flux
	double U_face(int flag_xy,int flag_f,int i,int j,double **Phia,double V);	//interpolated velocity at cell face
	//flag_xy->(x,y)=(0,1); flag_f->(W/S,E/N)=(0,1)
	void F_adv(int flag,double **Fa,double **A_xa,double **A_ya);	//volume fraction advection equation
	void LS_adv(int flag);	//LS advection equation
	void reinit(double **Phia,double **Fa);	//LS reinitialization algorithm
	public:
			CLSVOF(); ~CLSVOF();
			void ini(double w,double h,double xc,double yc);
			void den_ini();	//initialize density field at n+1 time step
			void solve();	//CLSVOF advection algorithm
			void prop_updt();	//update density and advection field for next time step
			void mass_err();	//calculate mass error
			void lsvf_write(int t);	//tecplot file output
};
CLSVOF::CLSVOF()
{
	tag=new int*[J+2];
	LS_tag=new int*[J+2];
	Phi1=new double*[J+2];
	Phi2=new double*[J+2];
	XT=new VECTOR*[J+1];
	psi=new poly2*[J+1];
	F1=new double*[J+1];
	F2=new double*[J+1];
	A_x1=new double*[J+1];
	A_x2=new double*[J+1];
	A_y1=new double*[J+1];
	A_y2=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		tag[j]=new int[I+2];
		LS_tag[j]=new int[I+2];
		Phi1[j]=new double[I+2];
		Phi2[j]=new double[I+2];
		if(j<J+1)
		{
			XT[j]=new VECTOR[I+1];
			psi[j]=new poly2[I+1];
			F1[j]=new double[I+1];
			F2[j]=new double[I+1];
			A_x1[j]=new double[I+1];
			A_x2[j]=new double[I+1];
			A_y1[j]=new double[I+1];
			A_y2[j]=new double[I+1];
		}
	}
	cout<<"CLSVOF: MEMORY ALLOCATED"<<endl;
}
CLSVOF::~CLSVOF()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] tag[j];
		delete[] LS_tag[j];
		delete[] Phi1[j];
		delete[] Phi2[j];
		if(j<J+1)
		{
			delete[] XT[j];
			delete[] psi[j];
			delete[] F1[j];
			delete[] F2[j];
			delete[] A_x1[j];
			delete[] A_x2[j];
			delete[] A_y1[j];
			delete[] A_y2[j];
		}
	}
	delete[] tag;
	delete[] LS_tag;
	delete[] Phi1;
	delete[] Phi2;
	delete[] XT;
	delete[] psi;
	delete[] F1;
	delete[] F2;
	delete[] A_x1;
	delete[] A_x2;
	delete[] A_y1;
	delete[] A_y2;
	cout<<"CLSVOF: MEMORY RELEASED"<<endl;
}
void CLSVOF::ini(double w,double h,double xc,double yc)
{
	beta=3.0/dx;
	cout<<"CLSVOF: beta*dx = "<<beta*dx<<endl;
	om[0]=4.0/9.0; om[1]=om[2]=5.0/18.0;
	eta[0]=0.0; eta[1]=sqrt(3.0/5.0); eta[2]=-sqrt(3.0/5.0);
	INI ms(Xm,Ym,CX,CY,F,Phi,beta,w,h,xc,yc);
	ms.LS(tag);
	ms.VF();	//initial exact volume fraction is calculated here
	updt_ghost(F);
	mass_act=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass_act+=F[j][i];
	for(int j=1;j<=J;j++)	//initialize density, viscosity, and advection field (including ghost nodes)
	{
		for(int i=1;i<=I;i++)
		{
			rho_n[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			A_x[j][i]=u[j][i];
			A_y[j][i]=v[j][i];
		}
	}
	for(int j=1;j<=J;j++)	//left and right ghost cell values
	{
		rho_n[j][0]=rho_n[j][1];
		mu[j][0]=mu[j][1];
		rho_n[j][I+1]=rho_n[j][I];
		mu[j][I+1]=mu[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost cell values
	{
		rho_n[0][i]=rho_n[1][i];
		mu[0][i]=mu[1][i];
		rho_n[J+1][i]=rho_n[J][i];
		mu[J+1][i]=mu[J][i];
	}
}
void CLSVOF::prop_updt()
{
	for(int j=0;j<=J+1;j++)	//including ghost nodes
	{
		for(int i=0;i<=I+1;i++)
		{
			rho_n[j][i]=rho_np1[j][i];
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			if((i>=1)&&(i<=I)&&(j>=1)&&(j<=J))	//only inner domain
			{
				A_x[j][i]=u[j][i];
				A_y[j][i]=v[j][i];
			}
		}
	}
}
void CLSVOF::updt_ghost(double **Phia)
{
	for(int j=1;j<=J;j++)	//left and right ghost nodes (Neumann bc)
	{
		Phia[j][0]=Phia[j][1];
		Phia[j][I+1]=Phia[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost nodes (Neumann bc)
	{
		Phia[0][i]=Phia[1][i];
		Phia[J+1][i]=Phia[J][i];
	}
}
void CLSVOF::den_ini()
{
	for(int j=0;j<=J+1;j++)	//calculate density field at n+1 step (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
}
void CLSVOF::cell_tag()
{
	for(int j=1;j<=J;j++)	//reinitialization of the cell tags
		for(int i=1;i<=I;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//tagging algorithm
	{
		for(int i=1;i<=I;i++)
		{
			if((F[j][i]>TRUNC_L)&&(F[j][i]<(1.0-TRUNC_U)))	//interfacial cell
			{
				tag[j][i]=1;	//tag interfacial cell
				tag[j][i-1]=1; tag[j][i+1]=1;	//tag neighbouring cells
				tag[j-1][i]=1; tag[j+1][i]=1;
				tag[j-1][i-1]=1; tag[j-1][i+1]=1;
				tag[j+1][i-1]=1; tag[j+1][i+1]=1;
			}
		}
	}
}
void CLSVOF::calc_dp()
{
	double K1x,K1y,K2x,K2y;
	double xs,ys;	//intermediate point
	double alp_x,alp_y;	//lever rule ratio
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			XT[j][i]=VECTOR();	//reinitialization
			if(tag[j][i]==1)	//calculation of departure point using Euler-Cauchy RK2 method
			{
				K1x=-0.5*dt*(u_EW[j][i]+u_EW[j][i-1]);
				K1y=-0.5*dt*(v_NS[j][i]+v_NS[j-1][i]);
				xs=CX[i]+0.5*K1x; ys=CY[j]+0.5*K1y;
				alp_x=(Xm[i]-xs)/dx; alp_y=(Ym[j]-ys)/dy;
				K2x=-0.5*dt*((1.0-alp_x)*(3.0*u_EW[j][i]-u_EW_nm1[j][i])+alp_x*(3.0*u_EW[j][i-1]-u_EW_nm1[j][i-1]));	//Adams-Bashforth scheme + Lever rule
				K2y=-0.5*dt*((1.0-alp_y)*(3.0*v_NS[j][i]-v_NS_nm1[j][i])+alp_y*(3.0*v_NS[j-1][i]-v_NS_nm1[j-1][i]));
				XT[j][i].x=CX[i]+K1x+K2x;	//departure point
				XT[j][i].y=CY[j]+K1y+K2y;
			}
		}
	}
}
void CLSVOF::poly_2nd(double **Fa,double **Phia)
{
	for(int j=1;j<=J;j++)	//reinitialization
	{
		for(int i=1;i<=I;i++)
		{
			psi[j][i].tag=0;
			psi[j][i].a_00=0.0;
			psi[j][i].a_10=0.0;
			psi[j][i].a_01=0.0;
			psi[j][i].a_11=0.0;
			psi[j][i].a_20=0.0;
			psi[j][i].a_02=0.0;
			psi[j][i].LSc=0.0;
		}
	}
	for(int j=1;j<=J;j++)	//psi calculation for interface cells
	{
		for(int i=1;i<=I;i++)
		{
			if(tag[j][i]==1)
			{
				if((Fa[j][i]>TRUNC_L)&&(Fa[j][i]<(1.0-TRUNC_U)))	//interfacial cell
				{
					psi[j][i].tag=1;
					psi[j][i].a_00=Phia[j][i];	//interface polynomial coefficients
					psi[j][i].a_10=0.5*(Phia[j][i+1]-Phia[j][i-1])/dx;
					psi[j][i].a_01=0.5*(Phia[j+1][i]-Phia[j-1][i])/dy;
					psi[j][i].a_11=0.25*(Phia[j+1][i+1]-Phia[j+1][i-1]-Phia[j-1][i+1]+Phia[j-1][i-1])/(dx*dy);
					psi[j][i].a_20=0.5*(Phia[j][i+1]-2.0*Phia[j][i]+Phia[j][i-1])/(pow(dx,2.0));
					psi[j][i].a_02=0.5*(Phia[j+1][i]-2.0*Phia[j][i]+Phia[j-1][i])/(pow(dy,2.0));
					psi[j][i].LSc=LS_corr(Fa[j][i],psi[j][i].a_00,psi[j][i].a_10,psi[j][i].a_01,psi[j][i].a_11,psi[j][i].a_20,psi[j][i].a_02);	//LS correction
				}
			}
		}
	}
}
double CLSVOF::LS_corr(double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02)
{
	double P[9],A[9];
	double gamma,temp,D=-1.0,C=2.0*(Fa-0.5);	//initial guess = -1.0
	double x1,y1;
	double func,func1;	//function and its derivative
	int cnt=0;	//no of iterations of NR method
	for(int l2=0;l2<3;l2++)	//calculation of gamma
	{
		for(int l1=0;l1<3;l1++)
		{
			x1=(0.5*dx)*eta[l1]; y1=(0.5*dy)*eta[l2];
			P[l2*3+l1]=(a_00+a_10*x1+a_01*y1+a_11*x1*y1+a_20*pow(x1,2.0)+a_02*pow(y1,2.0));	//calculation of interfacial polynomial
			temp=beta*P[l2*3+l1];
			if((l1==0)&&(l2==0)) gamma=temp;	//initialize gamma in the 1st iteration
			else if(gamma>temp) gamma=temp;
		}
	}
	gamma=-gamma+EPS;	//final value of gamma
	for(int l2=0;l2<3;l2++)	//calculation of A_g
		for(int l1=0;l1<3;l1++)
			A[l2*3+l1]=tanh(beta*P[l2*3+l1]+gamma);
	do	//Newton-Raphson method
	{
		func=func1=0.0;	//reinitialization
		temp=D;	//store value of previous iteration (variable is reused)
		for(int l2=0;l2<3;l2++)	//calculation of function and its derivatives
		{
			for(int l1=0;l1<3;l1++)
			{
				func+=om[l1]*om[l2]*(A[l2*3+l1]+D)/(1.0+A[l2*3+l1]*D);
				func1+=om[l1]*om[l2]*(1.0-pow(A[l2*3+l1],2.0))/(pow((1.0+A[l2*3+l1]*D),2.0));
			}
		}
		func-=C;
		D-=func/func1;
		cnt++;
		if(func>0.0) cout<<"CLSVOF: ERROR IN FUNC"<<endl;
		if(func1<0.0) cout<<"CLSVOF: ERROR IN FUNC1"<<endl;
	}
	while(abs(D-temp)>=EPS);
	return ((atanh(D)+gamma)/beta);
}
double CLSVOF::vol_frac_flux(int flag_xy,int flag_f,int i,int j)
{
	if(psi[j][i].tag==0) cout<<"NOT AN INTERFACIAL CELL!"<<endl;
	double P,x1,y1,sum=0.0;
	for(int l1=0;l1<3;l1++)	//Gauss quadrature
	{
		if(flag_xy==0)	//flux in x direction
		{
			if(flag_f==1) x1=0.5*dx;	//east face
			else if(flag_f==0) x1=-0.5*dx;	//west face
			y1=0.5*dy*eta[l1];
		}
		else if(flag_xy==1)	//flux in y direction
		{
			if(flag_f==1) y1=0.5*dy;	//north face
			else if(flag_f==0) y1=-0.5*dy;	//south face
			x1=0.5*dx*eta[l1];
		}
		P=psi[j][i].a_00+psi[j][i].a_10*x1+psi[j][i].a_01*y1+psi[j][i].a_11*x1*y1+psi[j][i].a_20*pow(x1,2.0)+psi[j][i].a_02*pow(y1,2.0);
		sum+=om[l1]*0.5*(1.0+tanh(beta*(P+psi[j][i].LSc)));
	}
	return sum;
}
double CLSVOF::U_face(int flag_xy,int flag_f,int i,int j,double **Phia,double V)	//QUICK + FOU schemes
{
	if(abs(V)<=EPS) return 0.0;
	int sum_tag;	//sum of tags
	if(flag_xy==0)	//interpolation in X direction
	{
		if(flag_f==0)	//west face
		{
			if((i>2)&&(i<I))	//may use QUICK scheme
			{
				sum_tag=tag[j][i]+tag[j][i+1]+tag[j][i-1]+tag[j][i-2];
				if(sum_tag==0) return (0.0625*(3.0*Phia[j][i-1]*(3.0+SGN(V))+3.0*Phia[j][i]*(3.0-SGN(V))
								-Phia[j][i-2]*(1.0+SGN(V))-Phia[j][i+1]*(1.0-SGN(V))));	//QUICK scheme
			}
			return (0.5*(Phia[j][i-1]*(1+SGN(V))+Phia[j][i]*(1-SGN(V))));	//FOU scheme
		}
		else if(flag_f==1)	//east face
		{
			if((j>1)&&(j<J-1))	//may use QUICK scheme
			{
				sum_tag=tag[j][i]+tag[j][i+1]+tag[j][i+2]+tag[j][i-1];
				if(sum_tag==0) return (0.0625*(3.0*Phia[j][i]*(3.0+SGN(V))+3.0*Phia[j][i+1]*(3.0-SGN(V))
								-Phia[j][i-1]*(1.0+SGN(V))-Phia[j][i+2]*(1.0-SGN(V))));	//QUICK scheme
			}
			return (0.5*(Phia[j][i]*(1+SGN(V))+Phia[j][i+1]*(1-SGN(V))));	//FOU scheme
		}
		else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
	}
	else if(flag_xy==1)	//interpolation in Y direction
	{
		if(flag_f==0)	//south face
		{
			if((j>2)&&(j<J))	//may use QUICK scheme
			{
				sum_tag=tag[j][i]+tag[j+1][i]+tag[j-1][i]+tag[j-2][i];
				if(sum_tag==0) return (0.0625*(3.0*Phia[j-1][i]*(3.0+SGN(V))+3.0*Phia[j][i]*(3.0-SGN(V))
								-Phia[j-2][i]*(1.0+SGN(V))-Phia[j+1][i]*(1.0-SGN(V))));	//QUICK scheme
			}
			return (0.5*(Phia[j-1][i]*(1+SGN(V))+Phia[j][i]*(1-SGN(V))));	//FOU scheme
		}
		else if(flag_f==1)	//north face
		{
			if((j>1)&&(j<J-1))	//may use QUICK scheme
			{
				sum_tag=tag[j][i]+tag[j+1][i]+tag[j+2][i]+tag[j-1][i];
				if(sum_tag==0) return (0.0625*(3.0*Phia[j][i]*(3.0+SGN(V))+3.0*Phia[j+1][i]*(3.0-SGN(V))
								-Phia[j-1][i]*(1.0+SGN(V))-Phia[j+2][i]*(1.0-SGN(V))));	//QUICK scheme
			}
			return (0.5*(Phia[j][i]*(1+SGN(V))+Phia[j+1][i]*(1-SGN(V))));	//FOU scheme
		}
		else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
	}
	else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
}
void CLSVOF::F_adv(int flag,double **Fa,double **A_xa,double **A_ya)
{
	double F_advx[2],F_advy[2];	//volume fraction advection fluxes
	double rho_advx[2],rho_advy[2];	//density advection fluxes
	double u_fx[2],u_fy[2];	//advected velocity for X momentum equation
	double v_fx[2],v_fy[2];	//advected velocity for Y momentum equation
	double rho_cell;	//cell density
	int iup,jup,flag_f;	//upwind cell index and face flag
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			//----------advected velocity for X and Y momentum equations----------------------------
			if(i==1) u_fx[0]=0.0;	//left boundary (no penetration)
			else u_fx[0]=U_face(0,0,i,j,A_xa,u_EW[j][i-1]);
			if(i==I) u_fx[1]=0.0;	//right boundary (no penetration)
			else u_fx[1]=U_face(0,1,i,j,A_xa,u_EW[j][i]);
			if(j==1) u_fy[0]=0.0;	//bottom boundary (no slip)
			else u_fy[0]=U_face(1,0,i,j,A_xa,v_NS[j-1][i]);
			if(j==J) u_fy[1]=A_xa[J][i];	//top boundary (outflow)
			else u_fy[1]=U_face(1,1,i,j,A_xa,v_NS[j][i]);

			if(i==1) v_fx[0]=0.0;	//left boundary (no slip)
			else v_fx[0]=U_face(0,0,i,j,A_ya,u_EW[j][i-1]);
			if(i==I) v_fx[1]=0.0;	//right boundary (no slip)
			else v_fx[1]=U_face(0,1,i,j,A_ya,u_EW[j][i]);
			if(j==1) v_fy[0]=0.0;	//bottom boundary (no penetration)
			else v_fy[0]=U_face(1,0,i,j,A_ya,v_NS[j-1][i]);
			if(j==J) v_fy[1]=A_ya[J][i];	//top boundary (outflow)
			else v_fy[1]=U_face(1,1,i,j,A_ya,v_NS[j][i]);
			//----------solution for single phase cells--------------------------------------------
			if(tag[j][i]==0)
			{
				if(F[j][i]<TRUNC_L)	//density flux for empty cells
				{
					rho_advx[0]=rho_0*u_EW[j][i-1];
					rho_advx[1]=rho_0*u_EW[j][i];
					rho_advy[0]=rho_0*v_NS[j-1][i];
					rho_advy[1]=rho_0*v_NS[j][i];
				}
				else if(F[j][i]>(1.0-TRUNC_U))	//density flux for completely filled cells
				{
					rho_advx[0]=rho_1*u_EW[j][i-1];
					rho_advx[1]=rho_1*u_EW[j][i];
					rho_advy[0]=rho_1*v_NS[j-1][i];
					rho_advy[1]=rho_1*v_NS[j][i];
				}
				else cout<<"CLSVOF: ERROR IN TAGGING ALGORITHM"<<endl;
				if(flag==1)	//step 1 of TVDRK3
				{
					F1[j][i]=F[j][i];
					A_x1[j][i]=rho_n[j][i]*A_x[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
									+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
									-rho_n[j][i]*A_x[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy));
					A_y1[j][i]=rho_n[j][i]*A_y[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
									+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
									-rho_n[j][i]*A_y[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy));
					A_x1[j][i]/=rho_n[j][i];	//extract intermediate velocity field
					A_y1[j][i]/=rho_n[j][i];
				}
				else if(flag==2)	//step 2 of TVDRK3
				{
					F2[j][i]=F[j][i];
					A_x2[j][i]=0.75*rho_n[j][i]*A_x[j][i]+0.25*(rho_n[j][i]*A_x1[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
									+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
									-rho_n[j][i]*A_x1[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
					A_y2[j][i]=0.75*rho_n[j][i]*A_y[j][i]+0.25*(rho_n[j][i]*A_y1[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
									+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
									-rho_n[j][i]*A_y1[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
					A_x2[j][i]/=rho_n[j][i];	//extract intermediate velocity field
					A_y2[j][i]/=rho_n[j][i];
				}
				else if(flag==3)	//step 3 of TVDRK3
				{
					F[j][i]=F[j][i];
					rho_np1[j][i]=rho_n[j][i];	//update density field
					A_x[j][i]=0.33334*rho_n[j][i]*A_x[j][i]
							+0.66666*(rho_n[j][i]*A_x2[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
							+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
							-rho_n[j][i]*A_x2[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
					A_y[j][i]=0.33334*rho_n[j][i]*A_y[j][i]
							+0.66666*(rho_n[j][i]*A_y2[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
							+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
							-rho_n[j][i]*A_y2[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
					A_x[j][i]/=rho_np1[j][i];	//extract intermediate velocity field
					A_y[j][i]/=rho_np1[j][i];
				}
				continue;
			}
			//----------volume fraction flux in X direction--------------------------------------------
			if(i==1)	//left boundary
			{
				F_advx[0]=0.0;	//wall boundary
				if(u_EW[j][i]>0.0) {iup=i; flag_f=1;}
				else {iup=i+1; flag_f=0;}
				if((abs(Fa[j][iup])<TRUNC_L)||(abs(u_EW[j][i])<=EPS)) F_advx[1]=0.0;
				else if(abs(1.0-Fa[j][iup])<=TRUNC_U) F_advx[1]=1.0;	//for completely filled cell
				else F_advx[1]=vol_frac_flux(0,flag_f,iup,j);
			}
			else if(i==I)	//right boundary
			{
				if(u_EW[j][i-1]>0.0) {iup=i-1; flag_f=1;}
				else {iup=i; flag_f=0;}
				if((abs(Fa[j][iup])<TRUNC_L)||(abs(u_EW[j][i-1])<=EPS)) F_advx[0]=0.0;
				else if(abs(1.0-Fa[j][iup])<=TRUNC_U) F_advx[0]=1.0;	//for completely filled cell
				else F_advx[0]=vol_frac_flux(0,flag_f,iup,j);
				F_advx[1]=0.0;	//wall boundary
			}
			else	//inner domain
			{
				if(u_EW[j][i-1]>0.0) {iup=i-1; flag_f=1;}
				else {iup=i; flag_f=0;}
				if((abs(Fa[j][iup])<TRUNC_L)||(abs(u_EW[j][i-1])<=EPS)) F_advx[0]=0.0;
				else if(abs(1.0-Fa[j][iup])<=TRUNC_U) F_advx[0]=1.0;	//for completely filled cell
				else F_advx[0]=vol_frac_flux(0,flag_f,iup,j);
				if(u_EW[j][i]>0.0) {iup=i; flag_f=1;}
				else {iup=i+1; flag_f=0;}
				if((abs(Fa[j][iup])<TRUNC_L)||(abs(u_EW[j][i])<=EPS)) F_advx[1]=0.0;
				else if(abs(1.0-Fa[j][iup])<=TRUNC_U) F_advx[1]=1.0;	//for completely filled cell
				else F_advx[1]=vol_frac_flux(0,flag_f,iup,j);
			}
		//--------------volume fraction flux in Y direction------------------------------------------
			if(j==1)	//bottom boundary
			{
				F_advy[0]=0.0;	//wall boundary
				if(v_NS[j][i]>0.0) {jup=j; flag_f=1;}
				else {jup=j+1; flag_f=0;}
				if((abs(Fa[jup][i])<TRUNC_L)||(abs(v_NS[j][i])<=EPS)) F_advy[1]=0.0;
				else if(abs(1.0-Fa[jup][i])<=TRUNC_U) F_advy[1]=1.0;	//for completely filled cell
				else F_advy[1]=vol_frac_flux(1,flag_f,i,jup);
			}
			else if(j==J)	//top boundary
			{
				if(v_NS[j-1][i]>0.0) {jup=j-1; flag_f=1;}
				else {jup=j; flag_f=0;}
				if((abs(Fa[jup][i])<TRUNC_L)||(abs(v_NS[j-1][i])<=EPS)) F_advy[0]=0.0;
				else if(abs(1.0-Fa[jup][i])<=TRUNC_U) F_advy[0]=1.0;	//for completely filled cell
				else F_advy[0]=vol_frac_flux(1,flag_f,i,jup);
				F_advy[1]=0.0;	//wall boundary
			}
			else	//inner domain
			{
				if(v_NS[j-1][i]>0.0) {jup=j-1; flag_f=1;}
				else {jup=j; flag_f=0;}
				if((abs(Fa[jup][i])<TRUNC_L)||(abs(v_NS[j-1][i])<=EPS)) F_advy[0]=0.0;
				else if(abs(1.0-Fa[jup][i])<=TRUNC_U) F_advy[0]=1.0;	//for completely filled cell
				else F_advy[0]=vol_frac_flux(1,flag_f,i,jup);
				if(v_NS[j][i]>0.0) {jup=j; flag_f=1;}
				else {jup=j+1; flag_f=0;}
				if((abs(Fa[jup][i])<TRUNC_L)||(abs(v_NS[j][i])<=EPS)) F_advy[1]=0.0;
				else if(abs(1.0-Fa[jup][i])<=TRUNC_U) F_advy[1]=1.0;	//for completely filled cell
				else F_advy[1]=vol_frac_flux(1,flag_f,i,jup);
			}
		//---------------density fluxes---------------------------------------------------------------
			rho_advx[0]=((rho_1-rho_0)*F_advx[0]+rho_0)*u_EW[j][i-1];
			rho_advx[1]=((rho_1-rho_0)*F_advx[1]+rho_0)*u_EW[j][i];
			rho_advy[0]=((rho_1-rho_0)*F_advy[0]+rho_0)*v_NS[j-1][i];
			rho_advy[1]=((rho_1-rho_0)*F_advy[1]+rho_0)*v_NS[j][i];
		//-------------------------------------------------------------------------------------------
			if(flag==1)	//step 1 of TVDRK3
			{
				F1[j][i]=F[j][i]+dt*((u_EW[j][i-1]*F_advx[0]-u_EW[j][i]*F_advx[1])/dx+(v_NS[j-1][i]*F_advy[0]-v_NS[j][i]*F_advy[1])/dy
								-F[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy));
				A_x1[j][i]=rho_n[j][i]*A_x[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
								+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
								-rho_n[j][i]*A_x[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy));
				A_y1[j][i]=rho_n[j][i]*A_y[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
								+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
								-rho_n[j][i]*A_y[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy));
				rho_cell=rho_1*F1[j][i]+rho_0*(1.0-F1[j][i]);
				A_x1[j][i]/=rho_cell;	//extract intermediate velocity field
				A_y1[j][i]/=rho_cell;
				if(F1[j][i]<=TRUNC_L)	//truncation of unrealistic values
				{
					if(F1[j][i]<-0.001) cout<<"i = "<<i<<", j = "<<j<<"F1 = "<<F1[j][i]<<endl;
					F1[j][i]=0.0;
				}
				else if(F1[j][i]>=(1.0-TRUNC_U))
				{
					if(abs(F1[j][i])>1.001) cout<<"i = "<<i<<", j = "<<j<<"F1 = "<<F1[j][i]<<endl;
					F1[j][i]=1.0;
				}
			}
			else if(flag==2)	//step 2 of TVDRK3
			{
				F2[j][i]=0.75*F[j][i]+0.25*(F1[j][i]+dt*((u_EW[j][i-1]*F_advx[0]-u_EW[j][i]*F_advx[1])/dx+(v_NS[j-1][i]*F_advy[0]
								-v_NS[j][i]*F_advy[1])/dy-F1[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				rho_cell=rho_1*F1[j][i]+rho_0*(1.0-F1[j][i]);
				A_x2[j][i]=0.75*rho_n[j][i]*A_x[j][i]+0.25*(rho_cell*A_x1[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
								+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
								-rho_cell*A_x1[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				A_y2[j][i]=0.75*rho_n[j][i]*A_y[j][i]+0.25*(rho_cell*A_y1[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
								+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
								-rho_cell*A_y1[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				rho_cell=rho_1*F2[j][i]+rho_0*(1.0-F2[j][i]);
				A_x2[j][i]/=rho_cell;	//extract intermediate velocity field
				A_y2[j][i]/=rho_cell;
				if(F2[j][i]<=TRUNC_L)	//truncation of unrealistic values
				{
					if(F2[j][i]<-0.001) cout<<"i = "<<i<<", j = "<<j<<"F2 = "<<F2[j][i]<<endl;
					F2[j][i]=0.0;
				}
				else if(F2[j][i]>=(1.0-TRUNC_U))
				{
					if(abs(F2[j][i])>1.001) cout<<"i = "<<i<<", j = "<<j<<"F2 = "<<F2[j][i]<<endl;
					F2[j][i]=1.0;
				}
			}
			else if(flag==3)	//step 3 of TVDRK3
			{
				F[j][i]=0.33334*F[j][i]+0.66666*(F2[j][i]+dt*((u_EW[j][i-1]*F_advx[0]-u_EW[j][i]*F_advx[1])/dx+(v_NS[j-1][i]*F_advy[0]
								-v_NS[j][i]*F_advy[1])/dy-F2[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				rho_cell=rho_1*F2[j][i]+rho_0*(1.0-F2[j][i]);
				A_x[j][i]=0.33334*rho_n[j][i]*A_x[j][i]+0.66666*(rho_cell*A_x2[j][i]+dt*((u_fx[0]*rho_advx[0]-u_fx[1]*rho_advx[1])/dx
								+(u_fy[0]*rho_advy[0]-u_fy[1]*rho_advy[1])/dy
								-rho_cell*A_x2[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				A_y[j][i]=0.33334*rho_n[j][i]*A_y[j][i]+0.66666*(rho_cell*A_y2[j][i]+dt*((v_fx[0]*rho_advx[0]-v_fx[1]*rho_advx[1])/dx
								+(v_fy[0]*rho_advy[0]-v_fy[1]*rho_advy[1])/dy
								-rho_cell*A_y2[j][i]*((u_EW[j][i-1]-u_EW[j][i])/dx+(v_NS[j-1][i]-v_NS[j][i])/dy)));
				rho_cell=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
				A_x[j][i]/=rho_cell;	//extract intermediate velocity field
				A_y[j][i]/=rho_cell;
				if(F[j][i]<=TRUNC_L)	//truncation of unrealistic values
				{
					if(F[j][i]<-0.001) cout<<"i = "<<i<<", j = "<<j<<"F = "<<F[j][i]<<endl;
					F[j][i]=0.0;
				}
				else if(F[j][i]>=(1.0-TRUNC_U))
				{
					if(abs(F[j][i])>1.001) cout<<"i = "<<i<<", j = "<<j<<"F = "<<F[j][i]<<endl;
					F[j][i]=1.0;
				}
				rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);	//update density field
			}
		}
	}
	if(flag==3) { updt_ghost(F); updt_ghost(rho_np1); }
}
void CLSVOF::LS_adv(int flag)
{
	double xt,yt;	//departure point in local coordinates
	double temp;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			LS_tag[j][i]=0;	//reinitialization
			if(tag[j][i]==1)
			{
				//--------------------update LS values iff departure cell contains the interface---------------------------
				if(psi[j][i].tag==0) continue;	//departure cell in not an interfacial cell
				xt=XT[j][i].x-CX[i]; yt=XT[j][i].y-CY[j];	//departure point in local coordinate system
				if(flag==1)	//step 1 of TVDRK3
				{
					LS_tag[j][i]=1;
					Phi1[j][i]=psi[j][i].LSc+psi[j][i].a_00+psi[j][i].a_10*xt+psi[j][i].a_01*yt
							+psi[j][i].a_11*xt*yt+psi[j][i].a_20*pow(xt,2.0)+psi[j][i].a_02*pow(yt,2.0);
				}
				else if(flag==2)	//step 2 of TVDRK3
				{
					LS_tag[j][i]=1;
					temp=psi[j][i].LSc+psi[j][i].a_00+psi[j][i].a_10*xt+psi[j][i].a_01*yt+psi[j][i].a_11*xt*yt
							+psi[j][i].a_20*pow(xt,2.0)+psi[j][i].a_02*pow(yt,2.0);	//variable is reused here
					Phi2[j][i]=0.75*Phi[j][i]+0.25*temp;
				}
				else if(flag==3)	//step 3 of TVDRK3
				{
					LS_tag[j][i]=1;
					temp=psi[j][i].LSc+psi[j][i].a_00+psi[j][i].a_10*xt+psi[j][i].a_01*yt+psi[j][i].a_11*xt*yt
							+psi[j][i].a_20*pow(xt,2.0)+psi[j][i].a_02*pow(yt,2.0);	//variable is reused here
					Phi[j][i]=0.33334*Phi[j][i]+0.66666*temp;
				}
			}
		}
	}
	//-----------------REINITIALIZATION OF THE LS FIELD---------------------------
	if(flag==1)	//step 1 of TVDRK3
		reinit(Phi1,F1);
	else if(flag==2)	//step 2 of TVDRK3
		reinit(Phi2,F2);
	else if(flag==3)	//step 3 of TVDRK3
		reinit(Phi,F);

}
void CLSVOF::reinit(double **Phia,double **Fa)
{
	double h=MAX2(dx,dy);
	double temp,a,b;
	//--------------------INITIALIZATION SCHEME (BASED ON ADVECTED LS FIELD)--------------------------------------------------
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if((LS_tag[j][i]==1)&&((Fa[j][i]<0.01)||(Fa[j][i]>0.94))) LS_tag[j][i]=0;	//LS re-tagging
			if(LS_tag[j][i]==0)	Phia[j][i]=(Fa[j][i]-0.5)*100.0;	//rest LS for untagged cells
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phia[j][0]=Phia[j][1]; LS_tag[j][0]=LS_tag[j][1];
		Phia[j][I+1]=Phia[j][I]; LS_tag[j][I+1]=LS_tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phia[0][i]=Phia[1][i]; LS_tag[0][i]=LS_tag[1][i];
		Phia[J+1][i]=Phia[J][i]; LS_tag[J+1][i]=LS_tag[J][i];
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(LS_tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phia[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phia[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phia[j][i])==1.0) a=MIN2(Phia[j][i+1],Phia[j][i-1]);
					else a=MAX2(Phia[j][i+1],Phia[j][i-1]);
				}
				if(j==0) b=Phia[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phia[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phia[j][i])==1.0) b=MIN2(Phia[j+1][i],Phia[j-1][i]);
					else b=MAX2(Phia[j+1][i],Phia[j-1][i]);
				}
				if(SGN(Phia[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phia[j][i]=MIN2(Phia[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phia[j][i]=MAX2(Phia[j][i],temp);
				}
			}
		}
	}
}
void CLSVOF::solve()
{
	cell_tag(); calc_dp();
	poly_2nd(F,Phi);
	F_adv(1,F,A_x,A_y); LS_adv(1);
	poly_2nd(F1,Phi1);
	F_adv(2,F1,A_x1,A_y1); LS_adv(2);
	poly_2nd(F2,Phi2);
	F_adv(3,F2,A_x2,A_y2); LS_adv(3);
}
void CLSVOF::lsvf_write(int t)
{
	string fname="ls_vol_frac_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::mass_err()
{
	double mass=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass+=F[j][i];
	cout<<"CLSVOF: MASS ERROR = "<<(mass-mass_act)<<endl;
}
