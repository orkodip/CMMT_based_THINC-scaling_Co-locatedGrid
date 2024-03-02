//INITIALIZATION OF THE VOLUME FRACTIONS AND LEVEL SETS FOR A RECTANGULAR COLUMN USING THINC METHOD
class INI	//initial volume fraction and level sets
{
	double *X,*Y,*CX,*CY,**Fa,**Phia;
	double dx,dy;	//grid size
	double beta,W,H,Xc,Yc;	//interface thickness parameter, height, width, and coordinates of the farthest corner point
	double om[3],eta[3];	//Gaussian weights and points
	void reinit(double **Phi,int **tag);	//determine the initial level set function from volume fractions
	public:
		INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double Beta,
						double w,double h,double XC,double YC);	//initialization of the input variables
		void LS(int **tag);	//determine the initial level sets
		void VF();	//determine the initial volume fractions
};
INI::INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double Beta,double w,double h,double XC,double YC)
{
	X=x; Y=y; CX=cx; CY=cy;
	Fa=fa; Phia=phia; beta=Beta;
	W=w; H=h; Xc=XC; Yc=YC;
	dx=X[2]-X[1]; dy=Y[2]-Y[1];
	om[0]=4.0/9.0; om[1]=om[2]=5.0/18.0;
	eta[0]=0.0; eta[1]=sqrt(3.0/5.0); eta[2]=-sqrt(3.0/5.0);
}
void INI::LS(int **tag)
{
	int i_cell,j_cell;	//index of cell containing the fartherst corner point
	for(int j=1;j<=J;j++)	//determine j index of the cell containing the farthest corner point
	{
		if((Y[j]>=Yc)&&(Y[j-1]<=Yc)) { j_cell=j; break; }
	}
	for(int i=1;i<=I;i++)	//determine i index of the cell containing the farthest corner point
	{
		if((X[i]>=Xc)&&(X[i-1]<=Xc)) { i_cell=i; break; }
	}
	for(int i=1;i<=i_cell;i++) { Phia[j_cell][i]=Yc-CY[j_cell]; tag[j_cell][i]=1; }
	for(int j=1;j<=j_cell;j++) { Phia[j][i_cell]=Xc-CX[i_cell]; tag[j][i_cell]=1; }
	reinit(Phia,tag);	//generate LS field of the entire domain
}
void INI::reinit(double **Phi,int **tag)
{
	double h=MAX2(dx,dy);
	double temp,a,b;
//--------------------INITIALIZATION SCHEME--------------------------------------------------
	for(int j=1;j<=J;j++)	//initialize level sets
	{
		for(int i=1;i<=I;i++)
		{
			if(tag[j][i]==0)
			{
				if((CX[i]<=Xc)&&(CY[j]<=Yc)) Phi[j][i]=100.0;	//inside fluid 1
				else Phi[j][i]=-100.0;	//outside fluid 1
			}
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phi[j][0]=Phi[j][1]; tag[j][0]=tag[j][1];
		Phi[j][I+1]=Phi[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phi[0][i]=Phi[1][i]; tag[0][i]=tag[1][i];
		Phi[J+1][i]=Phi[J][i]; tag[J+1][i]=tag[J][i];
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
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==I+1) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==J+1) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
void INI::VF()
{
	double sum,a_00,a_10,a_01,a_11,a_20,a_02,P;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			a_00=Phia[j][i];	//interface polynomial coefficients
			a_10=0.5*(Phia[j][i+1]-Phia[j][i-1])/dx;
			a_01=0.5*(Phia[j+1][i]-Phia[j-1][i])/dy;
			a_11=0.25*(Phia[j+1][i+1]-Phia[j+1][i-1]-Phia[j-1][i+1]+Phia[j-1][i-1])/(dx*dy);
			a_20=0.5*(Phia[j][i+1]-2.0*Phia[j][i]+Phia[j][i-1])/(pow(dx,2.0));
			a_02=0.5*(Phia[j+1][i]-2.0*Phia[j][i]+Phia[j-1][i])/(pow(dy,2.0));
			sum=0.0;
			for(int l2=0;l2<3;l2++)	//Gaussian quadrature to calculate actual volume fractions
			{
				for(int l1=0;l1<3;l1++)
				{
					P=(a_00+a_10*(0.5*dx*eta[l1])+a_01*(0.5*dy*eta[l2])+a_11*(0.5*dx*eta[l1])*(0.5*dy*eta[l2])
								+a_20*pow((0.5*dx*eta[l1]),2.0)+a_02*pow((0.5*dy*eta[l2]),2.0));	//calculation of interfacial polynomial
					sum+=om[l2]*om[l1]*tanh(beta*P);
				}
			}
			Fa[j][i]=0.5*(1.0+sum);
			if(abs(Fa[j][i])<=TRUNC_L) Fa[j][i]=0.0;	//truncation of insignificant values
			else if(abs(1.0-Fa[j][i])<=TRUNC_U) Fa[j][i]=1.0;
		}
	}
}
