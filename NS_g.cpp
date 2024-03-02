//NAVIER-STOKES SOLVER
//BOUNDARY CONDITIONS ARE AS FOLLOWS
//LEFT, RIGHT, AND BOTTOM BOUNDARY - NO SLIP AND NO PENETRATION
//TOP BOUNDARY - OUTFLOW (P=0)
//ARITHMETIC MEAN IS USED TO INTERPOLATE DENSITY
//HARMONIC MEAN IS USED TO INTERPOLATE VISCOSITY
class NS:public CLSVOF,public BICGSTAB
{
	protected:
	double **us,**vs;	//cell centered intermediate velocity field
	double **D_x,**D_y;	//cell centered diffusion term for x and y momentum equations

	void vel_bc();	//velocity boundary condition
	void P_bc();	//pressure boundary condition
	void updt_prev();	//store the solenoidal velocity field of nth time step
	void calc_curv();	//calculate curvature
	void diff();	//explicit viscous term
	void mom();	//solve the momentum equations
	void calc_fc();	//calculate face centered velocities
	void Press();	//solve the pressure Poisson equation
	void update();	//velocity and pressure update
	void continuity();	//calculate the continuity
	void max_CFL();	//calculate the max CFL number
	void max_GFN();	//determine the max grid Fourier number
	public:
			NS(); ~NS();
			void solve();
			void write_bin(int count);	//export intermediate file in binary format
			void read_bin(string fname);	//import intermediate data
};
NS::NS():BICGSTAB(I*J,1000)
{
	us=new double*[J+2];
	vs=new double*[J+2];
	D_x=new double*[J+1];
	D_y=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		us[j]=new double[I+2];
		vs[j]=new double[I+2];
		if(j<J+1)
		{
			D_x[j]=new double[I+1];
			D_y[j]=new double[I+1];
		}
	}
	cout<<"NS: MEMORY ALLOCATED"<<endl;
}
NS::~NS()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] us[j];
		delete[] vs[j];
		if(j<J+1)
		{
			delete[] D_x[j];
			delete[] D_y[j];
		}
	}
	delete[] us;
	delete[] vs;
	delete[] D_x;
	delete[] D_y;
	cout<<"NS: MEMORY RELEASED"<<endl;
}
void NS::vel_bc()
{
	for(int i=1;i<=I;i++)	//bottom and top boundary
	{
		v_NS[0][i]=0.0;	//no slip and no penetration
		u[0][i]=-u[1][i];
		v[0][i]=-v[1][i];

		u[J+1][i]=u[J][i];	//outflow
		v[J+1][i]=v[J][i];
	}
	for(int j=1;j<=J;j++)	//left and right boundary
	{
		u_EW[j][0]=0.0;	//no slip and no penetration
		u[j][0]=-u[j][1];
		v[j][0]=-v[j][1];

		u_EW[j][I]=0.0;	//no slip and no penetration
		u[j][I+1]=-u[j][I];
		v[j][I+1]=-v[j][I];
	}
	u[0][0]=-u[1][0]; u[0][I+1]=-u[1][I+1]; u[J+1][0]=-u[J+1][1]; u[J+1][I+1]=-u[J+1][I];	//corner points
	v[0][0]=-v[1][0]; v[0][I+1]=-v[1][I+1]; v[J+1][0]=-v[J+1][1]; v[J+1][I+1]=-v[J+1][I];
}
void NS::P_bc()
{
	for(int j=1;j<=J;j++)	//left and right boundary
	{
		P[j][0]=P[j][1];	//Neumann
		P[j][I+1]=P[j][I];	//Neumann
	}
	for(int i=1;i<=I;i++)	//bottom and top boundary
	{
		P[0][i]=P[1][i];	//Neumann
		P[J+1][i]=-P[J][i];	//Dirichlet
	}
}
void NS::updt_prev()
{
	for(int j=1;j<=J;j++)
		for(int i=0;i<=I;i++)
			u_EW_nm1[j][i]=u_EW[j][i];
	for(int j=0;j<=J;j++)
		for(int i=1;i<=I;i++)
			v_NS_nm1[j][i]=v_NS[j][i];
}
void NS::calc_curv()
{
	double phi_x,phi_y,phi_xx,phi_yy,phi_xy;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			phi_x=0.5*(Phi[j][i+1]-Phi[j][i-1])/dx;
			phi_y=0.5*(Phi[j+1][i]-Phi[j-1][i])/dy;
			phi_xx=(Phi[j][i+1]-2.0*Phi[j][i]+Phi[j][i-1])/pow(dx,2.0);
			phi_yy=(Phi[j+1][i]-2.0*Phi[j][i]+Phi[j-1][i])/pow(dy,2.0);
			phi_xy=0.25*(Phi[j+1][i+1]+Phi[j-1][i-1]-Phi[j+1][i-1]-Phi[j-1][i+1])/(dx*dy);
			if((pow((pow(phi_x,2.0)+pow(phi_y,2.0)),1.5))!=0)	//check dinominator
				KC[j][i]=-(pow(phi_y,2.0)*phi_xx-2.0*phi_x*phi_y*phi_xy+pow(phi_x,2.0)*phi_yy)/(pow((pow(phi_x,2.0)+pow(phi_y,2.0)),1.5));
			else KC[j][i]=0.0;
		}
	}
	for(int j=1;j<=J;j++)	//EW face traversal
	{
		for(int i=0;i<=I;i++)
		{
			if((i!=0)&&(i!=I))	//inner domain
				K_EW[j][i]=0.5*(KC[j][i]+KC[j][i+1]);
			else if(i==0)	//left boundary
				K_EW[j][0]=KC[j][1];
			else if(i==I)	//right boundary
				K_EW[j][I]=KC[j][I];
		}
	}
	for(int j=0;j<=J;j++)	//NS face traversal
	{
		for(int i=1;i<=I;i++)
		{
			if((j!=0)&&(j!=J))	//inner domain
				K_NS[j][i]=0.5*(KC[j][i]+KC[j+1][i]);
			else if(j==0)	//bottom boundary
				K_NS[0][i]=KC[1][i];
			else if(j==J)	//top boundary
				K_NS[J][i]=KC[J][i];
		}
	}
}
void NS::diff()
{
	double dx2i=pow(dx,-2),dy2i=pow(dy,-2),dxdyi=pow((dx*dy),-1);
	double mu_EW[2],mu_NS[2];	//viscosity at cell faces
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			mu_EW[0]=2.0*(mu[j][i]*mu[j][i-1])/(mu[j][i]+mu[j][i-1]); mu_EW[1]=2.0*(mu[j][i+1]*mu[j][i])/(mu[j][i+1]+mu[j][i]);
			mu_NS[0]=2.0*(mu[j][i]*mu[j-1][i])/(mu[j][i]+mu[j-1][i]); mu_NS[1]=2.0*(mu[j+1][i]*mu[j][i])/(mu[j+1][i]+mu[j][i]);
			D_x[j][i]=2.0*dx2i*(mu_EW[1]*(u[j][i+1]-u[j][i])-mu_EW[0]*(u[j][i]-u[j][i-1]))
					+dy2i*(mu_NS[1]*(u[j+1][i]-u[j][i])-mu_NS[0]*(u[j][i]-u[j-1][i]))
					+0.25*dxdyi*(mu_NS[1]*(v[j][i+1]+v[j+1][i+1]-v[j][i-1]-v[j+1][i-1])
									-mu_NS[0]*(v[j][i+1]+v[j-1][i+1]-v[j][i-1]-v[j-1][i-1]));
			D_y[j][i]=dx2i*(mu_EW[1]*(v[j][i+1]-v[j][i])-mu_EW[0]*(v[j][i]-v[j][i-1]))
					+2.0*dy2i*(mu_NS[1]*(v[j+1][i]-v[j][i])-mu_NS[0]*(v[j][i]-v[j-1][i]))
					+0.25*dxdyi*(mu_EW[1]*(u[j+1][i+1]+u[j+1][i]-u[j-1][i+1]-u[j-1][i])
									-mu_EW[0]*(u[j+1][i-1]+u[j+1][i]-u[j-1][i]-u[j-1][i-1]));
		}
	}
}
void NS::mom()
{
	diff();	//calculate diffusion term
	for(int j=1;j<=J;j++)	//explicit calculation
	{
		for(int i=1;i<=I;i++)
		{
			us[j][i]=A_x[j][i]+dt/rho_np1[j][i]*D_x[j][i];
			vs[j][i]=A_y[j][i]+dt/rho_np1[j][i]*D_y[j][i];
		}
	}
}
void NS::calc_fc()
{
	double rdx=dt/dx,rdy=dt/dy;
	double rho_fp;
	VECTOR n;	//interface normal vector
	double y_eff;	//effective height of the interface
	for(int j=1;j<=J;j++)	//calculation of face centered values excluding the boundaries
	{
		for(int i=1;i<=I;i++)
		{
			if(i<I)
			{
				rho_fp=2.0/(rho_np1[j][i]+rho_np1[j][i+1]);
				n.x=(Phi[j][i+1]-Phi[j][i])/dx;
				n.y=0.25*(Phi[j+1][i]+Phi[j+1][i+1]-Phi[j-1][i]-Phi[j-1][i+1])/dy;
				n=n.unit();
				y_eff=CY[j]-0.5*(Phi[j][i]+Phi[j][i+1])*n.y;
				u_EW[j][i]=0.5*(us[j][i]+us[j][i+1])+rdx*rho_fp*(sigma*K_EW[j][i]+(rho_0-rho_1)*Grav.y*y_eff)*(F[j][i+1]-F[j][i]);
			}
			if(j<J)
			{
				rho_fp=2.0/(rho_np1[j][i]+rho_np1[j+1][i]);
				n.x=0.25*(Phi[j][i+1]+Phi[j+1][i+1]-Phi[j][i-1]-Phi[j+1][i-1])/dx;
				n.y=(Phi[j+1][i]-Phi[j][i])/dy;
				n=n.unit();
				y_eff=Ym[j]-0.5*(Phi[j][i]+Phi[j+1][i])*n.y;
				v_NS[j][i]=0.5*(vs[j][i]+vs[j+1][i])+rdy*rho_fp*(sigma*K_NS[j][i]+(rho_0-rho_1)*Grav.y*y_eff)*(F[j+1][i]-F[j][i]);
			}
		}
	}
	for(int i=1;i<=I;i++)	//outflow top boundary
	{
		rho_fp=2.0/(rho_np1[J][i]+rho_np1[J+1][i]);
		n.x=0.25*(Phi[J][i+1]+Phi[J+1][i+1]-Phi[J][i-1]-Phi[J+1][i-1])/dx;
		n.y=(Phi[J+1][i]-Phi[J][i])/dy;
		n=n.unit();
		y_eff=Ym[J]-0.5*(Phi[J][i]+Phi[J+1][i])*n.y;
		v_NS[J][i]=vs[J][i]+rdy*rho_fp*(sigma*K_NS[J][i]+(rho_0-rho_1)*Grav.y*y_eff)*(F[J+1][i]-F[J][i]);
	}
}
void NS::Press()
{
	int cnt=0; R[cnt]=0;
	double bt=pow((dx/dy),2.0);
	double rdx=dx/dt,rdy=bt*dy/dt,rdx2=dx*dx/dt;
	double A_e,A_w,A_n,A_s;	//coefficients
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			A_w=2.0/(rho_np1[j][i]+rho_np1[j][i-1]);
			A_e=2.0/(rho_np1[j][i]+rho_np1[j][i+1]);
			A_s=bt*2.0/(rho_np1[j][i]+rho_np1[j-1][i]);
			A_n=bt*2.0/(rho_np1[j][i]+rho_np1[j+1][i]);
			X[(j-1)*I+i-1]=P[j][i];	//initial guess
			b[(j-1)*I+i-1]=-rdx2*((u_EW[j][i]-u_EW[j][i-1])/dx+(v_NS[j][i]-v_NS[j-1][i])/dy);
			if((i==1)&&(j==1))	//bottom left corner cell
			{
				C[cnt]=(j-1)*I+i-1;	//1,1 term
				A[cnt]=A_e+A_n;
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,1 term
				A[cnt]=-A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,2 term
				A[cnt]=-A_n;
				cnt++;
			}
			else if((i==1)&&(j>1)&&(j<J))	//left boundary
			{
				C[cnt]=(j-2)*I+i-1;	//1,j-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,j term
				A[cnt]=A_e+A_n+A_s;
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,j term
				A[cnt]=-A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//1,j+1 term
				A[cnt]=-A_n;
				cnt++;
			}
			else if((i==1)&&(j==J))	//top left corner cell
			{
				C[cnt]=(j-2)*I+i-1;	//1,J-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//1,J term
				A[cnt]=A_e+A_s+2.0*A_n;
				cnt++;
				C[cnt]=(j-1)*I+i;	//2,J term
				A[cnt]=-A_e;
				cnt++;
			}
			else if((j==J)&&(i>1)&&(i<I))	//top boundary
			{
				C[cnt]=(j-2)*I+i-1;	//i,J-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,J term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,J term
				A[cnt]=A_e+A_w+A_s+2.0*A_n;
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,J term
				A[cnt]=-A_e;
				cnt++;
			}
			else if((j==J)&&(i==I))	//top right corner cell
			{
				C[cnt]=(j-2)*I+i-1;	//I,J-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,J term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,J term
				A[cnt]=A_w+A_s+2.0*A_n;
				cnt++;
			}
			else if((i==I)&&(j>1)&&(j<J))	//right boundary
			{
				C[cnt]=(j-2)*I+i-1;	//I,j-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//I-1,j term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,j term
				A[cnt]=A_w+A_n+A_s;
				cnt++;
				C[cnt]=j*I+i-1;	//I,j+1 term
				A[cnt]=-A_n;
				cnt++;
			}
			else if((i==I)&&(j==1))	//bottom right corner cell
			{
				C[cnt]=(j-1)*I+i-2;	//I-1,1 term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//I,1 term
				A[cnt]=A_n+A_w;
				cnt++;
				C[cnt]=j*I+i-1;	//I,2 term
				A[cnt]=-A_n;
				cnt++;
			}
			else if((j==1)&&(i>1)&&(i<I))	//bottom boundary
			{
				C[cnt]=(j-1)*I+i-2;	//i-1,1 term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,1 term
				A[cnt]=A_e+A_w+A_n;
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,1 term
				A[cnt]=-A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,2 term
				A[cnt]=-A_n;
				cnt++;
			}
			else	//inside domain
			{
				C[cnt]=(j-2)*I+i-1;	//i,j-1 term
				A[cnt]=-A_s;
				cnt++;
				C[cnt]=(j-1)*I+i-2;	//i-1,j term
				A[cnt]=-A_w;
				cnt++;
				C[cnt]=(j-1)*I+i-1;	//i,j term
				A[cnt]=A_e+A_w+A_n+A_s;
				cnt++;
				C[cnt]=(j-1)*I+i;	//i+1,j term
				A[cnt]=-A_e;
				cnt++;
				C[cnt]=j*I+i-1;	//i,j+1 term
				A[cnt]=-A_n;
				cnt++;
			}
			R[(j-1)*I+i]=cnt;
		}
	}
	BICGSTAB::solve(C,R,A,X,b);	//SSOR preconditioned iterations
	for(int j=1;j<=J;j++)	//updation of the calculated P
		for(int i=1;i<=I;i++)
			P[j][i]=X[(j-1)*I+i-1];
}
void NS::update()
{
	double rdx=dt/dx,rdy=dt/dy;
	double rho_fp,rho_fm;
	VECTOR n;	//interface normal vector
	double y_eff_p,y_eff_m;	//effective height of the interface at the plus and minus face
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(i<I)
			{
				rho_fp=2.0/(rho_np1[j][i]+rho_np1[j][i+1]);
				u_EW[j][i]=u_EW[j][i]-rdx*(P[j][i+1]-P[j][i])*rho_fp;
			}
			rho_fp=2.0/(rho_np1[j][i]+rho_np1[j+1][i]);
			v_NS[j][i]=v_NS[j][i]-rdy*(P[j+1][i]-P[j][i])*rho_fp;
			rho_fp=2.0/(rho_np1[j][i]+rho_np1[j][i+1]);
			n.x=(Phi[j][i+1]-Phi[j][i])/dx;
			n.y=0.25*(Phi[j+1][i]+Phi[j+1][i+1]-Phi[j-1][i]-Phi[j-1][i+1])/dy;
			n=n.unit();
			y_eff_p=CY[j]-0.5*(Phi[j][i]+Phi[j][i+1])*n.y;
			rho_fm=2.0/(rho_np1[j][i]+rho_np1[j][i-1]);
			n.x=(Phi[j][i]-Phi[j][i-1])/dx;
			n.y=0.25*(Phi[j+1][i]+Phi[j+1][i-1]-Phi[j-1][i]-Phi[j-1][i-1])/dy;
			n=n.unit();
			y_eff_m=CY[j]-0.5*(Phi[j][i]+Phi[j][i-1])*n.y;
			u[j][i]=us[j][i]-0.5*rdx*((P[j][i+1]-P[j][i]-(sigma*K_EW[j][i]+(rho_0-rho_1)*Grav.y*y_eff_p)*(F[j][i+1]-F[j][i]))*rho_fp
							+(P[j][i]-P[j][i-1]-(sigma*K_EW[j][i-1]+(rho_0-rho_1)*Grav.y*y_eff_m)*(F[j][i]-F[j][i-1]))*rho_fm);
			rho_fp=2.0/(rho_np1[j][i]+rho_np1[j+1][i]);
			n.x=0.25*(Phi[j][i+1]+Phi[j+1][i+1]-Phi[j][i-1]-Phi[j+1][i-1])/dx;
			n.y=(Phi[j+1][i]-Phi[j][i])/dy;
			n=n.unit();
			y_eff_p=Ym[j]-0.5*(Phi[j][i]+Phi[j+1][i])*n.y;
			rho_fm=2.0/(rho_np1[j][i]+rho_np1[j-1][i]);
			n.x=0.25*(Phi[j][i+1]+Phi[j-1][i+1]-Phi[j][i-1]-Phi[j-1][i-1])/dx;
			n.y=(Phi[j][i]-Phi[j-1][i])/dy;
			n=n.unit();
			y_eff_m=Ym[j-1]-0.5*(Phi[j][i]+Phi[j-1][i])*n.y;
			v[j][i]=vs[j][i]-0.5*rdy*((P[j+1][i]-P[j][i]-(sigma*K_NS[j][i]+(rho_0-rho_1)*Grav.y*y_eff_p)*(F[j+1][i]-F[j][i]))*rho_fp
							+(P[j][i]-P[j-1][i]-(sigma*K_NS[j-1][i]+(rho_0-rho_1)*Grav.y*y_eff_m)*(F[j][i]-F[j-1][i]))*rho_fm);
		}
	}
}
void NS::continuity()
{
	double cont=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			cont+=(u_EW[j][i]-u_EW[j][i-1])/dx+(v_NS[j][i]-v_NS[j-1][i])/dy;	//inner domain
	cout<<"NS: count = "<<COUNT<<" cont = "<<abs(cont)<<endl;
}
void NS::solve()
{
	updt_prev();
	vel_bc();
	calc_curv();
	mom();
	calc_fc();
	Press(); P_bc();
	update();
	CLSVOF::prop_updt();
	COUNT++;
	//continuity(); mass_err(); max_CFL(); max_GFN();
	if((COUNT%10)==0) { continuity(); mass_err(); max_CFL(); }
}
void NS::max_CFL()
{
	double u_max=0.0,v_max=0.0;
	double cfl_x,cfl_y;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(abs(u[j][i])>u_max) u_max=abs(u[j][i]);
			if(abs(v[j][i])>v_max) v_max=abs(v[j][i]);
		}
	}
	cfl_x=u_max*dt/dx; cfl_y=v_max*dt/dy;
	if(cfl_x>cfl_y) cout<<"NS: CFL_cell = "<<cfl_x<<endl;
	else cout<<"NS: CFL_cell = "<<cfl_y<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			if(abs(u_EW[j][i])>u_max) u_max=abs(u_EW[j][i]);
			if(abs(v_NS[j][i])>v_max) v_max=abs(v_NS[j][i]);
		}
	}
	cfl_x=u_max*dt/dx; cfl_y=v_max*dt/dy;
	if(cfl_x>cfl_y) cout<<"NS: CFL_f = "<<cfl_x<<endl;
	else cout<<"NS: CFL_f = "<<cfl_y<<endl;
}
void NS::max_GFN()
{
	double max=0.0,gfn[3];
	double nu_EW,nu_NS;
	double comm=dt*(dx*dx+dy*dy)/(dx*dx*dy*dy);
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			gfn[0]=mu[j][i]/rho_np1[j][i]*comm;
			gfn[1]=0.5*(mu[j][i+1]/rho_np1[j][i+1]+mu[j][i]/rho_np1[j][i])*comm;
			gfn[2]=0.5*(mu[j+1][i]/rho_np1[j+1][i]+mu[j][i]/rho_np1[j][i])*comm;
			for(int k=0;k<3;k++) if(gfn[k]>max) max=gfn[k];
		}
	}
	cout<<"NS: MAX GFN = "<<max<<", comm = "<<comm<<endl;
}
void NS::write_bin(int count)
{
	string fname="inter_"+to_string(count);
	ofstream p_out(fname);
	MBASE::write_bin(p_out);
	p_out.close();
	cout<<"NS: INTERMEDIATE FILE OUTPUT SUCCESSFULL AT n = "<<count<<endl;
}
void NS::read_bin(string fname)
{
	ifstream p_in(fname);
	if(p_in.fail()) throw(0);	//file does not exist
	MBASE::read_bin(p_in);
	p_in.close();
	CLSVOF::den_ini();
	CLSVOF::prop_updt();
	cout<<"NS: SOLUTION INITIALIZED SUCCESSFULLY"<<endl;
}
