//ABSTRACT BASE CLASS
class MBASE
{
	protected:
	double *Xm,*Ym;	//the mesh variables
	double *CX,*CY;	//the node variables
	double **u,**v,**P;	//cell centered X velocity, Y velocity, pressure
	double **u_EW,**v_NS;	//advection velocities at the cell faces
	double **u_EW_nm1,**v_NS_nm1;	//advection velocities at the cell faces at n-1 time level
	double **F,**Phi;	//volume fraction and level set function
	double **KC,**K_EW,**K_NS;	//interface curvature at cell centers and cell faces
	double **rho_np1,**rho_n,**mu;	//cell centered density and dynamic viscosity
	double **A_x,**A_y;	//cell centered advection term for X and Y momentum equations
	int *C,*R;	//CSR matrices
	double *A,*X,*b;

	double rho_0,rho_1;	//density for fluid 0 and fluid 1
	double mu_0,mu_1;	//viscosity for fluid 0 and fluid 1
	double sigma;	//surface tension coefficient
	VECTOR Grav;	//acceleration due to gravity
	double dx,dy,dt;	//grid spacing, time step
	int COUNT;
	void grid_gen();	//uniform grid generator
	public:
			MBASE(); ~MBASE();
			void ini(int count,double Rho_0,double Rho_1,double Mu_0,double Mu_1,double SIG,double Dt);	//initialize the problem
			void grid_write();	//export grid
			void write(int t);	//export velocity and pressure field
			void write_adv(int t);	//export advection field
			void write_den_vis(int t);	//export density and viscosity field
			void write_curv(int t);	//export curvature field
			void write_bin(ofstream &p_out);	//export intermediate data in binary format
			void read_bin(ifstream &p_in);	//import intermediate data
};
MBASE::MBASE()
{
	Xm=new double[I+1];
	Ym=new double[J+1];
	CX=new double[I+1];
	CY=new double[J+1];
	P=new double*[J+2];
	u_EW=new double*[J+2];
	u_EW_nm1=new double*[J+2];
	K_EW=new double*[J+2];
	u=new double*[J+2];
	v=new double*[J+2];
	v_NS=new double*[J+1];
	v_NS_nm1=new double*[J+1];
	K_NS=new double*[J+1];
	F=new double*[J+2];
	Phi=new double*[J+2];
	rho_np1=new double*[J+2];
	rho_n=new double*[J+2];
	mu=new double*[J+2];
	A_x=new double*[J+1];
	A_y=new double*[J+1];
	KC=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		P[j]=new double[I+2];
		u[j]=new double[I+2];
		v[j]=new double[I+2];
		u_EW[j]=new double[I+1];
		u_EW_nm1[j]=new double[I+1];
		K_EW[j]=new double[I+1];
		F[j]=new double[I+2];
		Phi[j]=new double[I+2];
		rho_np1[j]=new double[I+2];
		rho_n[j]=new double[I+2];
		mu[j]=new double[I+2];
		if(j<J+1)
		{
			v_NS[j]=new double[I+2];
			v_NS_nm1[j]=new double[I+2];
			K_NS[j]=new double[I+2];
			A_x[j]=new double[I+1];
			A_y[j]=new double[I+1];
			KC[j]=new double[I+1];
		}
	}

	//--------------------------------CSR matrices------------------------------------------

	int IMAX=I,JMAX=J;	//order of the pentadiagonal matrix
	int NNZ=5*IMAX*JMAX;	//no of non-zero elements in the sparse pentadiagonal periodic matrix
	C=new int[NNZ];
	R=new int[IMAX*JMAX+1];
	A=new double[NNZ];
	X=new double[IMAX*JMAX];
	b=new double[IMAX*JMAX];

	cout<<"MBASE: MEMORY ALLOCATED"<<endl;
}
MBASE::~MBASE()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] P[j];
		delete[] u[j];
		delete[] v[j];
		delete[] u_EW[j];
		delete[] u_EW_nm1[j];
		delete[] K_EW[j];
		delete[] F[j];
		delete[] Phi[j];
		delete[] rho_np1[j];
		delete[] rho_n[j];
		delete[] mu[j];
		if(j<J+1)
		{
			delete[] v_NS[j];
			delete[] K_NS[j];
			delete[] v_NS_nm1[j];
			delete[] A_x[j];
			delete[] A_y[j];
			delete[] KC[j];
		}
	}
	delete[] Xm;
	delete[] Ym;
	delete[] CX;
	delete[] CY;
	delete[] P;
	delete[] u;
	delete[] v;
	delete[] u_EW;
	delete[] u_EW_nm1;
	delete[] K_EW;
	delete[] v_NS;
	delete[] v_NS_nm1;
	delete[] K_NS;
	delete[] F;
	delete[] Phi;
	delete[] rho_np1;
	delete[] rho_n;
	delete[] mu;
	delete[] A_x;
	delete[] A_y;
	delete[] KC;

	//--------------------------------CSR matrices------------------------------------------

	delete[] C;
	delete[] R;
	delete[] A;
	delete[] X;
	delete[] b;

	cout<<"MBASE: MEMORY RELEASED"<<endl;
}
void MBASE::ini(int count,double Rho_0,double Rho_1,double Mu_0,double Mu_1,double SIG,double Dt)
{
	rho_0=Rho_0; rho_1=Rho_1;
	mu_0=Mu_0; mu_1=Mu_1;
	sigma=SIG;
	Grav.x=0.0; Grav.y=-9.81;
	dx=WIDTH/I; dy=HEIGHT/J;
	dt=Dt;
	COUNT=count;
	cout<<"MBASE: Rho_0 = "<<rho_0<<", Rho_1 = "<<rho_1<<endl;
	cout<<"MBASE: Mu_0 = "<<mu_0<<", Mu_1 = "<<mu_1<<endl;
	cout<<"MBASE: Sigma = "<<sigma<<", CURVATURE IS CALCULATED FROM LS FIELD"<<endl;
	cout<<"MBASE: dt = "<<dt<<endl;
	grid_gen();
}
void MBASE::grid_gen()	//uniform grid generation done here
{
	Xm[0]=Ym[0]=CX[0]=CY[0]=0.0;
	for(int i=1;i<=I;i++)	//generation of mesh and cell centers
	{
		Xm[i]=Xm[i-1]+dx;
		CX[i]=0.5*(Xm[i]+Xm[i-1]);
	}
	for(int j=1;j<=J;j++)
	{
		Ym[j]=Ym[j-1]+dy;
		CY[j]=0.5*(Ym[j]+Ym[j-1]);
	}
	cout<<"MBASE: GRID GENERATED"<<endl;
}
void MBASE::grid_write()
{
	ofstream p_out("mesh.dat");
	p_out<<"TITLE = \"MESH\""<<endl;
	p_out<<"FILETYPE = GRID"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK"<<endl;
	for(int j=0;j<=J;j++)	//print X co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=0;j<=J;j++)	//print Y co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Ym[j];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: GRID WRITE SUCCESSFULL"<<endl;
}
void MBASE::write(int t)
{
	string fname="uvp_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"FLOW AND PRESSURE FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"u\",\"v\",\"P\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<u[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<v[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<P[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_adv(int t)
{
	string fname="adv_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"ADVECTION FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"adv_x\",\"adv_y\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_x[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<A_y[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: ADVECTION FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_curv(int t)
{
	string fname="curv_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"CURVATURE FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"KC\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<KC[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: CURVATURE FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_den_vis(int t)
{
	string fname="den_vis_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"DENSITY AND VISCOSITY FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"rho_np1\",\"mu\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<rho_n[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<mu[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: DENSITY AND VISCOSITY FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_bin(ofstream &p_out)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_out.write((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			p_out.write((char *) &F[j][i],sizeof(F[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &Phi[j][i],sizeof(Phi[j][i]));
	cout<<"MBASE: INTERMEDIATE FILE OUTPUT SUCCESSFULL"<<endl;
}
void MBASE::read_bin(ifstream &p_in)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_in.read((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			p_in.read((char *) &F[j][i],sizeof(F[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &Phi[j][i],sizeof(Phi[j][i]));
	cout<<"MBASE: SOLUTION INITIALIZED SUCCESSFULLY"<<endl;
}
