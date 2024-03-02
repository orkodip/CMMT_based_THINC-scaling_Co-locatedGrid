//BiCG-stab WITH SSOR(m) PRECONDITIONING
class BICGSTAB
{
	double *r0,*r,*p,*s,*m,*n;	//vectors required
	double *y,*t;	//solution of preconditioner
	double alpha,beta,gama,omega,sum,norm;	//variables required
	int N,MAX;	//number of elements in the vectors, and max no of iterations

	void SSOR(int *C,int *R,double *A,double *xp,double *bp,int M=1);	//SSOR(M) preconditioner
	protected:
			void solve(int *C,int *R,double *A,double *X,double *b);	//BiCG-stab algorithm
	public:
			BICGSTAB(int n,int max); ~BICGSTAB();
};
BICGSTAB::BICGSTAB(int n1,int max)
{
	N=n1; MAX=max;
	r0=new double[N];
	r=new double[N];
	p=new double[N];
	s=new double[N];
	m=new double[N];
	n=new double[N];
	y=new double[N];
	t=new double[N];
	cout<<"BICGSTAB: MEMORY ALLOCATED"<<endl;
}
BICGSTAB::~BICGSTAB()
{
	delete[] r0;
	delete[] r;
	delete[] p;
	delete[] s;
	delete[] m;
	delete[] n;
	delete[] y;
	delete[] t;
	cout<<"BICGSTAB: MEMORY RELEASED"<<endl;
}
void BICGSTAB::solve(int *C,int *R,double *A,double *X,double *b)
{
	int M;	//iteration counter
	norm=0.0;	//reinitialization
	for(int i=0;i<N;i++)	//matrix vector multiplication
	{
		sum=0.0;
		for(int d=R[i];d<R[i+1];d++) sum+=A[d]*X[C[d]];
		r[i]=b[i]-sum;	//initial residual vector
		r0[i]=r[i]; p[i]=r[i];
		norm+=r[i]*r[i];
	}
	gama=norm;
	norm=sqrt(norm);	//initial residual norm
	if(abs(norm)<=TOL) return;
	for(M=0;M<MAX;M++)	//convergence loop
	{
		SSOR(C,R,A,y,p,8);	//compute solution of 1st preconditioner
		norm=0.0;	//reinitialization
		for(int i=0;i<N;i++)	//matrix vector multiplication
		{
			m[i]=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++) m[i]+=A[d]*y[C[d]];
			norm+=m[i]*r0[i];
		}
		alpha=gama/norm;
		norm=0.0;	//reinitialization
		for(int d=0;d<N;d++)
		{
			s[d]=r[d]-alpha*m[d]; norm+=s[d]*s[d];
			t[d]=0.0;	//reinitialization
		}
		norm=sqrt(norm);
		if(norm<=TOL)	//convergence check1
		{
			for(int d=0;d<N;d++) X[d]+=alpha*y[d];
			break;
		}
		SSOR(C,R,A,t,s,8);	//compute solution of 2nd preconditioner
		sum=norm=0.0;	//reinitialization
		for(int i=0;i<N;i++)	//matrix vector multiplication
		{
			n[i]=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++) n[i]+=A[d]*t[C[d]];
			sum+=n[i]*s[i]; norm+=n[i]*n[i];
		}
		omega=sum/norm;
		norm=0.0;	//reinitialization
		for(int d=0;d<N;d++)
		{
			X[d]+=alpha*y[d]+omega*t[d];
			r[d]=s[d]-omega*n[d];
			norm+=r[d]*r[d];
		}
		norm=sqrt(norm);
		//cout<<"BICGSTAB: "<<M+1<<" "<<norm<<endl;
		if(norm<=TOL) break;	//convergence check2
		sum=gama;	//sum is used to store gama_M
		gama=0.0;	//reinitialization
		for(int d=0;d<N;d++) gama+=r[d]*r0[d];	//calculation of gama_M+1
		if(abs(gama)<=TOL)	//restart condition
		{
			gama=0.0;	//reinitialization
			for(int d=0;d<N;d++)
			{
				r0[d]=r[d];
				p[d]=r[d];
				gama+=r[d]*r[d];
			}
		}
		else
		{
			beta=alpha*gama/(omega*sum);
			for(int d=0;d<N;d++)
			{
				p[d]=r[d]+beta*(p[d]-omega*m[d]);
				y[d]=0.0;	//reinitialization
			}
		}
	}
	//cout<<"BICGSTAB: "<<M+1<<endl;
}
void BICGSTAB::SSOR(int *C,int *R,double *A,double *xp,double *bp,int M)
{
	double diag,sum;	//diagonal element of the matrix
	for(int m=0;m<M;m++)
	{
		for(int i=0;i<N;i++)	//matrix vector multiplication, forward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++)
			{
				if(C[d]==i) diag=A[d];
				else sum+=A[d]*xp[C[d]];
			}
			xp[i]=(bp[i]-sum)/diag;
		}
		for(int i=N-1;i>=0;i--)	//matrix vector multiplication, backward sweep
		{
			sum=0.0;	//reinitialization
			for(int d=R[i];d<R[i+1];d++)
			{
				if(C[d]==i) diag=A[d];
				else sum+=A[d]*xp[C[d]];
			}
			xp[i]=(bp[i]-sum)/diag;
		}
	}
}
