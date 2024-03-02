//POST-PROCESSING CODE
class POST:public NS
{
	double Xl[161],Yl[161],T[161];	//non-dimensional X, Y locations and time
	public:
			void calc_xloc(int NITER);
			void calc_yloc(int NITER);
			void write();
};
void POST::calc_xloc(int NITER)
{
	int i_f;
	int n=NITER/20;
	for(int i=2;i<=I;i++)
	{
		if((F[1][i-1]>=0.5)&&(F[1][i]<=0.5))	//cells in between which interface is present
		{ i_f=i; break; }
	}
	Xl[n]=CX[i_f-1]+(CX[i_f]-CX[i_f-1])/(F[1][i_f]-F[1][i_f-1])*(0.5-F[1][i_f-1]);
	Xl[n]/=0.05715;
	T[n]=NITER*dt*sqrt(9.81/0.05715);
}
void POST::calc_yloc(int NITER)
{
	int j_f;
	int n=NITER/20;
	for(int j=J-1;j>=1;j--)
	{
		if((F[j][1]>=0.5)&&(F[j+1][1]<=0.5))	//cells in between which interface is present
		{ j_f=j; break; }
	}
	Yl[n]=CY[j_f+1]+(CY[j_f]-CY[j_f+1])/(F[j_f][1]-F[j_f+1][1])*(0.5-F[j_f+1][1]);
	Yl[n]/=0.05715;
}
void POST::write()
{
	Xl[0]=1.0; Yl[0]=1.0; T[0]=0.0;	//initial location
	ofstream p_out("post.dat");
	p_out<<"TITLE = \"NON DIMENSIONAL TIME AND POSITION\""<<endl;
	p_out<<"VARIABLES = \"NT\",\"NX\",\"NY\""<<endl;
	p_out<<"ZONE I=161, DATAPACKING=POINT"<<endl;
	for(int j=0;j<=160;j++)
		p_out<<T[j]<<"\t"<<Xl[j]<<"\t"<<Yl[j]<<endl;
	p_out.close();
	cout<<"POST: FILE OUTPUT SUCCESSFULL."<<endl;
}
