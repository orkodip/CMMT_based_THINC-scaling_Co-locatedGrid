/*HEADER FILE

CONTENTS

SIGN FUNCTION
MAXIMUM OF TWO NUMBERS
MINIMUM OF TWO NUMBERS
CHECK WHETHER A NUMBER IS GREATER THAN ZERO OR NOT
MAXIMUM OF THREE NUMBERS
MINIMUM OF THREE NUMBERS
VECTOR CLASS (USES THE CONCEPT OF OPERATOR OVERLOADING)
MINMOD FUNCTION OF TWO NUMBERS

*/

#define SGN(A) ((A)>=0.0 ? 1.0 : -1.0)

double MAX2 (double a, double b)
{
	if (a>=b)
		return a;
	else
		return b;
}

double MIN2 (double a, double b)
{
	if (a<=b)
		return a;
	else
		return b;
}

double CHK(double a)
{
	if (a>=0.0)
		return 1.0;
	else
		return 0.0;
}

double MAX3(double a,double b,double c)
{
	if((a>=b)&&(a>=c)) return a;
	else if((b>=a)&&(b>=c)) return b;
	else return c;
}

double MIN3(double a,double b,double c)
{
	if((a<=b)&&(a<=c)) return a;
	else if((b<=a)&&(b<=c)) return b;
	else return c;
}

class VECTOR
{
	public:
			double x,y,z;	//x,y,z components of the vector
			VECTOR() {x=0.0; y=0.0; z=0.0;}	//null vector initialization
			VECTOR(double a,double b,double c=0.0) {x=a; y=b; z=c;}	//initialize the vector with prescribed values
			double mag() {return sqrt(pow(x,2.0)+pow(y,2.0)+pow(z,2.0));}	//return magnitude of the vector
			VECTOR unit();	//return unit vector
			VECTOR operator+(VECTOR V);	//vector addition
			VECTOR operator-();	//change the direction of the vector
			VECTOR operator-(VECTOR V);	//vector substraction
			friend double dot(VECTOR V1,VECTOR V2);	//dot product
			friend VECTOR operator*(double m,VECTOR V);	//scalar multiplication (type 1)
			friend VECTOR operator*(VECTOR V,double m);	//scalar multiplication (type 2)
			friend VECTOR operator*(VECTOR V1,VECTOR V2);	//cross product
};
VECTOR VECTOR::unit()
{
	double M=mag();
	VECTOR uV(x/M,y/M,z/M);	//calculate the unit vector
	return uV;
}
VECTOR VECTOR::operator+(VECTOR V)
{
	VECTOR add;
	add.x=x+V.x;
	add.y=y+V.y;
	add.z=z+V.z;
	return add;
}
VECTOR VECTOR::operator-()
{
	VECTOR s(-x,-y,-z);
	return s;
}
VECTOR VECTOR::operator-(VECTOR V)
{
	VECTOR sub;
	sub.x=x-V.x;
	sub.y=y-V.y;
	sub.z=z-V.z;
	return sub;
}
double dot(VECTOR V1,VECTOR V2) {return (V1.x*V2.x+V1.y*V2.y+V1.z*V2.z);}
VECTOR operator*(double m,VECTOR V)
{
	VECTOR mul(m*V.x,m*V.y,m*V.z);
	return mul;
}
VECTOR operator*(VECTOR V,double m)
{
	VECTOR mul(m*V.x,m*V.y,m*V.z);
	return mul;
}
VECTOR operator*(VECTOR V1,VECTOR V2)
{
	VECTOR cp;
	cp.x=V1.y*V2.z-V2.y*V1.z;	//calculate the cross product
	cp.y=-(V1.x*V2.z-V2.x*V1.z);
	cp.z=V1.x*V2.y-V2.x*V1.y;
	return cp;
}

double MINMOD(double a,double b)
{
	if(abs(a)<=abs(b)) return a;
	else return b;
}
