#include<iostream>
#include<iomanip>
#include<cmath>
#include "linearvector.h"
Vector::Vector()
{
	this->length=0;
	value.reset();
}
Vector::Vector(const Vector &A)
{
	this->clear();
	this->length=A.length;
	value=std::make_unique<double[]>(this->length);
	for(unsigned int i=0;i<this->length;i++)
		this->value[i]=A[i];
}
Vector::Vector(const unsigned int &length,const double value[])
{
	this->clear();
	this->length=length;
	if(length)
	{
		this->value=std::make_unique<double[]>(length);
		for(unsigned int i=0;i<length;i++)
			this->value[i]=value[i];
	}
	else
		std::cerr<<"[ERROR]Can't generate an empty vector!"<<std::endl;
}
Vector Vector::operator=(const Vector &A)
{
	this->clear();
	this->length=A.length;
	this->value=std::make_unique<double[]>(this->length);
	for(unsigned int i=0;i<this->length;i++) this->value[i]=A.value[i];
	return A;
}
double& Vector::operator[](const unsigned int &x) const
{
	if(x<this->length) return this->value[x];
	else
	{
		std::cerr<<"[ERROR]Can't get the element: vector not big enough!"<<std::endl;
		return this->value[0];
	}
}
std::ostream &operator<<(std::ostream &output,const Vector &A)
{
	if(!A.empty())
		for(unsigned int i=0;i<A.length;i++)
		{
			for(unsigned int j=0;j<A.length;j++)
				if(fabs(A.value[i])<1e-6) output<<std::setw(10)<<std::fixed<<std::setprecision(2)<<0.0<<' ';
				else output<<std::setw(10)<<std::fixed<<std::setprecision(2)<<A[i]<<' ';
			output<<std::endl;
		}
	else output<<"[WARNING]The vector is empty."<<std::endl;
	return output;
}
std::istream &operator>>(std::istream &input,Vector &A)
{
	A.clear();
	unsigned int length;
	input>>length;
	if(length)
	{
		A.length=length;
		double tmp;
		for(unsigned int i=0;i<length;i++)
		{
			input>>tmp;
			A[i]=tmp;
		}
	}
	else std::cerr<<"[ERROR]Can't input the matrix: vector of size zero is not supported!"<<std::endl;
	return input;
}
Vector operator+(const Vector &A,const Vector &B)
{
	Vector re;
	if(A.length!=B.length) std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
	else
	{
		re.length=A.length;
		for(unsigned int i=0;i<re.length;i++)
			re[i]=A[i]+B[i];
	}
	return re;
}
Vector operator-(const Vector &A,const Vector &B)
{
	Vector re;
	if(A.length!=B.length) std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
	else
	{
		re.length=A.length;
		for(unsigned int i=0;i<re.length;i++)
			re[i]=A[i]-B[i];
	}
	return re;
}
Vector operator*(const double &k,const Vector &A)
{
	Vector re;
	if(A.empty()) std::cerr<<"[WARNING]The vector is empty."<<std::endl;
	else
	{
		re=A;
		for(unsigned int i=0;i<re.length;i++)
			re[i]*=k;
	}
	return re;
}
double included_angle(const Vector &A,const Vector &B)
{
	if(A.length!=B.length)
	{
		std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
		return 0;
	}
	else return acos(A*B/A.mod()/B.mod());
}
double distance(const Vector &A,const Vector &B)
{
	if(A.length!=B.length)
	{
		std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
		return 0;
	}
	else return (A-B).mod();
}
double operator*(const Vector &A,const Vector &B)
{
	double re=0;
	if(A.length!=B.length) std::cerr<<"[ERROR]Illegal vector scalar product!"<<std::endl;
	else
		for(unsigned int i=0;i<A.length;i++)
			re+=A[i]*B[i];
	return re;
}
unsigned int Vector::size() const
{
	return this->length;
}
void Vector::clear()
{
	if(!this->empty())
	{
		this->length=0;
		this->value.reset();
	}
}
bool Vector::empty() const
{
	if(this->size()) return false;
	else return true;
}
double Vector::mod() const
{
	return sqrt((*this)*(*this));
}
void Vector::unitize()
{
	double div=this->mod();
	for(unsigned int i=0;i<this->length;i++)
		this->value[i]/=div;
}
