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
Vector::Vector(const unsigned int &length)
{
	this->clear();
	this->length=length;
	if(length)
	{
		this->value=std::make_unique<double[]>(length);
		for(unsigned int i=0;i<length;i++)
			this->value[i]=0;
	}
	else
		std::cerr<<"[ERROR]Can't generate an empty vector!"<<std::endl;
}
Vector::~Vector()
{
	this->clear();
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
	{
		if(fabs(A[0])<1e-4) output<<'('<<std::setw(9)<<std::fixed<<std::setprecision(2)<<0.0;
		else output<<'('<<std::setw(9)<<std::fixed<<std::setprecision(2)<<A[0];
		for(unsigned int i=1;i<A.length;i++)
		{
			if(fabs(A[i])<1e-4) output<<','<<std::setw(10)<<std::fixed<<std::setprecision(2)<<0.0;
			else output<<','<<std::setw(10)<<std::fixed<<std::setprecision(2)<<A[i];
		}
		output<<')'<<std::endl;
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
	Vector re(A.length);
	if(A.length!=B.length) std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
	else
	{
		for(unsigned int i=0;i<re.length;i++)
			re[i]=A[i]+B[i];
	}
	return re;
}
Vector operator-(const Vector &A,const Vector &B)
{
	Vector re(A.length);
	if(A.length!=B.length) std::cerr<<"[ERROR]The two vectors are different in size!"<<std::endl;
	else
	{
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
bool Vector::iszero() const
{
	for(unsigned int i=0;i<this->size();i++)
		if(fabs((*this)[i])>1e-4) return false;
	return true;
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

VectorGroup::VectorGroup()
{
	this->height=0;
	this->value.clear();
	this->value.shrink_to_fit();
}
VectorGroup::VectorGroup(const VectorGroup &VG)
{
	this->height=0;
	this->value=VG.value;
}
VectorGroup::VectorGroup(const unsigned int &n,const Vector V[])
{
	this->clear();
	if(n)
	{
		unsigned int height=V[0].size();
		for(unsigned int i=1;i<n;i++)
			if(V[i].size()!=height)
				std::cerr<<"[ERROR]The vectors are not the same in size!"<<std::endl;
		this->height=height;
		for(unsigned int i=0;i<n;i++)
			this->value.push_back(V[i]);
	}
	else std::cerr<<"[ERROR]Can't generate an empty vector!"<<std::endl;
}
VectorGroup::VectorGroup(const std::vector<Vector> &V)
{
	this->clear();
	if(V.size())
	{
		unsigned int height=V[0].size();
		for(unsigned int i=1;i<V.size();i++)
			if(V[i].size()!=height)
				std::cerr<<"[ERROR]The vectors are not the same in size!"<<std::endl;
		this->height=height;
		this->value=V;
	}
	else std::cerr<<"[ERROR]Can't generate an empty vector!"<<std::endl;
}
VectorGroup::~VectorGroup()
{
	this->clear();
}
VectorGroup VectorGroup::operator=(const VectorGroup &VG)
{
	std::cout<<"FUUUCK"<<std::endl;
	this->clear();
	if(!VG.empty())
	{
		this->height=VG.value[0].size();
		this->value=VG.value;
	}
	return VG;
}
double &VectorGroup::operator()(const unsigned int &row,const unsigned int &column) const
{
	return this->value[row][column];
}
std::ostream &operator<<(std::ostream &output,const VectorGroup &VG)
{
	if(!VG.empty())
	{
		for(unsigned int i=0;i<VG.size();i++)
			output<<VG.value[i];
		output<<std::endl;
	}
	else output<<"[WARNING]The vector group is empty."<<std::endl;
	return output;
}
unsigned int VectorGroup::size() const
{
	return this->value.size();
}
void VectorGroup::clear()
{
	this->height=0;
	this->value.clear();
	this->value.shrink_to_fit();
}
bool VectorGroup::empty() const
{
	if(this->size()) return false;
	else return true;
}
unsigned int VectorGroup::get_height() const
{
	return this->height;
}
void VectorGroup::add(Vector &x)
{
	if(this->empty())
	{
		this->height=x.size();
		this->value.push_back(x);
	}
	else
	{
		if(this->height!=x.size())
			std::cerr<<"[ERROR]The vector added in has invalid size!"<<std::endl;
		else
			this->value.push_back(x);
	}
}
void VectorGroup::append(const VectorGroup &VG)
{
	if(VG.empty())
		return ;
	if(this->empty())
	{
		(*this)=VG;
		return ;
	}
	if(this->get_height()!=VG.get_height())
	{
		std::cerr<<"[ERROR]The VectorGroup appended in has invalid size!"<<std::endl;
		std::cerr<<this->get_height()<<' '<<VG.get_height()<<std::endl<<std::endl;
		return ;
	}
	Vector tmp(this->get_height());
	for(unsigned int i=0;i<VG.size();i++)
	{
		for(unsigned int j=0;j<VG.get_height();j++)
			tmp[j]=VG(i,j);
		this->add(tmp);
	}
}
bool VectorGroup::orthogonalize()
{
	if(this->empty()) return false;
	for(unsigned int i=1;i<this->value.size();i++)
	{
		Vector tmp=this->value[i];
		for(unsigned int j=0;j<i;j++)
			this->value[i]=this->value[i]-(tmp*this->value[j])/(this->value[j]*this->value[j])*this->value[j];
		if(this->value[i].iszero()) return false;
	}
	return true;
}
void VectorGroup::unitize()
{
	for(unsigned int i=0;i<this->value.size();i++)
		this->value[i].unitize();
}
bool VectorGroup::schmidt()
{
	if(!this->orthogonalize()) return false;
	this->unitize();
	return true;
}
