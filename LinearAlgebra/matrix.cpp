#include<iostream>
#include<iomanip>
#include<cmath>
#include<climits>
#include"matrix.h"
#include"linearvector.h"
Matrix::Matrix()
{
	this->row=this->column=0;
	this->value.reset();
}
Matrix::Matrix(const Matrix &A)
{
	this->clear();
	this->row=A.row;
	this->column=A.column;
	this->value=std::make_unique<double[]>(this->row*this->column);
	for(unsigned int i=0;i<this->row*this->column;i++) this->value[i]=A.value[i];
}
Matrix::Matrix(const unsigned int &row,const unsigned int &column,const double value[])
{
	this->clear();
	this->row=row,this->column=column;
	if(row&&column)
	{
		this->value=std::make_unique<double[]>(row*column);
		for(unsigned int i=0;i<row*column;i++)
			this->value[i]=value[i];
	}
	else
		std::cerr<<"[ERROR]Can't generate an empty matrix!"<<std::endl;
}
Matrix::Matrix(const unsigned int &row,const unsigned int &column)
{
	this->clear();
	this->row=row,this->column=column;
	if(row&&column)
	{
		this->value=std::make_unique<double[]>(row*column);
		for(unsigned int i=0;i<row*column;i++)
			this->value[i]=0;
	}
	else
		std::cerr<<"[ERROR]Can't generate an empty matrix!"<<std::endl;
}
Matrix::Matrix(const VectorGroup &VG)
{
	this->clear();
	this->row=VG.get_height();
	this->column=VG.size();
	for(unsigned int i=0;i<VG.get_height();i++)
		for(unsigned int j=0;j<VG.size();j++)
			(*this)(i,j)=VG(i,j);
}
double& Matrix::operator()(const unsigned int &row,const unsigned int &column) const
{
	return this->value[row*this->column+column];
}
double& Matrix::at(const unsigned int &row,const unsigned int &column) const
{
	if(row>=this->row||column>=this->column)
	{
		throw std::out_of_range("The matrix is not big enough.");
	}
	return this->value[row*this->column+column];
}
Matrix Matrix::operator=(const Matrix &A)
{
	this->clear();
	this->row=A.row;
	this->column=A.column;
	this->value=std::make_unique<double[]>(this->row*this->column);
	for(unsigned int i=0;i<this->row*this->column;i++) this->value[i]=A.value[i];
	return A;
}
Matrix Matrix::row_echelon() const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	re=*this;
	for(unsigned int i=0;i<re.get_column();i++)
	{
		unsigned int mx=i;
		for(unsigned int j=i+1;j<re.get_row();j++)
			if(fabs(re(j,i))>fabs(re(mx,i)))
				mx=j;
		if(fabs(re(mx,i))<1e-4) continue;
		if(mx!=i)
		{
			for(unsigned int j=i;j<re.get_column();j++)
				std::swap(re(i,j),re(mx,j));
		}
		for(unsigned int j=i+1;j<re.get_row();j++)
		{
			double div=re(j,i)/re(i,i);
			for(unsigned int k=i;k<re.get_column();k++)
				re(j,k)=re(j,k)-re(i,k)*div;
		}
	}
	return re;
}
Matrix Matrix::reduced_row_echelon() const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	re=this->row_echelon();
	for(unsigned int i=0;i<re.get_row();i++)
	{
		unsigned int main_element=-1;
		for(unsigned int j=0;j<re.get_column();j++)
			if(fabs(re(i,j))>1e-4)
			{
				main_element=j;
				break;
			}
		if(main_element==-1) continue;
		double div=re(i,main_element);
		for(unsigned int j=main_element;j<re.get_column();j++)
			re(i,j)=re(i,j)/div;
	}
	for(unsigned int i=re.get_row()-1;i>0;i--)
	{
		unsigned int main_element=-1;
		for(unsigned int j=0;j<re.get_column();j++)
			if(fabs(re(i,j))>1e-4)
			{
				main_element=j;
				break;
			}
		if(main_element==-1) continue;
		for(unsigned int j=0;j<i;j++)
			if(fabs(re(j,main_element))>1e-4)
			{
				double div=re(j,main_element)/re(i,main_element);
				for(unsigned int k=main_element;k<re.get_column();k++)
					re(j,k)=re(j,k)-div*re(i,k);
			}
	}
	return re;
}
std::ostream &operator<<(std::ostream &output,const Matrix &A)
{
	if(!A.empty())
		for(unsigned int i=0;i<A.row;i++)
		{
			for(unsigned int j=0;j<A.column;j++)
				if(fabs(A(i,j))<1e-4) output<<std::setw(10)<<std::fixed<<std::setprecision(2)<<0.0<<' ';
				else output<<std::setw(10)<<std::fixed<<std::setprecision(2)<<A(i,j)<<' ';
			output<<std::endl;
		}
	else output<<"[WARNING]The matrix is empty."<<std::endl;
	return output;
}
std::istream &operator>>(std::istream &input,Matrix &A)
{
	A.clear();
	unsigned int row,column;
	input>>row>>column;
	if(row&&column)
	{
		A.set_size(row,column);
		double tmp;
		for(unsigned int i=0;i<row;i++)
			for(unsigned int j=0;j<column;j++)
			{
				input>>tmp;
				A(i,j)=tmp;
			}
	}
	else std::cerr<<"[ERROR]Can't input the matrix: matrix of 0*0 is not supported!"<<std::endl;
	return input;
}
Matrix operator+(const Matrix &A,const Matrix &B)
{
	Matrix re;
	if(A.get_row()!=B.get_row()||A.get_column()!=B.get_column()) std::cerr<<"[ERROR]The two matrices are different in size!"<<std::endl;
	else
	{
		re.set_size(A.get_row(),A.get_column());
		for(unsigned int i=0;i<re.row;i++)
			for(unsigned int j=0;j<re.column;j++)
				re(i,j)=A(i,j)+B(i,j);
	}
	return re;
}
Matrix operator-(const Matrix &A,const Matrix &B)
{
	Matrix re;
	if(A.get_row()!=B.get_row()||A.get_column()!=B.get_column()) std::cerr<<"[ERROR]The two matrices are different in size!"<<std::endl;
	else
	{
		re.set_size(A.get_row(),A.get_column());
		for(unsigned int i=0;i<re.row;i++)
			for(unsigned int j=0;j<re.column;j++)
				re(i,j)=A(i,j)-B(i,j);
	}
	return re;
}
Matrix operator*(const double &k,const Matrix &A)
{
	Matrix re;
	if(A.empty()) std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
	else
	{
		re.set_size(A.get_row(),A.get_column());
		for(unsigned int i=0;i<re.row;i++)
			for(unsigned int j=0;j<re.column;j++)
				re(i,j)=A(i,j)*k;
	}
	return re;
}
Matrix operator*(const Matrix &A,const Matrix &B)
{
	Matrix re;
	if(A.get_column()!=B.get_row()) std::cerr<<"[ERROR]Illegal matrix multiplication!"<<std::endl;
	else
	{
		re.set_size(A.get_row(),B.get_column());
		for(unsigned int i=0;i<A.get_row();i++)
			for(unsigned int j=0;j<B.get_column();j++)
				for(unsigned int k=0;k<A.get_column();k++)
					re(i,j)+=A(i,k)*B(k,j);
	}
	return re;
}
Matrix identity_matrix(const unsigned int &x)
{
	Matrix re;
	re.set_size(x,x);
	for(unsigned int i=0;i<x;i++) re(i,i)=1;
	return re;
}
unsigned int Matrix::get_row() const
{
	return this->row;
}
unsigned int Matrix::get_column() const
{
	return this->column;
}
unsigned int Matrix::size() const
{
	return this->row*this->column;
}
void Matrix::clear()
{
	if(!this->empty())
	{
		this->value.reset();
		this->row=this->column=0;
	}
}
bool Matrix::empty() const
{
	if(this->size()) return false;
	else return true;
}
bool Matrix::iszero() const
{
	for(unsigned int i=0;i<this->get_row();i++)
		for(unsigned int j=0;j<this->get_row();j++)
			if(fabs((*this)(i,j))>1e-4)
				return false;
	return true;
}
void Matrix::set_size(const unsigned int &row,const unsigned int &column)
{
	this->clear();
	this->row=row,this->column=column;
	if(row&&column) this->value=std::make_unique<double[]>(row*column);
	else std::cerr<<"[ERROR]Can't set size: empty matrix!"<<std::endl;
}
double Matrix::determinant() const
{
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return 0;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return 0;
	}
	Matrix tmp;
	tmp=*this;
	int flag=1;
	for(unsigned int i=0;i<tmp.get_column();i++)
	{
		unsigned int mx=i;
		for(unsigned int j=i+1;j<tmp.get_row();j++)
			if(fabs(tmp(j,i))>fabs(tmp(mx,i)))
				mx=j;
		if(fabs(tmp(mx,i))<1e-4) return 0;
		if(mx!=i)
		{
			flag*=-1;
			for(unsigned int j=i;j<tmp.get_column();j++)
				std::swap(tmp(i,j),tmp(mx,j));
		}
		for(unsigned int j=i+1;j<tmp.get_row();j++)
		{
			double div=tmp(j,i)/tmp(i,i);
			for(unsigned int k=i;k<tmp.get_column();k++)
				tmp(j,k)-=tmp(i,k)*div;
		}
	}
	double re=flag;
	for(unsigned int i=0;i<tmp.get_column();i++) re*=tmp(i,i);
	return re;
}
double Matrix::cominor(const unsigned int &row,const unsigned int &column) const
{
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return 0;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return 0;
	}
	Matrix tmp;
	tmp.set_size(this->get_row()-1,this->get_column()-1);
	for(unsigned int i=0;i<row;i++)
	{
		for(unsigned int j=0;j<column;j++)
			tmp(i,j)=(*this)(i,j);
		for(unsigned int j=column+1;j<this->get_column();j++)
			tmp(i,j-1)=(*this)(i,j);
	}
	for(unsigned int i=row+1;i<this->get_row();i++)
	{
		for(unsigned int j=0;j<column;j++)
			tmp(i-1,j)=(*this)(i,j);
		for(unsigned int j=column+1;j<this->get_column();j++)
			tmp(i-1,j-1)=(*this)(i,j);
	}
	return tmp.determinant();
}
double Matrix::trace() const
{
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return 0;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return 0;
	}
	double re=0;
	for(unsigned int i=0;i<this->get_column();i++)
		re+=(*this)(i,i);
	return re;
}
Matrix Matrix::power(unsigned int x) const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return re;
	}
	Matrix tmp;
	re=identity_matrix(this->get_column());
	tmp=identity_matrix(this->get_column());
	while(x)
	{
		if(x&1) re=re*tmp;
		tmp=tmp*tmp;
		x>>=1;
	}
	return re;
}
Matrix Matrix::transpose() const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	re.set_size(this->get_column(),this->get_row());
	for(unsigned int i=0;i<this->get_column();i++)
		for(unsigned int j=0;j<this->get_row();j++)
			re(i,j)=(*this)(j,i);
	return re;
}
Matrix Matrix::adjugate() const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return re;
	}
	re.set_size(this->get_row(),this->get_column());
	for(unsigned int i=0;i<this->get_column();i++)
		for(unsigned int j=0;j<this->get_row();j++)
			if((i+j)&1)
				re(i,j)=(*this)(i,j);
			else
				re(i,j)=this->cominor(i,j);
	return re;
}
Matrix Matrix::inverse() const
{
	Matrix re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	if(this->get_row()!=this->get_column())
	{
		std::cerr<<"[WARNING]The row of the matrix is not equal to the column."<<std::endl;
		return re;
	}
	if(this->determinant()==0)
	{
		std::cerr<<"[WARNING]The matrix has no inverse matrix."<<std::endl;
		return re;
	}
	re=(1/this->determinant())*this->adjugate();
	return re;
}
VectorGroup Matrix::solve_linear_equation() const
{
	VectorGroup re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	else if(this->iszero())
		return re;
	Matrix tmp=this->reduced_row_echelon();
	for(unsigned int i=0;i<tmp.get_column();i++)
	{
		bool flag=true;
		for(unsigned int j=0;j<tmp.get_row();j++)
			if(fabs(tmp(i,j))>1e-4)
			{
				flag=false;
				break;
			}
		if(flag)
		{
			Vector vtmp(tmp.get_column());
			vtmp[i]=1;
			re.add(vtmp);
		}
	}
	std::unique_ptr<bool[]>vis=std::make_unique<bool[]>(this->get_column());
	for(unsigned int i=0;i<tmp.get_column();i++)
		vis[i]=false;
	for(unsigned int i=0;i<tmp.get_row();i++)
		for(unsigned int j=0;j<tmp.get_column();j++)
			if(fabs(tmp(i,j))>1e-4)
			{
				vis[j]=true;
				break;
			}
	for(unsigned int i=0;i<tmp.get_column();i++)
	{
		if(vis[i]) continue;
		bool flag=false;
		unsigned int last=0;
		for(unsigned int j=0;j<tmp.get_row();j++)
			if(fabs(tmp(i,j))>1e-4)
			{
				flag=true;
				last=j;
			}
		if(flag)
		{
			Vector vtmp(this->get_column());
			for(unsigned int j=0;j<this->get_row();j++)
				if(fabs(tmp(i,j))>1e-4)
				{
					if(j==last)
						vtmp[j]=1;
					else
						vtmp[j]=tmp(j,i)/tmp(last,i);
				}
			re.add(vtmp);
		}
	}
	if(re.empty())
	{
		Vector vtmp(tmp.get_column());
		re.add(vtmp);
	}
	return re;
}
VectorGroup Matrix::max_linear_independent_group() const
{
	VectorGroup re;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return re;
	}
	else if(this->iszero())
		return re;
	Matrix tmp=this->transpose();
	Matrix red=tmp.reduced_row_echelon();
	tmp=tmp.transpose();
	for(unsigned int i=0;i<red.get_column();i++)
	{
		bool flag=false;
		for(unsigned int j=0;j<red.get_row();j++)
			if(fabs(red(i,j))>1e-4)
			{
				flag=true;
				break;
			}
		if(flag)
		{
			Vector vtmp(tmp.get_row());
			for(unsigned int j=0;j<tmp.get_row();j++)
				vtmp[j]=tmp(j,i);
			re.add(vtmp);
		}
	}
	return re;
}
