#include<iostream>
#include<iomanip>
#include<cmath>
#include<climits>
#include<memory>
#include"matrix.h"
unsigned int Matrix::get_pos(const unsigned int &row,const unsigned int &column) const
{
	if(!this->empty()) return row*this->column+column;
	else return 0;
}
Matrix::Matrix()
{
	this->row=this->column=0;
	this->value.reset();
}
Matrix::Matrix(unsigned int row,unsigned int column,double value[])
{
	this->clear();
	this->row=row,this->column=column;
	if(row&&column)
	{
		this->value=std::make_shared<double[]>(row*column,0);
		for(unsigned int i=0;i<row*column;i++)
			this->value[i]=value[i];
	}
	else
		std::cerr<<"[ERROR]Can't get the element: empty matrix!"<<std::endl;
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
			if(fabs(re.get_element(j,i))>fabs(re.get_element(mx,i)))
				mx=j;
		if(fabs(re.get_element(mx,i))<1e-6) continue;
		if(mx!=i)
		{
			for(unsigned int j=i;j<re.get_column();j++)
				std::swap(re.value[re.get_pos(i,j)],re.value[re.get_pos(mx,i)]);
		}
		for(unsigned int j=i+1;j<re.get_row();j++)
		{
			double div=re.get_element(j,i)/re.get_element(i,i);
			for(unsigned int k=i;k<re.get_column();k++)
				re.set_element(j,k,re.get_element(j,k)-get_element(i,k)*div);
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
			if(fabs(re.get_element(i,j))>1e-6)
			{
				main_element=j;
				break;
			}
		if(main_element==-1) continue;
		for(unsigned int j=re.get_column()-1;j>=main_element;j--)
			re.set_element(i,j,re.get_element(i,j)/re.get_element(i,main_element));
	}
	for(unsigned int i=re.get_row()-1;i!=UINT_MAX;i--)
	{
		unsigned int main_element=-1;
		for(unsigned int j=0;j<re.get_column();j++)
			if(fabs(re.get_element(i,j))>1e-6)
			{
				main_element=j;
				break;
			}
		if(main_element==-1) continue;
		for(unsigned int j=0;j<i;j++)
			if(fabs(re.get_element(j,main_element))>1e-6)
				for(unsigned int k=re.get_column()-1;k>=main_element;k--)
					re.set_element(j,k,re.get_element(j,k)-re.get_element(j,main_element)/get_element(i,main_element)*re.get_element(j,k));
	}
	return re;
}
std::ostream& operator<<(std::ostream &output,const Matrix &A)
{
	if(!A.empty())
		for(unsigned int i=0;i<A.row;i++)
		{
			for(unsigned int j=0;j<A.column;j++)
				output<<std::fixed<<std::setprecision(2)<<A.get_element(i,j)<<' ';
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
				A.set_element(i,j,tmp);
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
				re.set_element(i,j,A.get_element(i,j)+B.get_element(i,j));
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
				re.set_element(i,j,A.get_element(i,j)-B.get_element(i,j));
	}
	return re;
}
Matrix operator*(const double &k,const Matrix &A)
{
	Matrix re;
	if(A.empty()) std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
	else
	{
		re=A;
		for(unsigned int i=0;i<re.row;i++)
			for(unsigned int j=0;j<re.column;j++)
				re.set_element(i,j,re.get_element(i,j)*k);
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
					re.set_element(i,j,re.get_element(i,j)+A.get_element(i,k)*B.get_element(k,j));
	}
	return re;
}
Matrix identity_matrix(const unsigned int &x)
{
	Matrix re;
	re.set_size(x,x);
	for(unsigned int i=0;i<x;i++) re.set_element(i,i,1);
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
double Matrix::get_element(const unsigned int &row,const unsigned int &column) const
{
	if(row>this->row||column>this->column)
	{
		std::cerr<<"[ERROR]Can't get the element: matrix not big enough!"<<std::endl;
		return 0;
	}
	else if(this->empty())
	{
		std::cerr<<"[ERROR]Can't get the element: empty matrix!"<<std::endl;
		return 0;
	}
	return this->value[this->get_pos(row,column)];
}
void Matrix::set_element(const unsigned int &row,const unsigned int &column,const double &new_value)
{
	if(row>this->row||column>this->column)
	{
		std::cerr<<"[ERROR]Can't set the element: matrix not big enough!"<<std::endl;
		return ;
	}
	else if(this->empty())
	{
		std::cerr<<"[ERROR]Can't set the element: empty matrix!"<<std::endl;
		return ;
	}
	this->value[this->get_pos(row,column)]=new_value;
}
void Matrix::set_size(const unsigned int &row,const unsigned int &column)
{
	this->clear();
	this->row=row,this->column=column;
	if(row&&column) this->value=std::make_shared<double[]>(row*column,0);
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
			if(fabs(tmp.get_element(j,i))>fabs(tmp.get_element(mx,i)))
				mx=j;
		if(fabs(tmp.get_element(mx,i))<1e-6) return 0;
		if(mx!=i)
		{
			flag*=-1;
			for(unsigned int j=i;j<tmp.get_column();j++)
				std::swap(tmp.value[tmp.get_pos(i,j)],tmp.value[tmp.get_pos(mx,j)]);
		}
		for(unsigned int j=i+1;j<tmp.get_row();j++)
		{
			double div=tmp.get_element(j,i)/tmp.get_element(i,i);
			for(unsigned int k=i;k<tmp.get_column();k++)
				tmp.set_element(j,k,tmp.get_element(j,k)-get_element(i,k)*div);
		}
	}
	double re=flag;
	for(unsigned int i=0;i<tmp.get_column();i++) re*=tmp.get_element(i,i);
	tmp.clear();
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
			tmp.set_element(i,j,this->get_element(i,j));
		for(unsigned int j=column+1;j<this->get_column();j++)
			tmp.set_element(i,j-1,this->get_element(i,j));
	}
	for(unsigned int i=row+1;i<this->get_row();i++)
	{
		for(unsigned int j=0;j<column;j++)
			tmp.set_element(i-1,j,this->get_element(i,j));
		for(unsigned int j=column+1;j<this->get_column();j++)
			tmp.set_element(i-1,j-1,this->get_element(i,j));
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
		re+=this->get_element(i,i);
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
			re.set_element(i,j,this->get_element(j,i));
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
				re.set_element(i,j,-this->cominor(i,j));
			else
				re.set_element(i,j,this->cominor(i,j));
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
Matrix Matrix::solve_linear_equation() const
{
	Matrix tmp;
	if(this->empty())
	{
		std::cerr<<"[WARNING]The matrix is empty."<<std::endl;
		return tmp;
	}
	tmp=this->reduced_row_echelon();
	for(unsigned int i=tmp.get_row()-1;i!=UINT_MAX;i--)
	{
		unsigned int cnt=0;
		for(unsigned int j=0;j<tmp.get_column();j++)
			if(fabs(tmp.get_element(i,j))>1e-6)
				cnt++;
		if(cnt>1)
		{
			break;
		}
		else if(i==0)
		{
			Matrix re;
			re.set_size(tmp.get_row(),1);
			return re;
		}
	}
	///TODO
	return tmp;
}
