#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include"QuetionB.h"
#include"matrix.h"
#include"linearvector.h"
void QuestionB::init_matrix(Matrix &A)
{
    A=Matrix(5,5);
    A(0,1)=-1;
    A(1,0)=1;
    A(2,2)=1;
    A(3,3)=1;
    A(4,4)=1;
}
void QuestionB::subquetion_1(const Matrix &A)
{
    std::cout<<"A="<<std::endl<<A<<std::endl;
}
void QuestionB::subquetion_234(const Matrix &A)
{
    std::vector<double>V;
    std::vector<VectorGroup>VG;
    unsigned int cnt=0,pos=0;
    double lim=0;
    for(unsigned int i=0;i<A.get_row();i++)
        for(unsigned int j=0;j<A.get_column();j++)
            lim+=fabs(A(i,j));
    std::mt19937 gen(time(NULL));
    std::uniform_real_distribution<double> dis(-lim-1,lim+1);
    for(unsigned int i=0;i<2000&&cnt<5;i++)
    {
        double l=dis(gen),r=dis(gen);
        if(l>r) std::swap(l,r);
        double dl=(A-l*identity_matrix(5)).determinant();
        double dr=(A-r*identity_matrix(5)).determinant();
        while(dl*dr>0)
        {
            l=dis(gen),r=dis(gen);
            if(l>r) std::swap(l,r);
            dl=(A-l*identity_matrix(5)).determinant();
            dr=(A-r*identity_matrix(5)).determinant();
        }
        while(r-l>1e-12)
        {
            if(fabs(dl)<1e-12)
            {
                r=l;
                break;
            }
            if(fabs(dr)<1e-12)
            {
                l=r;
                break;
            }
            double mid=(l+r)/2;
            dl=(A-l*identity_matrix(5)).determinant();
            dr=(A-r*identity_matrix(5)).determinant();
            double dm=(A-mid*identity_matrix(5)).determinant();
            if(fabs(dm)<1e-12)
            {
                l=r=mid;
                break;
            }
            if(dm*dl<0) r=mid;
            else l=mid;
        }
        bool flag=true;
        for(unsigned int i=0;i<V.size();i++)
        {
            if(fabs(V[i]-l)<1e-3)
            {
                flag=false;
                break;
            }
        }
        if(flag)
        {
            VectorGroup tmp=(A-l*identity_matrix(5)).solve_linear_equation();
            VG.push_back(tmp);
            cnt+=tmp.size();
            for(unsigned int j=0;j<tmp.size();j++)
                V.push_back(l);
        }
    }
    std::cout<<"f(x)=(x-1)^3*(x^2+1)"<<std::endl<<std::endl;
    if(!V.size()) std::cout<<"A has no characteristic value."<<std::endl<<std::endl;
    else if(V.size()==1) std::cout<<"The only characteristic value of A is:"<<std::endl<<V[0]<<std::endl<<std::endl;
    else
    {
        std::cout<<"The characteristic values of A are:"<<std::endl<<V[0];
        for(unsigned int i=1;i<V.size();i++)
            if(fabs(V[i]-V[0])>1e-3)
                std::cout<<','<<V[i];
        std::cout<<std::endl<<std::endl;
    }
    std::cout<<"The basic solution system of the equition ("<<V[0]<<"E-A)x=0:"<<std::endl<<VG[0]<<std::endl;
    for(unsigned int i=0;i<V.size();i++)
        if(fabs(V[i]-V[0])>1e-3)
            std::cout<<"The basic solution system of the equition ("<<V[i]<<"E-A)x=0:"<<std::endl<<VG[i]<<std::endl;
}
/*
void QuestionB::subquetion_3(const std::vector<double> &V)
{
    if(!V.size()) std::cout<<"A has no characteristic value."<<std::endl<<std::endl;
    else if(V.size()==1) std::cout<<"The only characteristic value of A is:"<<std::endl<<V[0]<<std::endl<<std::endl;
    else
    {
        std::cout<<"The characteristic values of A are:"<<std::endl<<V[0];
        for(unsigned int i=1;i<V.size();i++)
            std::cout<<','<<V[i];
        std::cout<<std::endl<<std::endl;
    }
}
void QuestionB::subquetion_4(const std::vector<double> &V,const VectorGroup &VG)
{
    for(unsigned int i=0;i<V.size();i++)
    {
        std::cout<<"The basic solution system of the equition ("<<V[i]<<"E-A)x=0:"<<std::endl<<result<<std::endl;
    }
}
*/
int QuestionB::main()
{
    Matrix A;
    init_matrix(A);
    subquetion_1(A);
    subquetion_234(A);
    return 0;
}
