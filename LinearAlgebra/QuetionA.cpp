#include<iostream>
#include<vector>
#include"QuetionA.h"
long long Pow(long long x,int y)
{
    long long tmp=x,re=1;
    while(y)
    {
        if(y&1) re*=tmp;
        y>>=1,tmp*=tmp;
    }
    return re;
}
void QuestionA::init_matrix(Matrix &A)
{
    unsigned int cnt=0;
    double a[105]={0};
    std::vector<long long>V[10];
    V[0].push_back(5);
    V[0].push_back(2);
    V[0].push_back(0);
    V[0].push_back(0);
    V[0].push_back(2);
    V[0].push_back(1);
    V[0].push_back(9);
    V[0].push_back(1);
    V[0].push_back(0);
    V[0].push_back(3);
    V[0].push_back(0);
    V[0].push_back(8);
    cnt=12;
    for(unsigned int i=1;i<10&&cnt<100;i++)
    {
        for(unsigned int j=(i&1)?1:0;j<V[i-1].size();j+=2)
        {
            long long tmp=Pow(V[i-1][j],i+1);
            if(!tmp)
            {
                V[i].push_back(tmp);
                continue;
            }
            int sz=0;
            while(tmp)
            {
                V[i].push_back(tmp%10);
                tmp/=10,sz++;
            }
            for(unsigned int head=V[i].size()-sz,tail=V[i].size()-1;head<tail;head++,tail--)
                std::swap(V[i][head],V[i][tail]);
        }
        cnt+=V[i].size();
    }
    for(unsigned int i=0,pos=0;i<10&&pos<100;i++)
        for(unsigned int j=0;j<V[i].size()&&pos<100;j++)
            a[pos++]=V[i][j];
    A=Matrix(10,10,a);
}
int QuestionA::main()
{
    Matrix A;
    init_matrix(A);
    std::cout<<"A="<<std::endl<<A<<std::endl;
    std::cout<<"det(A)= "<<A.determinant()<<std::endl<<std::endl;
    std::cout<<"the reduced row-echelon matrix of A is:"<<std::endl<<A.reduced_row_echelon()<<std::endl;
    return 0;
}
