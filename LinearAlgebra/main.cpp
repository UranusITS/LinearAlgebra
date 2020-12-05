#include<iostream>
#include"matrix.h"
using namespace std;
double a[9]={1,0,0,0,2,0,0,0,3};
int main()
{
	Matrix A=Matrix(3,3,a);
	cout<<A;
	cout<<A.determinant();
	return 0;
}
