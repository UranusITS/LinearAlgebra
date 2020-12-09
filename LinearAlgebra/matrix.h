#pragma once
#include<iostream>
#include<memory>
#include"linearvector.h"

#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
private:
	unsigned int row,column;
	std::unique_ptr<double[]>value;
public:
	Matrix();
	Matrix(const Matrix &);
	Matrix(Matrix &&)=default;
	Matrix(const unsigned int &,const unsigned int &,const double[]);
	Matrix(const unsigned int &,const unsigned int &);
	Matrix(const VectorGroup &);
	~Matrix()=default;
	double& operator()(const unsigned int &,const unsigned int &) const;
	double& at(const unsigned int &,const unsigned int &) const;
	Matrix operator=(const Matrix &);
	Matrix row_echelon() const;
	Matrix reduced_row_echelon() const;
	friend std::ostream& operator<<(std::ostream &,const Matrix &);
	friend std::istream& operator>>(std::istream &,Matrix &);
	friend Matrix operator+(const Matrix &,const Matrix &);
	friend Matrix operator-(const Matrix &,const Matrix &);
	friend Matrix operator*(const double &,const Matrix &);
	friend Matrix operator*(const Matrix &,const Matrix &);
	friend Matrix identity_matrix(const unsigned int &);
	unsigned int get_row() const;
	unsigned int get_column() const;
	unsigned int size() const;
	void clear();
	bool empty() const;
	bool iszero() const;
	void set_size(const unsigned int &,const unsigned int &);
	double determinant() const;
	double cominor(const unsigned int &,const unsigned int &) const;
	double trace() const;
	Matrix power(unsigned int) const;
	Matrix transpose() const;
	Matrix adjugate() const;
	Matrix inverse() const;
	VectorGroup solve_linear_equation() const;
	VectorGroup max_linear_independent_group() const;
};

#endif
