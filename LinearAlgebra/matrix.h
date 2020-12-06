#include<iostream>
#include<memory>
#pragma once

#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
private:
	unsigned int row,column;
	std::unique_ptr<double[]>value;
	unsigned int get_pos(const unsigned int &,const unsigned int &) const;
public:
	Matrix();
	Matrix(const Matrix &);
	Matrix(unsigned int,unsigned int,double[]);
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
	double get_element(const unsigned int &,const unsigned int &) const;
	void set_element(const unsigned int &,const unsigned int &,const double &);
	void set_size(const unsigned int &,const unsigned int &);
	double determinant() const;
	double cominor(const unsigned int &,const unsigned int &) const;
	double trace() const;
	Matrix power(unsigned int) const;
	Matrix transpose() const;
	Matrix adjugate() const;
	Matrix inverse() const;
	Matrix solve_linear_equation() const;
};

#endif
