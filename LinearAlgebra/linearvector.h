#include<memory>
#pragma once

#ifndef LINEARVECTOR_H
#define LINEARVECTOR_H

class Vector
{
private:
	unsigned int length;
	std::unique_ptr<double[]>value;
public:
	Vector();
	Vector(const Vector &);
	Vector(const unsigned int &,const double[]);
	Vector operator=(const Vector &);
	double& operator[](const unsigned int &) const;
	friend std::ostream &operator<<(std::ostream &,const Vector &);
	friend std::istream &operator>>(std::istream &,Vector &);
	friend Vector operator+(const Vector &,const Vector &);
	friend Vector operator-(const Vector &,const Vector &);
	friend Vector operator*(const double &,const Vector &);
	friend double operator*(const Vector &,const Vector &);
	friend double included_angle(const Vector &,const Vector &);
	friend double distance(const Vector &,const Vector &);
	unsigned int size() const;
	void clear();
	bool empty() const;
	double mod() const;
	void unitize();
};

#endif
