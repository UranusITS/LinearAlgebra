#pragma once
#include<memory>
#include<vector>

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
	Vector(Vector &&)=default;
	Vector(const unsigned int &,const double[]);
	Vector(const unsigned int &);
	~Vector();
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
	bool iszero() const;
	double mod() const;
	void unitize();
};

class VectorGroup
{
private:
	unsigned int height;
	std::vector<Vector>value;
public:
	VectorGroup();
	VectorGroup(const VectorGroup &);
	VectorGroup(VectorGroup &&)=default;
	VectorGroup(const unsigned int &,const Vector[]);
	VectorGroup(const std::vector<Vector> &);
	~VectorGroup();
	VectorGroup operator=(const VectorGroup &);
	double& operator()(const unsigned int &,const unsigned int &) const;
	friend std::ostream &operator<<(std::ostream &,const VectorGroup &);
	unsigned int size() const;
	void clear();
	bool empty() const;
	unsigned int get_height() const;
	void add(Vector &);
	bool orthogonalize();
	void unitize();
	bool schmidt();
};

#endif
