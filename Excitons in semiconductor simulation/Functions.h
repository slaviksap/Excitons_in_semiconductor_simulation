#pragma once
#include<random>
#include<iostream>
#include<cstdlib>
#include<algorithm>
#include<ctime>
#include"vector3.h"
using namespace std;
mt19937 random(time(0));
uniform_real_distribution<> gamma(0, 1);

const double Pi = 3.1415926535897;

double uniDistrib()
{
	return gamma(random);
}

double expDistrib(double a)
{
	return -a * log(uniDistrib());
}

double norm(Vector3 u)
{
	return sqrt(pow(u.getX(), 2) + pow(u.getY(), 2) + pow(u.getZ(), 2));
}

Vector3 randUnitVector()
{
	double cosPhi = 1 - 2 * uniDistrib();
	double sinPhi = sqrt(1 - cosPhi * cosPhi);
	double psi = 2 * Pi * uniDistrib();
	Vector3 v;
	v.setX(cos(psi) * sinPhi);
	v.setY(sin(psi) * sinPhi);
	v.setZ(cosPhi);
	return v;
}

int sgn(double x)
{
	if (x > 0)
		return 1;
	if (x < 0)
		return -1;
	return 0;
}

double min(double x1, double x2)
{
	if (x1 < x2)
		return x1;
	return x2;
}

double min(double x1, double x2, double x3)
{
	if (x1 < x2)
	{
		if (x1 < x3)
			return x1;
		return x3;
	}
	if (x2 < x3)
		return x2;
	return x3;
}

double min(double x1, double x2, double x3, double x4)
{
	return min(min(x1, x2), min(x3, x4));
}

double max(vector<double> v)
{
	double max = v[0];
	for (auto x : v)
	{
		if (x > max)
			max = x;
	}
	return max;
}

bool is_inside_sphere(Vector3 u, Vector3 c, double r)
{
	if (norm(u - c) <= r)
		return true;
	return false;
}

double vectorsAzimut(Vector3& u)
{
	if (u.getX() >= 0)
		return atan(u.getY() / u.getX());
	if (u.getY() >= 0)
		return Pi + atan(u.getY() / u.getX());
	return atan(u.getY() / u.getX()) - Pi;
}

double squared(double x)
{
	return x * x;
}

double cylinderDistance(Vector3& v, double cylX, double cylY, double cylR)
{
	return sqrt(squared(cylX - v.getX()) + squared(cylY - v.getY())) - cylR;
}

bool isInsideCylinder(Vector3& v, double cylX, double cylY, double cylR)
{
	if (sqrt(squared(cylX - v.getX()) + squared(cylY - v.getY())) < cylR)
		return true;
	return false;
}