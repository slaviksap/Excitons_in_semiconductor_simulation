#pragma once
class Vector3
{
private:
	double x;
	double y;
	double z;
public:
	Vector3()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	Vector3(double vx, double vy, double vz)
	{
		x = vx;
		y = vy;
		z = vz;
	}
	friend const Vector3 operator+(const Vector3& left, const Vector3& right)
	{
		return Vector3(left.x + right.x, left.y + right.y, left.z + right.z);
	}
	friend const Vector3 operator-(const Vector3& left, const Vector3& right)
	{
		return Vector3(left.x - right.x, left.y - right.y, left.z - right.z);
	}
	friend const Vector3 operator-(const Vector3& right)
	{
		return Vector3(-right.x, -right.y, -right.z);
	}
	friend const Vector3 operator*(const Vector3& left, const double right)
	{
		return Vector3(left.x * right, left.y * right, left.z * right);
	}
	friend const Vector3 operator*(const double left, const Vector3& right)
	{
		return Vector3(right.x * left, right.y * left, right.z * left);
	}
	friend const Vector3 operator/(const Vector3& left, const double right)
	{
		return Vector3(left.x / right, left.y / right, left.z / right);
	}
	friend Vector3& operator+=(Vector3& left, const Vector3& right)
	{
		left = left + right;
		return left;
	}
	friend Vector3& operator-=(Vector3& left, const Vector3& right)
	{
		left = left - right;
		return left;
	}
	friend bool operator==(const Vector3& left, const Vector3& right)
	{
		return left.x == right.x && left.y == right.y && left.z == right.z;
	}
	double getX()
	{
		return x;
	}
	double getY()
	{
		return y;
	}
	double getZ()
	{
		return z;
	}
	void setX(double value)
	{
		x = value;
	}
	void setY(double value)
	{
		y = value;
	}
	void setZ(double value)
	{
		z = value;
	}
	friend double scalProd(Vector3& left, Vector3& right)
	{
		return left.getX() * right.getX() + left.getY() * right.getY() + left.getZ() * right.getZ();
	}
};
