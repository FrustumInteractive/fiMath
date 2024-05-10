/* 3D Vector Class (double version) *
 * c Roger Dass 2011                *
 *                                  *
 *                                  */

#ifndef VEC3D_H
#define VEC3D_H

#include <math.h>

namespace FI {
namespace MATH {

// double version // TODO: make vec3 templated

class vec3d
{
public:

	inline vec3d()
	{
	}

	inline vec3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
	{
	}

	inline vec3d operator - () const
	{
		return vec3d(-x, -y, -z);
	}

	inline void operator -= (const vec3d &_v)
	{
		x -= _v.x;
		y -= _v.y;
		z -= _v.z;
	}

	inline void operator += (const vec3d &_v)
	{
		x += _v.x;
		y += _v.y;
		z += _v.z;
	}

	inline void operator *= (double _mul)
	{
		x *= _mul;
		y *= _mul;
		z *= _mul;
	}

	inline void operator *= (const vec3d &_v)
	{
		x *= _v.x;
		y *= _v.y;
		z *= _v.z;
	}

	inline void operator /= (double _div)
	{
		double mul = 1.0f / _div;
		x *= mul;
		y *= mul;
		z *= mul;
	}

	inline vec3d operator - (const vec3d &_v) const
	{
		return vec3d(x - _v.x, y - _v.y, z - _v.z);
	}

	inline vec3d operator + (const vec3d &_v) const
	{
		return vec3d(x + _v.x, y + _v.y, z + _v.z);
	}

	inline vec3d operator * (const vec3d &_v) const
	{
		return vec3d(x * _v.x, y * _v.y, z * _v.z);
	}

	inline vec3d operator * (double _m) const
	{
		return vec3d(x * _m, y * _m, z * _m);
	}

	inline vec3d operator / (const vec3d &_v) const
	{
		return vec3d(x / _v.x, y / _v.y, z / _v.z);
	}

	inline vec3d operator / (double _d) const
	{
		double m = 1.0f / _d;
		return vec3d(x * m, y * m, z * m);
	}

	inline vec3d operator | (const vec3d &_d) const
	{
		return vec3d(y * _d.z - z * _d.y,
					 z * _d.x - x * _d.z,
					 x * _d.y - y * _d.x);
	}

	inline bool operator == (const vec3d &_v) const
	{
		if (x == _v.x && y == _v.y && z == _v.z)
			return true;
		return false;
	}

	inline bool operator != (const vec3d &_v) const
	{
		if (x != _v.x || y != _v.y || z != _v.z)
			return true;
		return false;
	}

	inline double operator [] (int _i) const
	{
		const double *val = &x;
		return val[_i];
	}

	inline double length() const
	{
		double len = x * x + y * y + z * z;
		return (double) sqrt(len);
	}

	inline double squaredLength() const
	{
		return x * x + y * y + z * z;
	}

	inline double dot(const vec3d &_v) const
	{
		return x * _v.x + y * _v.y + z * _v.z;
	}

	inline void normalize()
	{
		double ln = length();
		if (!ln)
			return;
		double div = 1.0 / ln;
		x *= div;
		y *= div;
		z *= div;
	}

	inline void positive()
	{
		if (x < 0.0) x = -x;
		if (y < 0.0) y = -y;
		if (z < 0.0) z = -z;
	}

	inline void setLength(const double &len)
	{
		normalize();
		x *= len;
		y *= len;
		z *= len;
	}

	inline void zero()
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	double x=0, y=0, z=0;

};

inline vec3d operator * (double _m, const vec3d &_v)
{
	return vec3d(_v.x * _m, _v.y * _m, _v.z * _m);
}

} //NS MATH
} //NS FI

#endif