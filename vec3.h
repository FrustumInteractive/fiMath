/*   3D Vector Class		*
 *   c Roger Dass 2011	    *
 *						    *
 *						    */

#ifndef VEC3_H
#define VEC3_H

#include <math.h>

namespace FI {
namespace MATH {

#define vec3f vec3

class vec3 
{
public:
	
	vec3(){};
	~vec3(){};

	
	inline vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z)
	{
	}
	
	inline vec3 operator - () const
	{
		return vec3(-x, -y, -z);
	}
	
	inline void operator -= (const vec3 &_v)
	{
		x -= _v.x;
		y -= _v.y;
		z -= _v.z;
	}
	
	inline void operator += (const vec3 &_v)
	{
		x += _v.x;
		y += _v.y;
		z += _v.z;
	}
	
	inline void operator *= (float _mul)
	{
		x *= _mul;
		y *= _mul;
		z *= _mul;
	}
	
	inline void operator *= (const vec3 &_v)
	{
		x *= _v.x;
		y *= _v.y;
		z *= _v.z;
	}
	
	inline void operator /= (float _div)
	{
		float mul = 1.0f / _div;
		x *= mul;
		y *= mul;
		z *= mul;
	}
	
	inline vec3 operator - (const vec3 &_v) const
	{
		return vec3(x - _v.x, y - _v.y, z - _v.z);
	}
	
	inline vec3 operator + (const vec3 &_v) const
	{
		return vec3(x + _v.x, y + _v.y, z + _v.z);
	}
	
	inline vec3 operator * (const vec3 &_v) const
	{
		return vec3(x * _v.x, y * _v.y, z * _v.z);
	}
	
	inline vec3 operator * (float _m) const
	{
		return vec3(x * _m, y * _m, z * _m);
	}
	
	inline vec3 operator / (const vec3 &_v) const
	{
		return vec3(x / _v.x, y / _v.y, z / _v.z);
	}
	
	inline vec3 operator / (float _d) const
	{
		float m = 1.0f / _d;
		return vec3(x * m, y * m, z * m);
	}
	
	inline vec3 operator | (const vec3 &_d) const
	{
		return vec3(y * _d.z - z * _d.y,
					 z * _d.x - x * _d.z,
					 x * _d.y - y * _d.x);
	}
	
	inline bool operator == (const vec3 &_v) const
	{
		if (x == _v.x && y == _v.y && z == _v.z)
			return true;
		return false;
	}
	
	inline bool operator != (const vec3 &_v) const
	{
		if (x != _v.x || y != _v.y || z != _v.z)
			return true;
		return false;
	}
	
	inline float operator [] (int _i) const
	{
		const float *val = &x;
		return val[_i];
	}
	
	inline float length() const
	{
		float len = x * x + y * y + z * z;
		return (float) sqrtf(len);
	}
	
	inline float squaredLength() const
	{
		return x * x + y * y + z * z;
	}
	
	inline float dot(const vec3 &_v) const
	{
		return x * _v.x + y * _v.y + z * _v.z;
	}
	
	inline void normalize()
	{
		float ln = length();
		if (!ln)
			return;
		float div = 1.0f / ln;
		x *= div;
		y *= div;
		z *= div;
	}
	
	inline void positive()
	{
		if (x < 0.f) x = -x;
		if (y < 0.f) y = -y;
		if (z < 0.f) z = -z;
	}
	
	inline void setLength(const float &len)
	{
		normalize();
		x *= len;
		y *= len;
		z *= len;
	}
    
    inline void zero()
    {
        x = 0.0f;
        y = 0.0f;
        z = 0.0f;
    }
	
	float x, y, z;
    
  
};

inline vec3 operator * (float _m, const vec3 &_v)
{
	return vec3(_v.x * _m, _v.y * _m, _v.z * _m);
}


// double version

class vec3d 
{
public:

	inline vec3d()
	{
        x = 0.0;
        y = 0.0;
        z = 0.0;
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
	
	double x, y, z;
    
  
};

inline vec3d operator * (double _m, const vec3d &_v)
{
	return vec3d(_v.x * _m, _v.y * _m, _v.z * _m);
}


// --- Utility Funtions  
inline float dot(const vec3& a, const vec3& b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

inline vec3 cross( const vec3& a, const vec3& b )
{
	vec3 out( 
	a.y*b.z - b.y*a.z, 
	b.x*a.z - a.x*b.z, 
	a.x*b.y - a.y*b.x ); 
	return out;
}

inline vec3 normal( const vec3& v )
{
	vec3 out;
	float ln = v.length();
	if (!ln)
		return out;
	float div = 1.0f / ln;
	out = v*div;
	return out;
}
// ----

} //NS MATH
} //NS FI

#endif



