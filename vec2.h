#ifndef VEC2_H
#define VEC2_H


#include <math.h>


#define Mag(Normal) (sqrt(Normal.x*Normal.x + Normal.y*Normal.y + Normal.z*Normal.z))

namespace FI { namespace MATH {

class vec2 
{
public:
	
	inline vec2()
	{
        x = 0.0f;
        y = 0.0f;
	}
	
	inline vec2(float _x, float _y) : x(_x), y(_y)
	{
	}
	
	inline vec2 operator - () const
	{
		return vec2(-x, -y);
	}
	
	inline void operator -= (const vec2 &_v)
	{
		x -= _v.x;
		y -= _v.y;
	}
	
	inline void operator += (const vec2 &_v)
	{
		x += _v.x;
		y += _v.y;
	}
	
	inline void operator *= (float _mul)
	{
		x *= _mul;
		y *= _mul;
	}
	
	inline void operator *= (const vec2 &_v)
	{
		x *= _v.x;
		y *= _v.y;
	}
	
	inline void operator /= (float _div)
	{
		float mul = 1.0f / _div;
		x *= mul;
		y *= mul;
	}
	
	inline vec2 operator - (const vec2 &_v) const
	{
		return vec2(x - _v.x, y - _v.y);
	}
	
	inline vec2 operator + (const vec2 &_v) const
	{
		return vec2(x + _v.x, y + _v.y);
	}
	
	inline vec2 operator * (const vec2 &_v) const
	{
		return vec2(x * _v.x, y * _v.y);
	}
	
	inline vec2 operator * (float _m) const
	{
		return vec2(x * _m, y * _m);
	}
	
	inline vec2 operator / (const vec2 &_v) const
	{
		return vec2(x / _v.x, y / _v.y);
	}
	
	inline vec2 operator / (float _d) const
	{
		float m = 1.0f / _d;
		return vec2(x * m, y * m);
	}
	
	/*inline vec2 operator | (const vec2 &_d) const
	{
		return vec3(y * _d.z - z * _d.y,
					z * _d.x - x * _d.z,
					x * _d.y - y * _d.x);
	}*/
	
	inline bool operator == (const vec2 &_v) const
	{
		if (x == _v.x && y == _v.y)
			return true;
		return false;
	}
	
	inline bool operator != (const vec2 &_v) const
	{
		if (x != _v.x || y != _v.y)
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
		float len = x * x + y * y;
		return sqrtf(len);
	}
	
	inline float squaredLength() const
	{
		return x * x + y * y;
	}
	
	inline float dot(const vec2 &_v) const
	{
		return x * _v.x + y * _v.y;
	}
	
	inline void normalize()
	{
		float ln = length();
		if (!ln)
			return;
		float div = 1.0f / ln;
		x *= div;
		y *= div;
	}
	
	inline void positive()
	{
		if (x < 0) x = -x;
		if (y < 0) y = -y;
	}
	
	inline void setLength(const float &len)
	{
		normalize();
		x *= len;
		y *= len;
	}
    
    inline void zero()
    {
        x = 0.0f;
        y = 0.0f;
    }
	
	float x, y;
};

inline vec2 operator * (float _m, const vec2 &_v)
{
	return vec2(_v.x * _m, _v.y * _m);
}

inline vec2 reflect( const vec2 &incident, const vec2 &normal )
{
	return incident - 2.0f * normal.dot(incident) * normal;
}

inline vec2 refract( const vec2 &incident, const vec2 &normal, const float &n ) //n is ratio
{
	float ndoti = normal.dot(incident);
	float k = 1.0f - n*n * (1.0f - ndoti * ndoti);
	if (k < 0.0f) return vec2(0.0f, 0.0f);
	else return n * incident - (n * ndoti + sqrtf(k)) * normal;
}

} //ns MATH
} //ns FI

#endif


