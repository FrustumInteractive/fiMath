/*
 *  Matrix.h
 *  Matrix
 *
 *  Created by Roger Dass on 10-10-13.
 *  Copyright 2010 Roger Dass. All rights reserved.
 *
 */

#ifndef _FI_MATRIX_H_
#define _FI_MATRIX_H_

#include "vec3.h"
#include "memory.h"

namespace FI { namespace MATH {

// RAW FLOAT MATRIX FUNTIONS
namespace Matrix16f
{
	void mult(float *dst, const float *a, const float *b);
	void mult(float dst[4], const vec3 &a, const float *b);
	void mult(float dst[4], const float *a, const vec3 &b);
}

// Column Major
class Matrix
{
public:

	Matrix();

	Matrix (const Matrix &mat) {
		memcpy(m, mat.data(), sizeof(float)*16);
	}

	Matrix (const float i[16]) {
		memcpy(m, i, sizeof(float)*16);
	}

	void perspectiveLH( float width, float height, float znear, float zfar );
	void perspectiveRH( float width, float height, float znear, float zfar );
	void perspectiveFovLH( float fovY, float aspect, float znear, float zfar );
	void perspectiveFovRH( float fovY, float aspect, float znear, float zfar );
	void orthoOffCenterRH( float left, float right, float bottom, float top, float znear, float zfar);
	void orthoRH( float width, float height, float znear, float zfar); // (0,0) at center
	void orthoOffCenterLH( float left, float right, float bottom, float top, float znear, float zfar);
	void orthoLH( float width, float height, float znear, float zfar); // (0,0) at center
	
	void lookAtLH( const vec3& eye, const vec3& target, const vec3& up);
	void lookAtRH( const vec3& eye, const vec3& target, const vec3& up);

	void rotate( float angle, const vec3& axis ); // x,y,z must be unit vector
	void rotate( float angle, float axis_x, float axis_y, float axis_z );
	void translate( const vec3& p );
	void translate( float x, float y, float z );
	void scale( const vec3& s );
	void scale( float x, float y, float z );
	
	void identity();
	void transpose();
	void invert();

	void RHToLH();
	void LHToRH();

	void printCM(const char *prefix=0) const;
	void printRM(const char *prefix=0) const;

	const float *data() const {return m;}
	float determinant();

	Matrix inverse();

	// Operator overloads
	Matrix operator*( const Matrix& mat ) const;

	vec3 operator*( const vec3& v ) const;

	// assignment operator
	Matrix& operator=(const Matrix& mat);
	Matrix& operator=(const float a[16]);

	// equality comparison. doesn't modify object. therefore const.
	bool operator==(const Matrix& mat) const;

	float m[16];
};

} //ns MATH
} //NS FI

#endif /*_FI_MATRIX_H_*/
