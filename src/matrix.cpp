#include "math/matrix.h"
#include <memory.h>
#include <math.h>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace FI;
using namespace MATH;

// COLUMN MAJOR

//Helpers
const float EPSILON = 0.00001f;

// compute cofactor of 3x3 minor matrix without sign
// input params are 9 elements of the minor matrix
// NOTE: The caller must know its sign.
float getCofactor(float m0, float m1, float m2,
                  float m3, float m4, float m5,
                  float m6, float m7, float m8)
{
    return m0 * (m4 * m8 - m5 * m7) -
           m1 * (m3 * m8 - m5 * m6) +
           m2 * (m3 * m7 - m4 * m6);
}


Matrix::Matrix() :
	m {
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f 
	}
{
}

float Matrix::determinant()
{
	return m[0] * getCofactor(m[5],m[6],m[7], m[9],m[10],m[11], m[13],m[14],m[15]) -
		m[1] * getCofactor(m[4],m[6],m[7], m[8],m[10],m[11], m[12],m[14],m[15]) +
		m[2] * getCofactor(m[4],m[5],m[7], m[8],m[9], m[11], m[12],m[13],m[15]) - 
		m[3] * getCofactor(m[4],m[5],m[6], m[8],m[9], m[10], m[12],m[13],m[14]);
}

void Matrix::invert()
{
	// get cofactors of minor matrices
	float cofactor0 = getCofactor(m[5],m[6],m[7], m[9],m[10],m[11], m[13],m[14],m[15]);
	float cofactor1 = getCofactor(m[4],m[6],m[7], m[8],m[10],m[11], m[12],m[14],m[15]);
	float cofactor2 = getCofactor(m[4],m[5],m[7], m[8],m[9], m[11], m[12],m[13],m[15]);
	float cofactor3 = getCofactor(m[4],m[5],m[6], m[8],m[9], m[10], m[12],m[13],m[14]);

	// get determinant
	float determinant = m[0] * cofactor0 - m[1] * cofactor1 + m[2] * cofactor2 - m[3] * cofactor3;
	if(fabs(determinant) <= EPSILON) {
		identity();
		return;
	}

	// get rest of cofactors for adj(M)
	float cofactor4 = getCofactor(m[1],m[2],m[3], m[9],m[10],m[11], m[13],m[14],m[15]);
	float cofactor5 = getCofactor(m[0],m[2],m[3], m[8],m[10],m[11], m[12],m[14],m[15]);
	float cofactor6 = getCofactor(m[0],m[1],m[3], m[8],m[9], m[11], m[12],m[13],m[15]);
	float cofactor7 = getCofactor(m[0],m[1],m[2], m[8],m[9], m[10], m[12],m[13],m[14]);

	float cofactor8 = getCofactor(m[1],m[2],m[3], m[5],m[6], m[7],  m[13],m[14],m[15]);
	float cofactor9 = getCofactor(m[0],m[2],m[3], m[4],m[6], m[7],  m[12],m[14],m[15]);
	float cofactor10= getCofactor(m[0],m[1],m[3], m[4],m[5], m[7],  m[12],m[13],m[15]);
	float cofactor11= getCofactor(m[0],m[1],m[2], m[4],m[5], m[6],  m[12],m[13],m[14]);

	float cofactor12= getCofactor(m[1],m[2],m[3], m[5],m[6], m[7],  m[9], m[10],m[11]);
	float cofactor13= getCofactor(m[0],m[2],m[3], m[4],m[6], m[7],  m[8], m[10],m[11]);
	float cofactor14= getCofactor(m[0],m[1],m[3], m[4],m[5], m[7],  m[8], m[9], m[11]);
	float cofactor15= getCofactor(m[0],m[1],m[2], m[4],m[5], m[6],  m[8], m[9], m[10]);

	// build inverse matrix = adj(M) / det(M)
	// adjugate of M is the transpose of the cofactor matrix of M
	float invDeterminant = 1.0f / determinant;
	m[0] =  invDeterminant * cofactor0;
	m[1] = -invDeterminant * cofactor4;
	m[2] =  invDeterminant * cofactor8;
	m[3] = -invDeterminant * cofactor12;

	m[4] = -invDeterminant * cofactor1;
	m[5] =  invDeterminant * cofactor5;
	m[6] = -invDeterminant * cofactor9;
	m[7] =  invDeterminant * cofactor13;

	m[8] =  invDeterminant * cofactor2;
	m[9] = -invDeterminant * cofactor6;
	m[10]=  invDeterminant * cofactor10;
	m[11]= -invDeterminant * cofactor14;

	m[12]= -invDeterminant * cofactor3;
	m[13]=  invDeterminant * cofactor7;
	m[14]= -invDeterminant * cofactor11;
	m[15]=  invDeterminant * cofactor15;
}

void Matrix::perspectiveLH( float width, float height, float znear, float zfar )
{
	float znmzf = znear - zfar;
	float zfmzn = zfar - znear;
	float dzn = 2.0f * znear;

	m[0]=dzn/width;		m[4]=0;			m[8]=0;					m[12]=0;
	m[1]=0;				m[5]=dzn/height;	m[9]=0;					m[13]=0;
	m[2]=0;				m[6]=0;			m[10]=zfar/zfmzn;		m[14]=1;
	m[3]=0;				m[7]=0;			m[11]=znear*zfar/znmzf;	m[15]=0;
}

void Matrix::perspectiveRH( float width, float height, float znear, float zfar )
{
	float znmzf = znear - zfar;
	float dzn = 2.0f * znear;

	m[0]=dzn/width;		m[4]=0;			m[8]=0;					m[12]=0;
	m[1]=0;				m[5]=dzn/height;	m[9]=0;				m[13]=0;
	m[2]=0;				m[6]=0;			m[10]=zfar/znmzf;		m[14]=-1;
	m[3]=0;				m[7]=0;			m[11]=znear*zfar/znmzf;	m[15]=0;
}

void Matrix::Matrix::perspectiveFovLH( float fovY, float aspect, float znear, float zfar ) //verified CM
{
	float h = 1.0f/tanf(fovY/2.0f);
	float w = h / aspect;
	float zfmzn = zfar - znear;

	// the following looks like RM, but it's sequenced as CM
	float i[16] = {
		w, 					0, 					0, 					0,
		0, 					h, 					0, 					0,
		0, 					0, 					zfar/zfmzn,			1,
		0,					0,					-znear*zfar/zfmzn,	0
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::perspectiveFovRH( float fovY, float aspect, float zn, float zf )
{
	float h = 1.0f/tanf(fovY/2.0f);
	float w = h / aspect;
	float znmzf = zn-zf;
	// the following looks like RM, but it's sequenced as CM
	float i[16] = {
		w, 					0, 					0, 					0,
		0, 					h, 					0, 					0,
		0, 					0, 					zf/znmzf,			-1,
		0,					0,					zn*zf/znmzf,		0
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::orthoLH( float w, float h, float zn, float zf) //row major?
{
	float i[16] = {
		1/w,				0, 					0, 					0,
		0, 					1/h,				0, 					0,
		0, 					0, 					-2/(zf-zn),			zn/(zn-zf),
		0,					0,					-(zf+zn)/(zf-zn),	1
	};

	memcpy(m, i, sizeof(float)*16);;
}

void Matrix::orthoOffCenterLH( float l, float r, float b, float t, float zn, float zf)
{
	float i[16] = {
		2.0f/(r-l),		0, 				0, 				0,
		0, 				2.0f/(t-b),		0, 				0,
		0, 				0, 				1.0f/(zf-zn),	0,
		(l+r)/(l-r),	(t+b)/(b-t),	zn/(zn-zf),		1
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::orthoRH( float w, float h, float zn, float zf) //row major?
{
	float i[16] = {
		2/w,		0, 			0, 				0,
		0, 			2/h,		0, 				0,
		0, 			0, 			1.0f/(zn-zf),	zn/(zn-zf),
		0,			0,			0,				1
	};
	
	memcpy(m, i, sizeof(float)*16);
}

void Matrix::orthoOffCenterRH( float l, float r, float b, float t, float zn, float zf)
{
	float i[16] = {
		2.0f/(r-l),			0, 					0, 					0,
		0, 					2.0f/(t-b),			0, 					0,
		0, 					0, 					1.0f/(zn-zf),		0,
		(l+r)/(l-r),		(t+b)/(b-t),		zn/(zn-zf),			1
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::lookAtLH( const vec3& eye, const vec3& target, const vec3& up) //verified CM
{
	vec3 zaxis = target - eye;
	zaxis.normalize();

	vec3 xaxis = up | zaxis; //cross prod
	xaxis.normalize();

	vec3 yaxis = zaxis | xaxis;

	float i[16] = {
		xaxis.x,			yaxis.x,			zaxis.x,			0,
		xaxis.y,			yaxis.y,			zaxis.y,			0,
		xaxis.z,			yaxis.z,			zaxis.z,			0,
		-xaxis.dot(eye),	-yaxis.dot(eye),	-zaxis.dot(eye),	1
	};
	
	memcpy(m, i, sizeof(float)*16);
}

void Matrix::lookAtRH( const vec3& eye, const vec3& target, const vec3& up) //verified CM
{
	vec3 zaxis = normal(eye - target);
	vec3 xaxis = normal(cross(up, zaxis));
	vec3 yaxis = cross(zaxis, xaxis);

	float i[16] = {
		-xaxis.x,			yaxis.x,			-zaxis.x,			0,
		-xaxis.y,			yaxis.y,			-zaxis.y,			0,
		-xaxis.z,			yaxis.z,			-zaxis.z,			0,
		xaxis.dot(eye),		-yaxis.dot(eye),	-zaxis.dot(eye),	1
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::rotate( float angle, float x, float y, float z )
{
	rotate(angle, vec3(x,y,z));
}

void Matrix::rotate( float angle, const vec3& axis ) //verified CM
{
	const float c = cosf(angle);
	const float s = sinf(angle);

	float x = axis.x;
	float y = axis.y;
	float z = axis.z;

	float omc = 1.0f - c;
	float xy = x*y;
	float xz = x*z;
	float yz = y*z;
	float xs = x*s;
	float ys = y*s;
	float zs = z*s;

	float i[16] = {
		x*x*omc+c,		xy*omc-zs,		xz*omc+ys,		0,
		xy*omc+zs,		y*y*omc+c,		yz*omc-xs,		0,
		xz*omc-ys,		yz*omc+xs,		z*z*omc+c,		0,
		0,				0,				0,				1
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::translate( float x, float y, float z )
{
	translate(vec3(x,y,z));
}

void Matrix::translate( const vec3& p )	// verified CM
{
	float i[16] = {
		1,		0, 		0, 		0,
		0,		1,		0,		0,
		0,		0,		1,		0,
		p.x,	p.y,	p.z,	1
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::scale( float x, float y, float z )
{
	scale(vec3(x,y,z));
}

void Matrix::scale( const vec3& s )
{
	float i[16] =
	{	
		s.x,	0.0f,	0.0f,	0.0f,
		0.0f,	s.y,	0.0f,	0.0f,
		0.0f,	0.0f,	s.z,	0.0f,
		0.0f,	0.0f,	0.0f,	1.0f	
	};

	memcpy(m, i, sizeof(float)*16);
}

void Matrix::identity()
{
	m[0] = 1.0f;	m[4] = 0.0f;	m[8] = 0.0f;	m[12] = 0.0f;
	m[1] = 0.0f;	m[5] = 1.0f;	m[9] = 0.0f;	m[13] = 0.0f;
	m[2] = 0.0f;	m[6] = 0.0f;	m[10] = 1.0f;	m[14] = 0.0f;
	m[3] = 0.0f;	m[7] = 0.0f;	m[11] = 0.0f;	m[15] = 1.0f;	
}

void Matrix::transpose()
{
	float out[16];
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			out[j*4+i] = m[i*4+j];
	memcpy(m, out, 64);
	// todo: speed up alot by std::move (avoid copy)
}

void Matrix::RHToLH() //not sure if this is CM
{
	float rx=m[0], ry=m[4], rz=m[8],  rw=m[12];
	float lx=m[1], ly=m[5], lz=m[9],  lw=m[13];
	float ux=m[2], uy=m[6], uz=m[10], uw=m[14];
	float px=m[3], py=m[7], pz=m[11], pw=m[15];

	m[0]=rx;	m[4]=rz;	m[8]=ry;	m[12]=rw;
	m[1]=ux;	m[5]=uz;	m[9]=uy;	m[13]=uw;
	m[2]=lx;	m[6]=lz;	m[10]=ly;	m[14]=lw;
	m[3]=px;	m[7]=pz;	m[11]=py;	m[15]=pw;
}

void Matrix::LHToRH()	// not sure if this is CM
{
	float rx=m[0], ry=m[4], rz=m[8],  rw=m[12];
	float ux=m[1], uy=m[5], uz=m[9],  uw=m[13];
	float lx=m[2], ly=m[6], lz=m[10], lw=m[14];
	float px=m[3], py=m[7], pz=m[11], pw=m[15];

	m[0]=rx;	m[4]=rz;	m[8]=ry;	m[12]=rw;
	m[1]=lx;	m[5]=lz;	m[9]=ly;	m[13]=lw;
	m[2]=ux;	m[6]=uz;	m[10]=uy;	m[14]=uw;
	m[3]=px;	m[7]=pz;	m[11]=py;	m[15]=pw;
}


void Matrix::printRM(const char *prefix) const
{
    const char *pre;
    if(prefix)
    {
        pre = prefix;
    }
    
    int s = 12;
    cout << right << setprecision(2) << setiosflags(ios::fixed) << setfill(' ');
    cout << pre << setw(s)<<m[0] << setw(s)<<m[1] << setw(s)<<m[2] << setw(s)<<m[3] << endl;
    cout << pre << setw(s)<<m[4] << setw(s)<<m[5] << setw(s)<<m[6] << setw(s)<<m[7] << endl;
    cout << pre << setw(s)<<m[8] << setw(s)<<m[9] << setw(s)<<m[10]<< setw(s)<<m[11]<< endl;
    cout << pre << setw(s)<<m[12]<< setw(s)<<m[13]<< setw(s)<<m[14]<< setw(s)<<m[15]<< endl;
}

void Matrix::printCM(const char *prefix) const
{
    const char *pre = "";
    if(prefix)
    {
        pre = prefix;
    }
    
    int s=12;
    cout << right << setprecision(2) << setiosflags(ios::fixed) << setfill(' ');
    cout << pre << setw(s)<<m[0] << setw(s)<<m[4] << setw(s)<<m[8] << setw(s)<<m[12] << endl;
    cout << pre << setw(s)<<m[1] << setw(s)<<m[5] << setw(s)<<m[9] << setw(s)<<m[13] << endl;
    cout << pre << setw(s)<<m[2] << setw(s)<<m[6] << setw(s)<<m[10] << setw(s)<<m[14] << endl;
    cout << pre << setw(s)<<m[3] << setw(s)<<m[7] << setw(s)<<m[11] << setw(s)<<m[15] << endl;
}

Matrix Matrix::inverse() {
	Matrix out;
	out=m; out.invert();
	return out;
}

// Operator overloads
Matrix Matrix::operator*( const Matrix& mat ) const
{
	Matrix out;

	Matrix16f::mult(out.m, m, mat.m);
	return out;
}

vec3 Matrix::operator*( const vec3& v ) const
{
	float o[4];

	Matrix16f::mult(o, m, v);
	return vec3(o[0], o[1], o[2]);
}

// assignment operator
Matrix& Matrix::operator=(const Matrix& mat)
{
	const float *i = mat.data();
	m[0]=i[0];   m[1]=i[1];   m[2]=i[2];   m[3]=i[3];
	m[4]=i[4];   m[5]=i[5];   m[6]=i[6];   m[7]=i[7];
	m[8]=i[8];   m[9]=i[9];   m[10]=i[10]; m[11]=i[11];
	m[12]=i[12]; m[13]=i[13]; m[14]=i[14]; m[15]=i[15];
	return *this;
}

Matrix& Matrix::operator=(const float a[16])
{
	m[0]=a[0];   m[1]=a[1];   m[2]=a[2];   m[3]=a[3];
	m[4]=a[4];   m[5]=a[5];   m[6]=a[6];   m[7]=a[7];
	m[8]=a[8];   m[9]=a[9];   m[10]=a[10]; m[11]=a[11];
	m[12]=a[12]; m[13]=a[13]; m[14]=a[14]; m[15]=a[15];
	return *this;
}

// equality comparison. doesn't modify object. therefore const.
bool Matrix::operator==(const Matrix& mat) const
{
	const float *i = mat.data();
	return
			(
					m[0]==i[0] && m[1]==i[1] && m[2]==i[2] && m[3]==i[3] &&
					m[4]==i[4] && m[5]==i[5] && m[6]==i[6] && m[7]==i[7] &&
					m[8]==i[8] && m[9]==i[9] && m[10]==i[10] && m[11]==i[11] &&
					m[12]==i[12] && m[13]==i[13] && m[14]==i[14] && m[15]==i[15]
			);
}


// -- MATRIX16F --

void Matrix16f::mult( float *dst, const float *a, const float *b )
{
	dst[0]  = a[0]*b[0]  +  a[1]*b[4]  +  a[2]*b[8]   + a[3]*b[12];
	dst[1]  = a[0]*b[1]  +  a[1]*b[5]  +  a[2]*b[9]   + a[3]*b[13];
	dst[2]  = a[0]*b[2]  +  a[1]*b[6]  +  a[2]*b[10]  + a[3]*b[14];
	dst[3]  = a[0]*b[3]  +  a[1]*b[7]  +  a[2]*b[11]  + a[3]*b[15];

	dst[4]  = a[4]*b[0]  +  a[5]*b[4]  +  a[6]*b[8]   + a[7]*b[12];
	dst[5]  = a[4]*b[1]  +  a[5]*b[5]  +  a[6]*b[9]   + a[7]*b[13];
	dst[6]  = a[4]*b[2]  +  a[5]*b[6]  +  a[6]*b[10]  + a[7]*b[14];
	dst[7]  = a[4]*b[3]  +  a[5]*b[7]  +  a[6]*b[11]  + a[7]*b[15];

	dst[8]  = a[8]*b[0]  +  a[9]*b[4]  +  a[10]*b[8]  + a[11]*b[12];
	dst[9]  = a[8]*b[1]  +  a[9]*b[5]  +  a[10]*b[9]  + a[11]*b[13];
	dst[10] = a[8]*b[2]  +  a[9]*b[6]  +  a[10]*b[10] + a[11]*b[14];
	dst[11] = a[8]*b[3]  +  a[9]*b[7]  +  a[10]*b[11] + a[11]*b[15];

	dst[12] = a[12]*b[0] +  a[13]*b[4] +  a[14]*b[8]  + a[15]*b[12];
	dst[13] = a[12]*b[1] +  a[13]*b[5] +  a[14]*b[9]  + a[15]*b[13];
	dst[14] = a[12]*b[2] +  a[13]*b[6] +  a[14]*b[10] + a[15]*b[14];
	dst[15] = a[12]*b[3] +  a[13]*b[7] +  a[14]*b[11] + a[15]*b[15];
}

void Matrix16f::mult(float *dst, const vec3 &a, const float *b)
{
	/* row major v * M
	x*a + y*e + z*i + w*m
	x*b + y*f + z*j + w*n
	x*c + y*g + z*k + w*o
	x*d + y*h + z*l + w*p
	*/

	float x = a.x;
	float y = a.y;
	float z = a.z;

	dst[0] = x*b[0] + y*b[1] + z*b[2] + 1.0f*b[3];
	dst[1] = x*b[4] + y*b[5] + z*b[6] + 1.0f*b[7];
	dst[2] = x*b[8] + y*b[9] + z*b[10] + 1.0f*b[11];
	dst[3] = x*b[12] + y*b[13] + z*b[14] + 1.0f*b[15];
}

void Matrix16f::mult(float *dst, const float *a, const vec3 &b)
{
	/* row major M * v
	x*a + y*b + z*c + w*d;
	x*e + y*f + z*g + w*h;
	x*i + y*j + z*k + w*l;
	x*m + y*n + z*o + w*p;
	*/

	float x = b.x;
	float y = b.y;
	float z = b.z;

	dst[0] = x*a[0] + y*a[4] + z*a[8] + 1.0f*a[12];
	dst[1] = x*a[1] + y*a[5] + z*a[9] + 1.0f*a[13];
	dst[2] = x*a[2] + y*a[6] + z*a[10] + 1.0f*a[14];
	dst[3] = x*a[3] + y*a[7] + z*a[11] + 1.0f*a[15];
}
