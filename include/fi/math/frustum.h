#ifndef _FRUSTUM_H
#define _FRUSTUM_H

namespace FI {
namespace MATH {

class Matrix;

// This will allow us to create an object to keep track of our frustum
class Frustum 
{

public:

	// Call this every time the camera moves to update the frustum
	//void CalculateFrustum( const Matrix &modelview, const Matrix &projection);
	void calculateFrustum( const Matrix *modelViewProjection);

	// This takes a 3D point and returns TRUE if it's inside of the frustum
	bool pointInFrustum(float x, float y, float z);

	// This takes a 3D point and a radius and returns TRUE if the sphere is inside of the frustum
	bool sphereInFrustum(float x, float y, float z, float radius);

	// This takes the center and half the length of the cube.
	bool cubeInFrustum( float x, float y, float z, float size );

private:

	// This holds the A B C and D values for each side of our frustum.
	float m_Frustum[6][4];
};

} //ns MATH
} //ns FI

using FI::MATH::Frustum;

#endif
