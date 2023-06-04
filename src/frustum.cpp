#include "fi/math/frustum.h"
#include "fi/math/matrix.h"

#include <math.h>
#include <memory.h>								

using namespace FI;
using namespace MATH;

// We create an enum of the sides so we don't have to call each side 0 or 1.
// This way it makes it more understandable and readable when dealing with frustum sides.
enum FrustumSide
{
	RIGHT	= 0,		// The RIGHT side of the frustum
	LEFT	= 1,		// The LEFT	 side of the frustum
	BOTTOM	= 2,		// The BOTTOM side of the frustum
	TOP		= 3,		// The TOP side of the frustum
	BACK	= 4,		// The BACK	side of the frustum
	FRONT	= 5			// The FRONT side of the frustum
}; 

// Like above, instead of saying a number for the ABC and D of the plane, we
// want to be more descriptive.
enum PlaneData
{
	A = 0,				// The X value of the plane's normal
	B = 1,				// The Y value of the plane's normal
	C = 2,				// The Z value of the plane's normal
	D = 3				// The distance the plane is from the origin
};

///////////////////////////////// NORMALIZE PLANE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This normalizes a plane (A side) from a given frustum.
/////
///////////////////////////////// NORMALIZE PLANE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

void NormalizePlane(float frustum[6][4], int side)
{
	// Here we calculate the magnitude of the normal to the plane (point A B C)
	// Remember that (A, B, C) is that same thing as the normal's (X, Y, Z).
	// To calculate magnitude you use the equation:  magnitude = sqrt( x^2 + y^2 + z^2)
	float magnitude = (float)sqrtf( frustum[side][A] * frustum[side][A] + 
								   frustum[side][B] * frustum[side][B] + 
								   frustum[side][C] * frustum[side][C] );

	// Then we divide the plane's values by it's magnitude.
	// This makes it easier to work with.
	frustum[side][A] /= magnitude;
	frustum[side][B] /= magnitude;
	frustum[side][C] /= magnitude;
	frustum[side][D] /= magnitude; 
}


///////////////////////////////// CALCULATE FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This extracts our frustum from the projection and modelview matrix.
/////
///////////////////////////////// CALCULATE FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/*
void Frustum::CalculateFrustum( const Matrix &modl, const Matrix &proj )
{    
	//float   proj[16];								// This will hold our projection matrix
	//float   modl[16];								// This will hold our modelview matrix
	float   clip[16];								// This will hold the clipping planes


	// Now that we have our modelview and projection matrix, if we combine these 2 matrices,
	// it will give us our clipping planes.  To combine 2 matrices, we multiply them.

	clip[ 0] = modl.m[ 0] * proj.m[ 0] + modl.m[ 1] * proj.m[ 4] + modl.m[ 2] * proj.m[ 8] + modl.m[ 3] * proj.m[12];
	clip[ 1] = modl.m[ 0] * proj.m[ 1] + modl.m[ 1] * proj.m[ 5] + modl.m[ 2] * proj.m[ 9] + modl.m[ 3] * proj.m[13];
	clip[ 2] = modl.m[ 0] * proj.m[ 2] + modl.m[ 1] * proj.m[ 6] + modl.m[ 2] * proj.m[10] + modl.m[ 3] * proj.m[14];
	clip[ 3] = modl.m[ 0] * proj.m[ 3] + modl.m[ 1] * proj.m[ 7] + modl.m[ 2] * proj.m[11] + modl.m[ 3] * proj.m[15];

	clip[ 4] = modl.m[ 4] * proj.m[ 0] + modl.m[ 5] * proj.m[ 4] + modl.m[ 6] * proj.m[ 8] + modl.m[ 7] * proj.m[12];
	clip[ 5] = modl.m[ 4] * proj.m[ 1] + modl.m[ 5] * proj.m[ 5] + modl.m[ 6] * proj.m[ 9] + modl.m[ 7] * proj.m[13];
	clip[ 6] = modl.m[ 4] * proj.m[ 2] + modl.m[ 5] * proj.m[ 6] + modl.m[ 6] * proj.m[10] + modl.m[ 7] * proj.m[14];
	clip[ 7] = modl.m[ 4] * proj.m[ 3] + modl.m[ 5] * proj.m[ 7] + modl.m[ 6] * proj.m[11] + modl.m[ 7] * proj.m[15];

	clip[ 8] = modl.m[ 8] * proj.m[ 0] + modl.m[ 9] * proj.m[ 4] + modl.m[10] * proj.m[ 8] + modl.m[11] * proj.m[12];
	clip[ 9] = modl.m[ 8] * proj.m[ 1] + modl.m[ 9] * proj.m[ 5] + modl.m[10] * proj.m[ 9] + modl.m[11] * proj.m[13];
	clip[10] = modl.m[ 8] * proj.m[ 2] + modl.m[ 9] * proj.m[ 6] + modl.m[10] * proj.m[10] + modl.m[11] * proj.m[14];
	clip[11] = modl.m[ 8] * proj.m[ 3] + modl.m[ 9] * proj.m[ 7] + modl.m[10] * proj.m[11] + modl.m[11] * proj.m[15];

	clip[12] = modl.m[12] * proj.m[ 0] + modl.m[13] * proj.m[ 4] + modl.m[14] * proj.m[ 8] + modl.m[15] * proj.m[12];
	clip[13] = modl.m[12] * proj.m[ 1] + modl.m[13] * proj.m[ 5] + modl.m[14] * proj.m[ 9] + modl.m[15] * proj.m[13];
	clip[14] = modl.m[12] * proj.m[ 2] + modl.m[13] * proj.m[ 6] + modl.m[14] * proj.m[10] + modl.m[15] * proj.m[14];
	clip[15] = modl.m[12] * proj.m[ 3] + modl.m[13] * proj.m[ 7] + modl.m[14] * proj.m[11] + modl.m[15] * proj.m[15];
	
	// Now we actually want to get the sides of the frustum.  To do this we take
	// the clipping planes we received above and extract the sides from them.

	// This will extract the RIGHT side of the frustum
	m_Frustum[RIGHT][A] = clip[ 3] - clip[ 0];
	m_Frustum[RIGHT][B] = clip[ 7] - clip[ 4];
	m_Frustum[RIGHT][C] = clip[11] - clip[ 8];
	m_Frustum[RIGHT][D] = clip[15] - clip[12];

	// Now that we have a normal (A,B,C) and a distance (D) to the plane,
	// we want to normalize that normal and distance.

	// Normalize the RIGHT side
	NormalizePlane(m_Frustum, RIGHT);

	// This will extract the LEFT side of the frustum
	m_Frustum[LEFT][A] = clip[ 3] + clip[ 0];
	m_Frustum[LEFT][B] = clip[ 7] + clip[ 4];
	m_Frustum[LEFT][C] = clip[11] + clip[ 8];
	m_Frustum[LEFT][D] = clip[15] + clip[12];

	// Normalize the LEFT side
	NormalizePlane(m_Frustum, LEFT);

	// This will extract the BOTTOM side of the frustum
	m_Frustum[BOTTOM][A] = clip[ 3] + clip[ 1];
	m_Frustum[BOTTOM][B] = clip[ 7] + clip[ 5];
	m_Frustum[BOTTOM][C] = clip[11] + clip[ 9];
	m_Frustum[BOTTOM][D] = clip[15] + clip[13];

	// Normalize the BOTTOM side
	NormalizePlane(m_Frustum, BOTTOM);

	// This will extract the TOP side of the frustum
	m_Frustum[TOP][A] = clip[ 3] - clip[ 1];
	m_Frustum[TOP][B] = clip[ 7] - clip[ 5];
	m_Frustum[TOP][C] = clip[11] - clip[ 9];
	m_Frustum[TOP][D] = clip[15] - clip[13];

	// Normalize the TOP side
	NormalizePlane(m_Frustum, TOP);

	// This will extract the BACK side of the frustum
	m_Frustum[BACK][A] = clip[ 3] - clip[ 2];
	m_Frustum[BACK][B] = clip[ 7] - clip[ 6];
	m_Frustum[BACK][C] = clip[11] - clip[10];
	m_Frustum[BACK][D] = clip[15] - clip[14];

	// Normalize the BACK side
	NormalizePlane(m_Frustum, BACK);

	// This will extract the FRONT side of the frustum
	m_Frustum[FRONT][A] = clip[ 3] + clip[ 2];
	m_Frustum[FRONT][B] = clip[ 7] + clip[ 6];
	m_Frustum[FRONT][C] = clip[11] + clip[10];
	m_Frustum[FRONT][D] = clip[15] + clip[14];

	// Normalize the FRONT side
	NormalizePlane(m_Frustum, FRONT);
}
*/
void Frustum::calculateFrustum( const Matrix *mat )
{    
	const float *m = mat->data();
	// Assume we have modelview and projection matrix combined.
	// this will give us our clipping planes. To combine 2 matrices, we multiply them.
	
	// Now we actually want to get the sides of the frustum.  To do this we take
	// the clipping planes we received above and extract the sides from them.

	// This will extract the RIGHT side of the frustum
	m_Frustum[RIGHT][A] = m[ 3] - m[ 0];
	m_Frustum[RIGHT][B] = m[ 7] - m[ 4];
	m_Frustum[RIGHT][C] = m[11] - m[ 8];
	m_Frustum[RIGHT][D] = m[15] - m[12];

	// Now that we have a normal (A,B,C) and a distance (D) to the plane,
	// we want to normalize that normal and distance.

	// Normalize the RIGHT side
	NormalizePlane(m_Frustum, RIGHT);

	// This will extract the LEFT side of the frustum
	m_Frustum[LEFT][A] = m[ 3] + m[ 0];
	m_Frustum[LEFT][B] = m[ 7] + m[ 4];
	m_Frustum[LEFT][C] = m[11] + m[ 8];
	m_Frustum[LEFT][D] = m[15] + m[12];

	// Normalize the LEFT side
	NormalizePlane(m_Frustum, LEFT);

	// This will extract the BOTTOM side of the frustum
	m_Frustum[BOTTOM][A] = m[ 3] + m[ 1];
	m_Frustum[BOTTOM][B] = m[ 7] + m[ 5];
	m_Frustum[BOTTOM][C] = m[11] + m[ 9];
	m_Frustum[BOTTOM][D] = m[15] + m[13];

	// Normalize the BOTTOM side
	NormalizePlane(m_Frustum, BOTTOM);

	// This will extract the TOP side of the frustum
	m_Frustum[TOP][A] = m[ 3] - m[ 1];
	m_Frustum[TOP][B] = m[ 7] - m[ 5];
	m_Frustum[TOP][C] = m[11] - m[ 9];
	m_Frustum[TOP][D] = m[15] - m[13];

	// Normalize the TOP side
	NormalizePlane(m_Frustum, TOP);

	// This will extract the BACK side of the frustum
	m_Frustum[BACK][A] = m[ 3] - m[ 2];
	m_Frustum[BACK][B] = m[ 7] - m[ 6];
	m_Frustum[BACK][C] = m[11] - m[10];
	m_Frustum[BACK][D] = m[15] - m[14];

	// Normalize the BACK side
	NormalizePlane(m_Frustum, BACK);

	// This will extract the FRONT side of the frustum
	m_Frustum[FRONT][A] = m[ 3] + m[ 2];
	m_Frustum[FRONT][B] = m[ 7] + m[ 6];
	m_Frustum[FRONT][C] = m[11] + m[10];
	m_Frustum[FRONT][D] = m[15] + m[14];

	// Normalize the FRONT side
	NormalizePlane(m_Frustum, FRONT);
}

// The code below will allow us to make checks within the frustum.  For example,
// if we want to see if a point, a sphere, or a cube lies inside of the frustum.
// Because all of our planes point INWARDS (The normals are all pointing inside the frustum)
// we then can assume that if a point is in FRONT of all of the planes, it's inside.

///////////////////////////////// POINT IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This determines if a point is inside of the frustum
/////
///////////////////////////////// POINT IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

bool Frustum::pointInFrustum( float x, float y, float z )
{
	// Go through all the sides of the frustum
	for(unsigned char i = 0; i < 6; i++ )
	{
		// Calculate the plane equation and check if the point is behind a side of the frustum
		if(m_Frustum[i][A] * x + m_Frustum[i][B] * y + m_Frustum[i][C] * z + m_Frustum[i][D] <= 0)
		{
			// The point was behind a side, so it ISN'T in the frustum
			return false;
		}
	}

	// The point was inside of the frustum (In front of ALL the sides of the frustum)
	return true;
}


///////////////////////////////// SPHERE IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This determines if a sphere is inside of our frustum by it's center and radius.
/////
///////////////////////////////// SPHERE IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

bool Frustum::sphereInFrustum( float x, float y, float z, float radius )
{
	// Go through all the sides of the frustum
	for(unsigned char i = 0; i < 6; i++ )	
	{
		// If the center of the sphere is farther away from the plane than the radius
		if( m_Frustum[i][A] * x + m_Frustum[i][B] * y + m_Frustum[i][C] * z + m_Frustum[i][D] <= -radius )
		{
			// The distance was greater than the radius so the sphere is outside of the frustum
			return false;
		}
	}
	
	// The sphere was inside of the frustum!
	return true;
}


///////////////////////////////// CUBE IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*
/////
/////	This determines if a cube is in or around our frustum by it's center and 1/2 it's length
/////
///////////////////////////////// CUBE IN FRUSTUM \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*

bool Frustum::cubeInFrustum( float x, float y, float z, float size )
{
	// Basically, what is going on is, that we are given the center of the cube,
	// and half the length.  Think of it like a radius.  Then we checking each point
	// in the cube and seeing if it is inside the frustum.  If a point is found in front
	// of a side, then we skip to the next side.  If we get to a plane that does NOT have
	// a point in front of it, then it will return false.

	// *Note* - This will sometimes say that a cube is inside the frustum when it isn't.
	// This happens when all the corners of the bounding box are not behind any one plane.
	// This is rare and shouldn't effect the overall rendering speed.

	for(unsigned char i = 0; i < 6; i++ )
	{
		if(m_Frustum[i][A] * (x - size) + m_Frustum[i][B] * (y - size) + m_Frustum[i][C] * (z - size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x + size) + m_Frustum[i][B] * (y - size) + m_Frustum[i][C] * (z - size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x - size) + m_Frustum[i][B] * (y + size) + m_Frustum[i][C] * (z - size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x + size) + m_Frustum[i][B] * (y + size) + m_Frustum[i][C] * (z - size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x - size) + m_Frustum[i][B] * (y - size) + m_Frustum[i][C] * (z + size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x + size) + m_Frustum[i][B] * (y - size) + m_Frustum[i][C] * (z + size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x - size) + m_Frustum[i][B] * (y + size) + m_Frustum[i][C] * (z + size) + m_Frustum[i][D] > 0)
		   continue;
		if(m_Frustum[i][A] * (x + size) + m_Frustum[i][B] * (y + size) + m_Frustum[i][C] * (z + size) + m_Frustum[i][D] > 0)
		   continue;

		// If we get here, it isn't in the frustum
		return false;
	}

	return true;
}



