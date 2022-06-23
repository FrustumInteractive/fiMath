/*
 *  plane.cpp
 *  Engine
 *
 *  Created by Roger on 10-09-22.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "math/plane.h"

using namespace FI;
using namespace MATH;

Plane::Plane(const vec3 &origin, const vec3 &normal) 
{
	this->normal = normal;
	this->origin = origin;
	equation[0] = normal.x;
	equation[1] = normal.y;
	equation[2] = normal.z;
	equation[3] = -(normal.x*origin.x+normal.y*origin.y+normal.z*origin.z);
}

// Construct from triangle:
Plane::Plane(const vec3 &p1, const vec3 &p2, const vec3 &p3)
{
	normal = (p2-p1) | (p3-p1);
	normal.normalize();
	origin = p1;
	equation[0] = normal.x;
	equation[1] = normal.y;
	equation[2] = normal.z;
	equation[3] = -(normal.x*origin.x+normal.y*origin.y+normal.z*origin.z);
}

bool Plane::isFrontFacingTo(const vec3 &direction) const
{
	double dot = normal.dot(direction);
	return (dot <= 0);
}

float Plane::signedDistanceTo(const vec3 &point) const
{
	return (point.dot(normal)) + equation[3];
}
