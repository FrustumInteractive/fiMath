#ifndef _PLANE_H
#define _PLANE_H

#include "math/vec3.h"

namespace FI {
namespace MATH {

class Plane
{
public:

        float equation[4];
        vec3 origin;
        vec3 normal;

        Plane(const vec3& origin, const vec3& normal);
        Plane(const vec3& p1, const vec3& p2, const vec3& p3);

        bool isFrontFacingTo(const vec3& direction) const;
        float signedDistanceTo(const vec3& point) const;
};

} //ns MATH
} //ns FI

#endif /* _PLANE_H */