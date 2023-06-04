#ifndef _INTERSECTION_H
#define _INTERSECTION_H

#include "math/vec3.h"

namespace FI {
namespace MATH {

bool rayIntersectsTriangle(const vec3 &p, const vec3 &d, vec3 v0, vec3 v1, vec3 v2, float radius, float &t)
{

        vec3 e1,e2,h,s,q;
        float a,f,u,v;

        e1 = v1 - v0;
        e2 = v2 - v0;
        h = d | e2;
        a = e1.dot(h);

        if (a > -0.00001f && a < 0.00001f)
                return(false);

        f = 1.f/a;
        s = p - v0;
        u = f * s.dot(h);

        if (u < 0.0f || u > 1.0f)
                return(false);

        q = s | e1;
        v = f * d.dot(q);
        if (v < 0.0f || u + v > 1.0f)
                return(false);
        // at this stage we can compute t to find out where
        // the intersection point is on the line
        t = f * e2.dot(q);
        if (t > 0.00001f && t < radius) // ray intersection
                return(true);
        else // this means that there is a line intersection
                // but not a ray intersection
                return (false);
}

} //ns MATH
} //ns FI

#endif /*_INTERSECTION_H*/