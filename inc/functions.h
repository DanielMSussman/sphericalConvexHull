#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "std_include.h"
#include <set>

#define HOSTDEVICE inline __attribute__((always_inline))

/*! \file functions.h */


//!The dot product between d-Dimensional vectors.
HOSTDEVICE scalar dot(const dVec &p1, const dVec &p2)
    {
    scalar ans = 0.0;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans+=p1.x[dd]*p2.x[dd];

    return ans;
    };

//! an integer to the dth power... the slow way
HOSTDEVICE int idPow(int i)
    {
    int ans = i;
    for (int dd = 1; dd < DIMENSION; ++dd)
        ans *= i;
    return ans;
    };

//!The dot product between d-Dimensional iVecs.
HOSTDEVICE int dot(const iVec &p1, const iVec &p2)
    {
    int ans = 0;
    for (int dd = 0; dd < DIMENSION; ++dd)
        ans+=p1.x[dd]*p2.x[dd];

    return ans;
    };

//!Only use if DIMENSIOn ==3
HOSTDEVICE dVec cross(dVec &p1, dVec &p2)
    {
    dVec ans;
    ans[0] = p1[1]*p2[2]-p1[2]*p2[1];
    ans[1] = p1[2]*p2[0]-p1[0]*p2[2];
    ans[2] = p1[0]*p2[1]-p1[1]*p2[0];
    return ans;
    }

//!The norm of a d-Dimensional vector
HOSTDEVICE scalar norm(const dVec &p)
    {
    return sqrt(dot(p,p));
    };

HOSTDEVICE void rodriguesRotation(dVec &v, dVec &k, scalar theta)
    {
    dVec ans(0);
    dVec cp = cross(k,v);
    ans = cos(theta)*v + sin(theta)*cp + (1.0-cos(theta))*dot(k,v)*k;
    v = ans;
    };

//!fit integers into non-negative domains
HOSTDEVICE int wrap(int x,int m)
    {
    int ans = x;
    if(x >= m)
        ans = x % m;
    while(ans <0)
        ans += m;
    return ans;
    }


//!compute the sign of a scalar, and return zero if x = 0
HOSTDEVICE int computeSign(scalar x)
    {
    return ((x>0)-(x<0));
    };

//!compute the sign of a scalar, and return zero if x = 0...but return a scalar to avoid expensive casts on the GPU
HOSTDEVICE scalar computeSignNoCast(scalar x)
    {
    if (x > 0.) return 1.0;
    if (x < 0.) return -1.0;
    if (x == 0.) return 0.;
    return 0.0;
    };


#undef HOSTDEVICE
#endif
