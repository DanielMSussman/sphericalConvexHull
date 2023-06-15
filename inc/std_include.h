#ifndef STDINCLUDE_H
#define STDINCLUDE_H

/*! \file std_include.h
a file to be included all the time... carries with it things DMS often uses.
It includes some handy debugging / testing functions, and includes too many
standard library headers
It also defines scalars as either floats or doubles, depending on
how the program is compiled
*/

#define HOSTDEVICE inline __attribute__((always_inline))
#define THRESHOLD 1e-18
#define EPSILON 1e-18

#include <cmath>
#include <algorithm>
#include <memory>
#include <ctype.h>
#include <random>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/time.h>
#include <string.h>
#include <stdexcept>
#include <cassert>

using namespace std;


//double precision value of pi
#define PI 3.141592653589793115997963468544185161590576171875

//decide whether to compute everything in floating point or double precision
#ifndef SCALARFLOAT
//double variables types
#define scalar double
//trig and special funtions
#define Cos cos
#define Sin sin
#define Floor floor
#define Ceil ceil

#else
//floats
#define scalar float
#define Cos cosf
#define Sin sinf
#define Floor floorf
#define Ceil ceilf
#endif

//!for cuda memory alignment, it is convenient to define a class which holds sets of four positions specified by their (theta, phi) coordinates
class quadAngularPosition
    {
    public:
        HOSTDEVICE quadAngularPosition(){};
        HOSTDEVICE quadAngularPosition(const double value)
            {
            for (int dd = 0; dd < 8; ++dd)
                x[dd] = value;
            };
        HOSTDEVICE quadAngularPosition(const quadAngularPosition &other)
            {
            for (int dd = 0; dd < 8; ++dd)
                x[dd] = other.x[dd];
            };
        double x[8];
        HOSTDEVICE double& operator[](int i){return x[i];};
        //mutating operators
        HOSTDEVICE quadAngularPosition& operator=(const quadAngularPosition &other)
            {
            for (int dd = 0; dd < 8; ++dd)
                this->x[dd] = other.x[dd];
            return *this;
            }
        HOSTDEVICE quadAngularPosition& operator-=(const quadAngularPosition &other)
            {
            for (int dd = 0; dd < 8; ++dd)
                this->x[dd] -= other.x[dd];
            return *this;
            }
        HOSTDEVICE quadAngularPosition& operator+=(const quadAngularPosition &other)
            {
            for (int dd = 0; dd < 8; ++dd)
                this->x[dd] += other.x[dd];
            return *this;
            }
    };

#include "dDimensionalVectorTypes.h"

//!Report somewhere that code needs to be written
static void unwrittenCode(const char *message, const char *file, int line)
    {
    printf("\nCode unwritten (file %s; line %d)\nMessage: %s\n",file,line,message);
    throw std::exception();
    }

//!A utility function for checking if a file exists
inline bool fileExists(const std::string& name)
    {
    ifstream f(name.c_str());
    return f.good();
    }
//A macro to say code needs to be written
#define UNWRITTENCODE(message) (unwrittenCode(message,__FILE__,__LINE__))
//spot-checking of code for debugging
#define DEBUGCODEHELPER printf("\nReached: file %s at line %d\n",__FILE__,__LINE__);

#undef HOSTDEVICE
#endif
