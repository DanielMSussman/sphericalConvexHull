#ifndef convexHullCGAL_H
#define convexHullCGAL_H

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

#include "indexer.h"

#include "dDimensionalVectorTypes.h"

using namespace std;
class convexHullCGALInterface
    {
    public:
        //!meant to be used as the vertex model initializer... computes DT, then infers vertex positions and neighbor relations
        void sphericalConvexHullForVertexModel(vector<dVec> cellpoints, int n, vector<int> &allNeighs,
                    vector<unsigned int> &numNeighs, Index2D &nidx, vector<dVec> &vertexPositions, vector<int> &vertexNeighs, 
                    vector<int> &vertexCellNeighs, vector<unsigned int> &numVertexNeighs, Index2D &vnidx);
    };
#endif
