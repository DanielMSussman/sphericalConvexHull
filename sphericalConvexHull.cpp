#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy
#include "convexHullCGAL.h"
#include "sphericalDomain.h"
using namespace std;


void printdVecVectorToScreen(vector<dVec> &a)
    {
    cout <<"{";
    for (int ii = 0; ii < a.size();++ii)
        {
        cout <<"{"<<a[ii].x[0]<<","<<a[ii].x[1]<<","<<a[ii].x[2] << "}";
        if(ii != a.size()-1)
            cout <<",";
        }
    cout<< "}"<<endl;
    }

void printIntVector(vector<int> &a)
    {
    cout <<"{";
    for (int ii = 0; ii < a.size();++ii)
        {
        cout << a[ii];
        if(ii != a.size()-1)
            cout <<",";
        }
    cout<< "}"<<endl;
    }

/*!
core of spherical vertexmodel
*/
int main(int argc, char*argv[])
{

    int nCells = 50;
    vector<dVec> cellCentroids;
    //for example: choose random cell positions and project them onto a sphere of given radius
    double radius = 3.7;
    sphericalDomain sphere(radius);
    for (int ii = 0; ii < nCells; ++ii)
        {
        dVec point;
        //choose a point uniformly on the sphere
        double u = drand48();
        double v = drand48();
        double phi = 2.0*PI*u;
        double theta = acos(2.0*v-1.);
        point.x[0] = radius * sin(theta)*cos(phi);
        point.x[1] = radius * sin(theta)*sin(phi);
        point.x[2] = radius * cos(theta);
 
        cellCentroids.push_back(point);
        };

    vector<int> cellCellNeighbors;
    vector<unsigned int> numberOfCellNeighbors;
    Index2D neighborIndexer;
    vector<dVec> vertexPositions;
    vector<int> vertexVertexNeighbors;
    vector<int> vertexCellNeighbors;
    vector<unsigned int> numberOfVertexNeighbors;
    Index2D vertexNeighborIndexer;
    convexHullCGALInterface convexHuller;
    convexHuller.sphericalConvexHullForVertexModel(cellCentroids,
                                                   nCells,
                                                   cellCellNeighbors,
                                                   numberOfCellNeighbors,
                                                   neighborIndexer,
                                                   vertexPositions,
                                                   vertexVertexNeighbors,
                                                   vertexCellNeighbors,
                                                   numberOfVertexNeighbors,
                                                   vertexNeighborIndexer);

   printdVecVectorToScreen(cellCentroids);
   printdVecVectorToScreen(vertexPositions);
    printIntVector(vertexVertexNeighbors);
    return 0;
};
