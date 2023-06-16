#include "std_include.h" // std library includes, definition of scalar, etc.. has a "using namespace std" in it, because I'm lazy
#include "convexHullCGAL.h"
#include "sphericalDomain.h"
using namespace std;

void printdVec(dVec &b)
    {
    cout <<"{"<<b.x[0]<<","<<b.x[1]<<","<<b.x[2] << "}";
    }

void readTextFile(string filename, vector<dVec> &a)
    {
    ifstream loadFile(filename);
    string fileLine;

    while(getline(loadFile,fileLine))
        {
        istringstream iss(fileLine);
        double x,y,z;
        dVec entry;
        //throw exception if the file isn't properly formatted
        if(!(iss >> x >>y>>z))
            throw std::exception();
        entry.x[0]=x;
        entry.x[1]=y;
        entry.x[2]=z;
        a.push_back(entry);
        }
    };

void printdVecVectorToScreen(vector<dVec> &a)
    {
    cout <<"{";
    for (int ii = 0; ii < a.size();++ii)
        {
        printdVec(a[ii]);
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
    /*
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
    */
    //load file
    readTextFile("../data/nuclei_centroids.txt",cellCentroids);
    nCells = cellCentroids.size();
    //assume that all points are on a sphere whose center is the origin and whose radius is given by the norm of the first element
    radius = sqrt(cellCentroids[0].x[0]*cellCentroids[0].x[0]+cellCentroids[0].x[1]*cellCentroids[0].x[1]+cellCentroids[0].x[2]*cellCentroids[0].x[2]);
    sphericalDomain sphere(radius);
    for (int ii = 0; ii < nCells; ++ii)
        sphere.putInBoxReal(cellCentroids[ii]);

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
