This repository extracts out the core "generate the spherical convex hull of a point set" functionality from the curvedSpaceVertexModel repository (github.com/sussmanLab/curvedSpaceVertexModel)

It makes use of CGAL (https://www.cgal.org/); once that it downloaded and installed, this code can be run by:

cd build
cmake ..
make
./sphericalConvexHull.out
