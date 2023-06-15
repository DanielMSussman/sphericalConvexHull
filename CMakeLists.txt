cmake_minimum_required(VERSION 3.9.0)
project(sphericalCells)

add_definitions(-DDIMENSION=3)

set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14")

find_package(CGAL REQUIRED)

add_subdirectory(model)

include_directories(inc  model)

set(myLibs
    model
    )

# list the names of cpp files corresponding to linked executables you'd like...NON-GUI
foreach(ARG sphericalConvexHull)
add_executable("${ARG}.out" "${ARG}.cpp" )
target_link_libraries("${ARG}.out" ${myLibs}
    CGAL::CGAL
    )
endforeach()
