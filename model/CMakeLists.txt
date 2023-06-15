if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUGFLAGUP)
endif(${CMAKE_BUILD_TYPE} MATCHES "Debug")

include_directories(${CMAKE_SOURCE_DIR}/inc)
add_library(model
#    simpleModel.cpp
#    sphericalModel.cpp
#    sphericalVoronoi.cpp
    convexHullCGAL.cpp
#    sphericalVertexModel.cpp
    )
#target_compile_options(model PRIVATE -frounding-math)
