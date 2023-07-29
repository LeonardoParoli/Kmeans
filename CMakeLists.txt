cmake_minimum_required(VERSION 3.6)
project(Kmeans)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_executable(Kmeans
        main.cpp
        Kmeans/KmeansInitializer.cpp
        Kmeans/KmeansInitializer.h
        Kmeans/sequential/KmeansSequentialSolver.cpp
        Kmeans/sequential/KmeansSequentialSolver.h
        Kmeans/Point.h
        Kmeans/Kluster.h
        Kmeans/Kluster.cpp
        Kmeans/parallelOPENMP/KmeansParallelOMPSolver.cpp
        Kmeans/parallelOPENMP/KmeansParallelOMPSolver.h
        )