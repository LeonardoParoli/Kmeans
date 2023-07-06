cmake_minimum_required(VERSION 3.25)
project(KmeansSequential)

set(CMAKE_CXX_STANDARD 17)

add_executable(KmeansSequential main.cpp Kmeans/KmeansInitializer.cpp Kmeans/KmeansInitializer.h
        Kmeans/KmeansSequentialSolver.cpp Kmeans/KmeansSequentialSolver.h Kmeans/Point.h
        Kmeans/Kluster.h Kmeans/Kluster.cpp)
