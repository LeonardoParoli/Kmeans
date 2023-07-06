#include <random>
#include "KmeansInitializer.h"

KmeansInitializer::KmeansInitializer(int numPoints, int numClusters, double coordinateRange, double clusterRadius){
    this->numPoints = numPoints;
    this->numClusters = numClusters;
    points = new Point[numPoints];
    realCentroids = new Point[numClusters];

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(-clusterRadius, clusterRadius);


    int pointCursor = 0;
    int centroidCursor = 0;
    for (int i = 0; i < numClusters; ++i) {
        float centerX = distribution(gen) + coordinateRange / 2.0;
        float centerY = distribution(gen) + coordinateRange / 2.0;
        float centerZ = distribution(gen) + coordinateRange / 2.0;
        realCentroids[centroidCursor].x = centerX;
        realCentroids[centroidCursor].y = centerY;
        realCentroids[centroidCursor].z = centerZ;
        centroidCursor++;
        for (int j = 0; j < numPoints / numClusters; ++j) {
            points[pointCursor].x = std::max(0.0, std::min(coordinateRange,distribution(gen) + centerX));
            points[pointCursor].y = std::max(0.0, std::min(coordinateRange,distribution(gen) + centerY));
            points[pointCursor].z = std::max(0.0, std::min(coordinateRange,distribution(gen) + centerZ));
            pointCursor++;
        }
    }
}

KmeansInitializer::~KmeansInitializer() {
    delete[] points;
    delete[] realCentroids;
}

Point * KmeansInitializer::getPoints() const {
    return points;
}

Point * KmeansInitializer::getRealCentroids() const {
    return realCentroids;
}

int KmeansInitializer::getNumPoints() const{
    return numPoints;
}

int KmeansInitializer::getNumClusters() const{
    return numClusters;
}


