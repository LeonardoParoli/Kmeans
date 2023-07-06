#include <iostream>
#include "KmeansSequentialSolver.h"

Point *KmeansSequentialSolver::getPoints() {
    return points;
}

Kluster *KmeansSequentialSolver::getClusters() {
    return clusters;
}

KmeansSequentialSolver::~KmeansSequentialSolver() {
    delete[] points;
    delete[] clusters;
    delete[] selectedCentroids;
}

KmeansSequentialSolver::KmeansSequentialSolver(Point *workPoints, int numPoints, int numClusters, Point *selectedCentroids) {
    this->points= workPoints;
    this->numPoints= numPoints;
    this->numClusters=numClusters;
    this->selectedCentroids = selectedCentroids;
    Kluster clusters[numClusters];
    for(int i =0; i < numClusters; i++){
        clusters[i] = Kluster();
        Point centroid = {selectedCentroids[i].x, selectedCentroids[i].y, selectedCentroids[i].z};
        clusters[i].setCentroid(&centroid);
    }
    this->clusters = clusters;
}

void KmeansSequentialSolver::solve() {
    //Calculate base SSE
    double maxSSE = 0.0;
    for (int i = 0; i < numPoints; i++) {
        Point point = points[i];
        double minDistance = std::numeric_limits<double>::max();
        for (int k = 0; k < numClusters; k++) {
            Point centroid = selectedCentroids[k];
            double distance = point.calculateDistance(centroid);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
        maxSSE += minDistance * minDistance;
    }
    double currentSSE = maxSSE;
    std::cout <<"Max SSE = " << maxSSE << "" << std::endl;

    // termination condition: if current SSE is < 1% maxSSE or iteration > 100k
    int iteration = 0;
    while (currentSSE >= (maxSSE / 100) && iteration < 100000) {
        //clear clusters
        for(int i = 0; i < numClusters; i++){
            clusters[i].clearPoints();
        }

        //assign points to nearest centroid
        for (int i = 0; i < numPoints; i++) {
            Point point = points[i];
            double minDistance = std::numeric_limits<double>::max();
            int closestCluster = -1;
            for( int j = 0; j < numClusters; j++){
                Point *centroid = clusters[j].getCentroid();
                double distance = point.calculateDistance(*centroid);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestCluster = j;
                }
            }
            clusters[closestCluster].addPoint(&point);
        }

        //Update centroids
        for(int i = 0; i < numClusters; i++){
            clusters[i].updateCentroid();
        }

        //Update currentSSE
        currentSSE = 0.0;
        for (int i = 0; i < numClusters; i++) {
            Kluster cluster = clusters[i];
            std::vector<Point*> *clusterPoints = cluster.getPoints();
            Point *centroid = cluster.getCentroid();
            for (Point *point: *clusterPoints) {
                double distance = point->calculateDistance(*centroid);
                currentSSE += distance * distance;
            }
        }

        //updating iteration
        iteration++;
        std::cout <<"Current SSE = " << currentSSE << "" << std::endl;
    }
}

Point *KmeansSequentialSolver::getSelectedCentroids() {
    return selectedCentroids;
}
