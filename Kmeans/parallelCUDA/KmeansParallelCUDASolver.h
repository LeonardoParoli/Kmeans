//
// Created by posta on 11/07/2023.
//

#ifndef KMEANS_KMEANSPARALLELCUDASOLVER_H
#define KMEANS_KMEANSPARALLELCUDASOLVER_H


#include "../Point.h"
#include "../Kluster.h"

class KmeansParallelCUDASolver {
    private:
        Point *points;
        Point *selectedCentroids;
        Kluster *clusters;
        int numPoints;
        int numClusters;

    public:
        KmeansParallelCUDASolver(Point *workPoints, int numPoints, int numClusters, Point *selectedCentroids);
        ~KmeansParallelCUDASolver();
        void solve();
        Point *getPoints();
        Point *getSelectedCentroids();
        Kluster *getClusters();
};


#endif //KMEANS_KMEANSPARALLELCUDASOLVER_H
