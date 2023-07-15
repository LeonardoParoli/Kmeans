//
// Created by posta on 11/07/2023.
//

#ifndef KMEANS_KMEANSPARALLELOMPSOLVER_H
#define KMEANS_KMEANSPARALLELOMPSOLVER_H


#include "../Point.h"
#include "../Kluster.h"

class KmeansParallelOMPSolver {
    private:
        Point *points;
        Point *selectedCentroids;
        Kluster *clusters;
        int numPoints;
        int numClusters;

    public:
        KmeansParallelOMPSolver(Point *workPoints, int numPoints, int numClusters, Point *selectedCentroids);
        ~KmeansParallelOMPSolver();
        void solve(bool printConsole);
        Point *getPoints();
        Point *getSelectedCentroids();
        Kluster *getClusters();
};


#endif //KMEANS_KMEANSPARALLELOMPSOLVER_H
