#ifndef KMEANSSEQUENTIAL_KMEANSSEQUENTIALSOLVER_H
#define KMEANSSEQUENTIAL_KMEANSSEQUENTIALSOLVER_H

#include <vector>
#include "Point.h"
#include "Kluster.h"

class KmeansSequentialSolver {
    private:
        Point *points;
        Point *selectedCentroids;
        Kluster *clusters;
        int numPoints;
        int numClusters;

    public:
        KmeansSequentialSolver(Point *workPoints, int numPoints, int numClusters, Point *selectedCentroids);
        ~KmeansSequentialSolver();
        void solve();
        Point *getPoints();
        Point *getSelectedCentroids();
        Kluster *getClusters();
};


#endif //KMEANSSEQUENTIAL_KMEANSSEQUENTIALSOLVER_H
