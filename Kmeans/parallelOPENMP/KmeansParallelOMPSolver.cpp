#include <iostream>
#include <omp.h>
#include <random>
#include "KmeansParallelOMPSolver.h"

Point *KmeansParallelOMPSolver::getPoints() {
    return points;
}

Kluster *KmeansParallelOMPSolver::getClusters() {
    return clusters;
}

KmeansParallelOMPSolver::~KmeansParallelOMPSolver() {
    delete[] clusters;
}

KmeansParallelOMPSolver::KmeansParallelOMPSolver(Point* workPoints, int numPoints, int numClusters, Point *selectedCentroids) {
    this->points= workPoints;
    this->numPoints= numPoints;
    this->numClusters=numClusters;
    this->selectedCentroids = selectedCentroids;
    Kluster *tempClusters= new Kluster[numClusters];
    for(int i =0; i < numClusters; i++){
        tempClusters[i] = Kluster();
        Point centroid = {selectedCentroids[i].x, selectedCentroids[i].y, selectedCentroids[i].z};
        tempClusters[i].setCentroid(&centroid);
    }
    this->clusters = tempClusters;
}


void KmeansParallelOMPSolver::solve(bool printConsole) {
    //Calculate base SSE
    double maxSSE = 0.0;
    #pragma omp parallel for default(none) reduction(+:maxSSE)
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
    maxSSE = maxSSE/numPoints;
    if(printConsole){
        std::cout <<"Max SSE = " << maxSSE << "" << std::endl;
    }

    auto* assignments = new int[numPoints];
    for(int i = 0; i < numPoints; i++){
        assignments[i]=-1;
    }
    auto* clusterSizes = new int[numClusters];
    for(int i = 0; i<numClusters; i++){
        clusterSizes[i] = 0;
    }
    //Starting Kmeans
    double currentSSE = 0;
    double previousSSE = maxSSE;
    double threshold = 0.01;
    int iteration = 0;
    auto *currentCentroids = new Point[numClusters];
    #pragma omp parallel for shared(currentCentroids) default(none)
    for(int j=0; j<numClusters; j++){
        currentCentroids[j] = {selectedCentroids[j].x,selectedCentroids[j].y,selectedCentroids[j].z};
    }
    while (std::abs(previousSSE - currentSSE) >= threshold && iteration < 10000) {
        previousSSE = currentSSE;
        //clear clusters
        #pragma omp parallel for shared(clusterSizes) default(none)
        for(int i = 0; i<numClusters; i++){
            clusterSizes[i] = 0;
        }
        //assign points to nearest centroid
        #pragma omp parallel for shared(assignments, clusterSizes, currentCentroids) default(none)
        for (int i = 0; i < numPoints; i++) {
            Point point = points[i];
            double minDistance = std::numeric_limits<double>::max();
            int closestCluster = -1;
            Point centroid;
            for (int j = 0; j < numClusters; j++) {
                centroid = currentCentroids[j];
                double distance = point.calculateDistance(centroid);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestCluster = j;
                }
            }
            assignments[i] = closestCluster;
            #pragma omp atomic
            clusterSizes[closestCluster]++;
        }

        auto* newCentroids = new Point[numClusters];
        #pragma omp parallel  for shared(newCentroids) default(none)
        for (int i = 0; i < numClusters; i++) {
            newCentroids[i] = {0.0, 0.0, 0.0};
        }
        #pragma omp parallel default(none) shared(assignments, points, newCentroids)
        {
            auto* privateAccumulators = new Point[numClusters];

            #pragma omp for
            for (int i = 0; i < numClusters; i++) {
                privateAccumulators[i] = {0.0, 0.0, 0.0};
            }

            #pragma omp for
            for (int i = 0; i < numPoints; i++) {
                int clusterIndex = assignments[i];
                #pragma omp atomic
                privateAccumulators[clusterIndex].x += points[i].x;
                #pragma omp atomic
                privateAccumulators[clusterIndex].y += points[i].y;
                #pragma omp atomic
                privateAccumulators[clusterIndex].z += points[i].z;
            }
            #pragma omp critical
            {
                for (int i = 0; i < numClusters; i++) {
                    newCentroids[i].x += privateAccumulators[i].x;
                    newCentroids[i].y += privateAccumulators[i].y;
                    newCentroids[i].z += privateAccumulators[i].z;
                }
            }
            delete[] privateAccumulators;
        }

        #pragma omp parallel for shared(newCentroids,currentCentroids,clusterSizes) default(none)
        for(int i=0; i<numClusters; i++){
            if(clusterSizes[i]==0){
                Point newCentroid = {0.0f,0.0f,0.0f};
                std::random_device rd;
                std::mt19937 rng(rd());
                std::uniform_real_distribution<double> dist(0.0, 1000.0);
                newCentroid.x = dist(rng);
                newCentroid.y = dist(rng);
                newCentroid.z = dist(rng);
                currentCentroids[i] = newCentroid;
            }else{
                newCentroids[i].x /= clusterSizes[i];
                newCentroids[i].y /= clusterSizes[i];
                newCentroids[i].z /= clusterSizes[i];
                currentCentroids[i] = newCentroids[i];
            }
        }

        //Update currentSSE
        currentSSE = 0.0;
        #pragma omp parallel for shared(currentCentroids,assignments) default(none) reduction(+:currentSSE)
        for(int i=0; i <numPoints; i++){
            Point point = points[i];
            double distance = point.calculateDistance(currentCentroids[assignments[i]]);
            currentSSE += distance * distance;
        }
        currentSSE = currentSSE/numPoints;
        if(printConsole){
            std::cout <<"Current SSE = " << currentSSE << "" << std::endl;
        }
        iteration++;
    }

    #pragma omp parallel for shared(assignments) default(none)
    for(int i=0; i<numPoints; i++){
        auto *point = new Point(points[i].x, points[i].y, points[i].z);
        #pragma omp critical
        {
            clusters[assignments[i]].addPoint(point);
        }
    }
}

/*void KmeansParallelOMPSolver::solve(bool printConsole) {
    int maxThreads = omp_get_max_threads();

    //Calculate base SSE
    double maxSSE = 0.0;
    #pragma omp parallel for default(none) reduction(+:maxSSE)
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
    maxSSE = maxSSE/numPoints;
    if(printConsole){
        std::cout <<"Max SSE = " << maxSSE << "" << std::endl;
    }


    //preparing thread arrays
    std::vector<std::vector<std::vector<Point>>> threadClustersArray(maxThreads);
    for(int i=0; i<maxThreads; i++){
        threadClustersArray[i].resize(numClusters);
    }
    threadClustersArray.resize(maxThreads);
    std::vector<Point*> threadCentroidsArray(maxThreads);
    threadCentroidsArray.resize(maxThreads);
    #pragma omp parallel for shared(threadCentroidsArray,selectedCentroids,maxThreads,numClusters) default(none)
    for(int i=0; i < maxThreads; i++){
        threadCentroidsArray[i] = new Point[numClusters];
        for(int j=0; j<numClusters; j++){
            Point centroid = {selectedCentroids[j].x,selectedCentroids[j].y,selectedCentroids[j].z};
            threadCentroidsArray[i][j] = centroid;
        }
    }
    std::vector<std::vector<Point*>> threadPoints(maxThreads);
    threadPoints.resize(maxThreads);
    int itemsPerArray = numPoints / maxThreads;  // Number of items per array
    int remainingItems = numPoints % maxThreads; // Number of remaining items
    int currentIndex = 0;
    #pragma omp parallel for shared(maxThreads,threadPoints, currentIndex, itemsPerArray, remainingItems, numPoints) default(none)
    for (int i = 0; i < maxThreads; ++i) {
        int numItems = itemsPerArray;
        if (remainingItems > 0) {
            numItems++;
            remainingItems--;
        }
        threadPoints[i].resize(numItems);
        for (int j = 0; j < numItems; ++j) {
            threadPoints[i][j] = &points[currentIndex];
            currentIndex++;
        }
    }

    //Starting Kmeans
    double currentSSE = 0;
    double previousSSE = maxSSE;
    double threshold = 0.01;
    int iteration = 0;
    Point currentCentroids[numClusters];
    std::vector<std::vector<Point>> currentClusters(numClusters);
    currentClusters.resize(numClusters);
    #pragma omp parallel for shared(currentCentroids) default(none)
    for(int j=0; j<numClusters; j++){
        currentCentroids[j] = {selectedCentroids[j].x,selectedCentroids[j].y,selectedCentroids[j].z};
    }
    while (std::abs(previousSSE - currentSSE) >= threshold && iteration < 10000) {
        previousSSE = currentSSE;
        //clear clusters
        #pragma omp parallel for collapse(2) shared(maxThreads,numClusters,threadClustersArray,threadCentroidsArray,currentCentroids) default(none)
        for(int i = 0; i < maxThreads; i++){
            for(int j = 0; j < numClusters; j++){
                threadClustersArray[i][j].clear();
                threadCentroidsArray[i][j] = {currentCentroids[j].x, currentCentroids[j].y,currentCentroids[j].z};
            }
        }

        for(int i=0; i < numClusters; i++){
            currentClusters[i].clear();
        }

        //assign points to nearest centroid
        #pragma omp parallel for shared(maxThreads,threadPoints,threadCentroidsArray,threadClustersArray) default(none)
        for(int i= 0; i<maxThreads; i++){
            std::vector<Point*> pointsBlock = threadPoints[i];
            for(auto point : pointsBlock){
                double minDistance = std::numeric_limits<double>::max();
                int closestCluster = -1;
                for( int j = 0; j < numClusters; j++) {
                    Point centroid = threadCentroidsArray[i][j];
                    double distance = point->calculateDistance(centroid);
                    if (distance < minDistance) {
                        minDistance = distance;
                        closestCluster = j;
                    }
                }
                threadClustersArray[i][closestCluster].push_back(*point);
            }
        }

        //reduction by thread
        for(int i = 0; i < maxThreads; i++){
            for( int j = 0; j < numClusters; j++) {
                for (auto point: threadClustersArray[i][j]) {
                    currentClusters[j].push_back(point);
                }
            }
        }

        //Update centroids
        #pragma omp parallel for shared(currentClusters, currentCentroids) default(none)
        for(int j=0; j < numClusters; j++){
            int clusterSize = currentClusters[j].size();
            if(clusterSize > 0){
                Point newCentroid = {0.0f,0.0f,0.0f};
                for(Point point : currentClusters[j]){
                    newCentroid.x += point.x;
                    newCentroid.y += point.y;
                    newCentroid.z += point.z;
                }
                newCentroid.x = newCentroid.x / clusterSize;
                newCentroid.y = newCentroid.y / clusterSize;
                newCentroid.z = newCentroid.z / clusterSize;
                currentCentroids[j] = newCentroid;
            }else{
                Point newCentroid = {0.0f,0.0f,0.0f};
                std::random_device rd;
                std::mt19937 rng(rd());
                std::uniform_real_distribution<double> dist(0.0, 1000.0);
                newCentroid.x = dist(rng);
                newCentroid.y = dist(rng);
                newCentroid.z = dist(rng);
                currentCentroids[j] = newCentroid;
            }
        }

        //Update currentSSE
        currentSSE = 0.0;
        #pragma omp parallel for shared(currentClusters,currentCentroids) default(none) reduction(+:currentSSE)
        for(int i = 0; i < numClusters; i++){
            for(auto point: currentClusters[i]){
                double distance = point.calculateDistance(currentCentroids[i]);
                currentSSE += distance * distance;
            }
        }
        currentSSE = currentSSE/numPoints;
        if(printConsole){
            std::cout <<"Current SSE = " << currentSSE << "" << std::endl;
        }
        iteration++;
    }


    for(int i = 0; i < numClusters; i++){
        clusters[i].setCentroid(&currentCentroids[i]);
        int size = currentClusters[i].size();
        for(int j=0; j < size; j++ ){
            auto* point = new Point(currentClusters[i][j].x, currentClusters[i][j].y, currentClusters[i][j].z);
            clusters[i].addPoint(point);
        }
    }
}*/

Point *KmeansParallelOMPSolver::getSelectedCentroids() {
    return selectedCentroids;
}