#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include "Kmeans/KmeansInitializer.h"
#include "Kmeans/KmeansSequentialSolver.h"

int main() {
    int numPoints = 5000;
    int numClusters = 5;
    double coordinateRange = 1000;
    double clusterRadius = 100;

    std::cout << "Initializing Kmeans dataset" << std::endl;
    KmeansInitializer initializer = KmeansInitializer(numPoints,numClusters,coordinateRange,clusterRadius);
    Point *points = initializer.getPoints();
    Point *realCentroids = initializer.getRealCentroids();

    std::string currentFilePath(__FILE__);
    std::size_t found = currentFilePath.find_last_of("/");
    std::string folderPath = currentFilePath.substr(0, found + 1);
    std::string filePath = folderPath + "/Initialization/initial_points.txt";
    std::ofstream outputInitialPointsFile(filePath);
    if (outputInitialPointsFile.is_open()) {
        for(int i=0; i < numPoints; i++ ){
            outputInitialPointsFile << "("<< points[i].x<<","<< points[i].y <<"," << points[i].z <<")" << std::endl;
        }
        outputInitialPointsFile.close();
    } else {
        std::cout << "Unable to open the file." << std::endl;
    }
    filePath = folderPath + "/Initialization/real_centroids.txt";
    std::ofstream outputRealCentroidsFile(filePath);
    if (outputRealCentroidsFile.is_open()) {
        for(int i=0; i < numClusters; i++ ){
            outputRealCentroidsFile << "("<< realCentroids[i].x<<","<< realCentroids[i].y <<"," << realCentroids[i].z <<")" << std::endl;
        }
        outputRealCentroidsFile.close();
    } else {
        std::cout << "Unable to open the file." << std::endl;
    }

    std::cout << "Selecting K centroids" << std::endl;
    Point workPoints[numPoints];
    Point selectedCentroids[numClusters];
    for (int i = 0; i < numPoints; i++) {
        workPoints[i] = points[i];
    }
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(workPoints, workPoints + numPoints, rng);
    for (int i = 0; i < numClusters; i++) {
        selectedCentroids[i] = workPoints[i];
    }

    std::cout << "Running K-Means..." << std::endl;
    KmeansSequentialSolver solver = KmeansSequentialSolver(workPoints,numPoints,numClusters,selectedCentroids);
    solver.solve();
    Kluster *clusters = solver.getClusters();

    std::cout << "Clustering done, saving results" << std::endl;
    filePath = folderPath + "/Results/clustered_points.txt";
    std::ofstream outputResultClustersFile(filePath);
    if (outputResultClustersFile.is_open()) {
        for(int i=0; i < numClusters; i++ ){
            Kluster cluster = clusters[i];
            Point centroid = *cluster.getCentroid();
            outputResultClustersFile << "Cluster " << i << "{" <<std::endl;
            outputResultClustersFile << "[" << centroid.x <<","<< centroid.y <<"," << centroid.z << "]" <<std::endl;
            outputResultClustersFile << "("<< points[i].x <<","<< points[i].y <<"," << points[i].z <<")" << std::endl;
            outputResultClustersFile << "}" <<std::endl;
        }
        outputResultClustersFile.close();
    } else {
        std::cout << "Unable to open the file." << std::endl;
    }

    std::cout << "Results saved, terminating..." << std::endl;

    return 0;
}
