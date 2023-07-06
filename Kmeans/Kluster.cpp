#include "Kluster.h"

Kluster::Kluster(){
    centroid = {-1,-1,-1};
}

Kluster::~Kluster() = default;

Point *Kluster::getCentroid() {
    return &centroid;
}

int Kluster::getSize() const{
    return points.size();
}

std::vector<Point*>* Kluster::getPoints() {
    return &points;
}

void Kluster::addPoint(Point *point) {
     points.push_back(point);
}

void Kluster::resetCluster() {
    points.clear();
    centroid = {-1,-1,-1};
}

void Kluster::clearPoints() {
    points.clear();
}

void Kluster::updateCentroid() {
    Point newCentroid = {0.0f,0.0f,0.0f};
    for(Point* point : points){
        newCentroid.x += point->x;
        newCentroid.y += point->y;
        newCentroid.z += point->z;
    }
    newCentroid.x = newCentroid.x / points.size();
    newCentroid.y = newCentroid.y / points.size();
    newCentroid.z = newCentroid.z / points.size();
    this->centroid = newCentroid;
}

void Kluster::setCentroid(Point *centroid) {
    this->centroid = *centroid;
}
