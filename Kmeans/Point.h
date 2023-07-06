#ifndef KMEANSSEQUENTIAL_POINT_H
#define KMEANSSEQUENTIAL_POINT_H

#include <cmath>

struct Point{
    float x;
    float y;
    float z;

    float calculateDistance(const Point b){
        float distance = 0.0;
        distance += pow(this->x - b.x, 2);
        distance += pow(this->y - b.y, 2);
        distance += pow(this->z - b.z, 2);
        return sqrt(distance);
    };
};

#endif //KMEANSSEQUENTIAL_POINT_H
