#pragma once
#include "geometry.h"

Vec3f getBarycenterCoord(Vec2i* coord, Vec2i targetPoint);
Vec3f getBarycenterCoord(Vec3f* coord, Vec3f targetPoint);
void rotatePoint(Vec3f& target_point, float angle, Vec3f axis);

class Ray
{
public:
    Vec3f origin;
    Vec3f direction;

    bool isIntersect(Vec3f* points, Vec3f& intersectPoint);
    bool isIntersect(Vec3f* points, Vec3f& intersectPoint, Vec3f& bc_coord);
    Ray(Vec3f origin, Vec3f direction);
    ~Ray();
};