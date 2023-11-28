#include "include/Ray.h"

Vec3f getBarycenterCoord(Vec2i* coord, Vec2i targetPoint)
{
    Vec3f barycenterCoord;
    barycenterCoord.raw[0] = (float(coord[1].y - coord[2].y)*float(targetPoint.x - coord[2].x) + float(coord[2].x - coord[1].x)*float(targetPoint.y - coord[2].y)) / (float(coord[1].y - coord[2].y)*float(coord[0].x - coord[2].x) + float(coord[2].x - coord[1].x)*float(coord[0].y - coord[2].y));
    barycenterCoord.raw[1] = (float(coord[2].y - coord[0].y)*float(targetPoint.x - coord[2].x) + float(coord[0].x - coord[2].x)*float(targetPoint.y - coord[2].y)) / (float(coord[1].y - coord[2].y)*float(coord[0].x - coord[2].x) + float(coord[2].x - coord[1].x)*float(coord[0].y - coord[2].y));
    barycenterCoord.raw[2] = 1 - barycenterCoord.raw[0] - barycenterCoord.raw[1];

    return barycenterCoord;
}

Vec3f getBarycenterCoord(Vec3f* coord, Vec3f targetPoint)
{
    Vec3f barycenterCoord;
    float f = ((coord[1].y - coord[2].y)*(coord[0].x - coord[2].x) + (coord[2].x - coord[1].x)*(coord[0].y - coord[2].y));
    if(std::abs(f) > 1e-10)
    {
        barycenterCoord.raw[0] = ((coord[1].y - coord[2].y)*(targetPoint.x - coord[2].x) + (coord[2].x - coord[1].x)*(targetPoint.y - coord[2].y)) / f;
        barycenterCoord.raw[1] = ((coord[2].y - coord[0].y)*(targetPoint.x - coord[2].x) + (coord[0].x - coord[2].x)*(targetPoint.y - coord[2].y)) / f;
        barycenterCoord.raw[2] = 1 - barycenterCoord.raw[0] - barycenterCoord.raw[1];
        return barycenterCoord;
    }
    
    f = ((coord[1].z - coord[2].z)*(coord[0].x - coord[2].x) + (coord[2].x - coord[1].x)*(coord[0].z - coord[2].z));
    if(std::abs(f) > 1e-10)
    {
        barycenterCoord.raw[0] = ((coord[1].z - coord[2].z)*(targetPoint.x - coord[2].x) + (coord[2].x - coord[1].x)*(targetPoint.z - coord[2].z)) / f;
        barycenterCoord.raw[2] = ((coord[2].z - coord[0].z)*(targetPoint.x - coord[2].x) + (coord[0].x - coord[2].x)*(targetPoint.z - coord[2].z)) / f;
        barycenterCoord.raw[1] = 1 - barycenterCoord.raw[0] - barycenterCoord.raw[2];
        return barycenterCoord;
    }

    f = ((coord[1].z - coord[2].z)*(coord[0].y - coord[2].y) + (coord[2].y - coord[1].y)*(coord[0].z - coord[2].z));
    if(std::abs(f) > 1e-10)
    {
        barycenterCoord.raw[1] = ((coord[1].z - coord[2].z)*(targetPoint.y - coord[2].y) + (coord[2].y - coord[1].y)*(targetPoint.z - coord[2].z)) / f;
        barycenterCoord.raw[2] = ((coord[2].z - coord[0].z)*(targetPoint.y - coord[2].y) + (coord[0].y - coord[2].y)*(targetPoint.z - coord[2].z)) / f;
        barycenterCoord.raw[0] = 1 - barycenterCoord.raw[1] - barycenterCoord.raw[2];
        return barycenterCoord;
    }

    return Vec3f(0, 0, 0);
}

void rotatePoint(Vec3f& target_point, float angle, Vec3f axis)
{
    float x = target_point.x;
    float y = target_point.y;
    float z = target_point.z;
    float u = axis.x;
    float v = axis.y;
    float w = axis.z;
    float cosA = cos(angle);
    float sinA = sin(angle);
    float xPrime = u*(u*x + v*y + w*z)*(1 - cosA) + x*cosA + (-w*y + v*z)*sinA;
    float yPrime = v*(u*x + v*y + w*z)*(1 - cosA) + y*cosA + (w*x - u*z)*sinA;
    float zPrime = w*(u*x + v*y + w*z)*(1 - cosA) + z*cosA + (-v*x + u*y)*sinA;
    target_point = Vec3f(xPrime, yPrime, zPrime);
}

Ray::Ray(Vec3f origin, Vec3f direction)
{
    this->origin = origin;
    this->direction = direction;
}

bool Ray::isIntersect(Vec3f* points, Vec3f& intersectPoint)
{
    // 这里明显不太对，之后要改一改

    Vec3f norm;
    norm.x = (points[1].y - points[0].y) * (points[2].z - points[0].z) - (points[1].z - points[0].z) * (points[2].y - points[0].y);
    norm.y = (points[1].z - points[0].z) * (points[2].x - points[0].x) - (points[1].x - points[0].x) * (points[2].z - points[0].z);
    norm.z = (points[1].x - points[0].x) * (points[2].y - points[0].y) - (points[1].y - points[0].y) * (points[2].x - points[0].x);

    float t = (norm.x * (points[0].x - this->origin.x) + norm.y * (points[0].y - this->origin.y) + norm.z * (points[0].z - this->origin.z)) / (norm.x * this->direction.x + norm.y * this->direction.y + norm.z * this->direction.z);
    Vec3f intersect = this->origin + this->direction * t;

    Vec3f bc_coord = getBarycenterCoord(points, intersect);

    if(bc_coord.x >= 0 && bc_coord.x <= 1 && bc_coord.y >= 0 && bc_coord.y <= 1 && bc_coord.z >= 0 && bc_coord.z <= 1 && t>0)
    {
        intersectPoint = intersect;
        return true;
    }
    else
    {
        return false;
    }

    return true;
}

bool Ray::isIntersect(Vec3f* points, Vec3f& intersectPoint, Vec3f& bc_coord)
{
    // 这里明显不太对，之后要改一改

    Vec3f norm;
    norm.x = (points[1].y - points[0].y) * (points[2].z - points[0].z) - (points[1].z - points[0].z) * (points[2].y - points[0].y);
    norm.y = (points[1].z - points[0].z) * (points[2].x - points[0].x) - (points[1].x - points[0].x) * (points[2].z - points[0].z);
    norm.z = (points[1].x - points[0].x) * (points[2].y - points[0].y) - (points[1].y - points[0].y) * (points[2].x - points[0].x);

    float t = (norm.x * (points[0].x - this->origin.x) + norm.y * (points[0].y - this->origin.y) + norm.z * (points[0].z - this->origin.z)) / (norm.x * this->direction.x + norm.y * this->direction.y + norm.z * this->direction.z);
    Vec3f intersect = this->origin + this->direction * t;

    Vec3f temp_bc_coord = getBarycenterCoord(points, intersect);

    if(temp_bc_coord.raw[0] >= 0 && temp_bc_coord.raw[0] <= 1 && temp_bc_coord.raw[1] >= 0 && temp_bc_coord.raw[1] <= 1 && temp_bc_coord.raw[2] >= 0 && temp_bc_coord.raw[2] <= 1 && t>0)
    {
        intersectPoint = intersect;
        bc_coord = temp_bc_coord;
        // bc_coord.raw[0] = temp_bc_coord.raw[0];
        // bc_coord.raw[1] = temp_bc_coord.raw[1];
        // bc_coord.raw[2] = temp_bc_coord.raw[2];        
        return true;
    }
    else
    {
        return false;
    }

    return true;
}

Ray::~Ray()
{
}
