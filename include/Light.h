#ifndef __LIGHT_H__
#define __LIGHT_H__
#include "geometry.h"

class Light 
{
public:
    Vec3f position;
    float intensity;
    
    Light(Vec3f position, float intensity)
    {
        this->position = position;
        this->intensity = intensity;
    }

    ~Light();
};

#endif