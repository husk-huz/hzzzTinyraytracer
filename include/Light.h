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

