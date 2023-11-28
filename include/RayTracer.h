#pragma once
#include "Ray.h"
#include "RayTraceable.h"
#include "tgaimage.h"
#include <vector>
#include "Light.h"

class RayTracer
{
private:
    Ray* ray;
    
public:
    std::vector<RayTraceable*> scene ;
    std::vector<Light*> lights;

    RayTracer(){}
    RayTracer(Ray* _ray, std::vector<RayTraceable*> _scene) : ray(_ray), scene(_scene) {}

    void addSon(RayTraceable* son)
    {
        son->father_tracer = this;
        this->scene.push_back(son);
        return;
    }

    TGAColor& Trace();

    TGAColor& Trace(Ray* new_ray, Vec3f& color, int depth);
};