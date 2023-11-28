#pragma once
#include "Ray.h"
#include "geometry.h"
#include "model.h"
#include <vector>
#include "Light.h"

class RayTracer;
class RayTraceable;

class RayTraceable
{
public:
    RayTracer* father_tracer;
    TGAColor color;
    enum RayTraceableType
    {
        NULL_TYPE,
        MODEL,
        SPHERE,
        INF_PLANE,
        BOX
    };
    RayTraceableType type;
    bool is_light = false;
    virtual bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights){};
    virtual bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Light*> lights){};
    virtual bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal){};
    float lightnessAdding(std::vector<RayTraceable*> tracing_objects, Model* light, float max_lightness, Vec3f point, Vec3f normal);
};

class RayTraceableModel : public RayTraceable
{
public:
    Model* model;
    RayTraceableModel(Model* model);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Light*> lights);
    bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
    bool GetIntersectFace(Ray ray, Vec3f& intersectPoint, Vec3f& normal, int& face_index);
};

class RayTraceableSphere : public RayTraceable
{
public:
    float r;
    Vec3f center;
    RayTraceableSphere(Vec3f center, float r);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Light*> lights);
    bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
};

class RayTraceableInfPlane : public RayTraceable
{
public:
    float a, b, c, d;
    RayTraceableInfPlane(float a, float b, float c, float d);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Light*> lights);
    bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
};

class RayTraceableBox : public RayTraceable
{
public:
/*
             0 ------- 1
             / |     /|
            /  |    / |
            2-------3 |
            | /6   |  /7
            |/     | /
            4-------5
*/
    std::vector<Vec3f> vertexes;
    RayTraceableBox(std::vector<Vec3f> vertexes);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Light*> lights);
    bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
    static bool SimpleFaceIntersect(Ray ray, std::vector<Vec3f> _points, Vec3f refer_norm, Vec3f& _intersectPoint, Vec3f& _normal);
    void rotate(float angle, Vec3f axis);
    void move(Vec3f offset);
};


class AbandonRayTraceable
{
private:
    bool getSphereIntersectPoint(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
    bool getInfPlaneIntersectPoint(Ray ray, Vec3f& intersectPoint, Vec3f& normal);

public:
    RayTracer* father_tracer;
    enum RayTraceableType
    {
        NULL_TYPE,
        MODEL,
        SPHERE,
        INF_PLANE
    };
    RayTraceableType type;
    TGAColor color;

    Model* model;

    float r;
    Vec3f center;

    float a, b, c, d;

    AbandonRayTraceable(float a, float b, float c, float d);
    AbandonRayTraceable(Vec3f center, float r);
    AbandonRayTraceable(Model* model);
    static float lightnessAdding(Model* light, float max_lightness, Vec3f point);
    static float lightnessAdding(Model* light, float max_lightness, Vec3f point, Vec3f normal);
    void ModelIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights);
    void SphereIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights);
    void InfPlaneIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights);
    bool ModelIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool SphereIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool InfPlaneIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool SimpleModelIntersect(Ray ray, Vec3f& intersectPoint);
    bool SimpleSphereIntersect(Ray ray, Vec3f& intersectPoint);
    bool SimpleInfPlaneIntersect(Ray ray, Vec3f& intersectPoint);
    void Intersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights);
    bool Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights);
    bool SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal);
    float lightnessAdding(std::vector<RayTraceable*> tracing_objects, Model* light, float max_lightness, Vec3f point, Vec3f normal);
};