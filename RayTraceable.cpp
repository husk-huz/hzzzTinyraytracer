#include "include/RayTraceable.h"
#include "include/RayTracer.h"
#include "include/model.h"

float RayTraceable::lightnessAdding(std::vector<RayTraceable*> tracing_objects, Model* light, float max_lightness, Vec3f point, Vec3f normal) 
{
    float lightness = 0.0f;
    float distance = (light->vert(0) - point).norm();
    Vec3f light_direction = (point - light->vert(0)).normalize();
    Ray light_ray(light->vert(0), light_direction);
    lightness = std::min(std::abs(max_lightness / (distance*distance)), max_lightness);
    lightness *= std::abs(normal*light_direction);
    for(RayTraceable* object: tracing_objects)
    {
        Vec3f intersectPoint;
        Vec3f objects_normal;
        if(object->SimpleIntersect(light_ray, intersectPoint, objects_normal))
        {
            float new_distance = (intersectPoint - light->vert(0)).norm();
            if(new_distance < distance - 1e-3)
            {
                return 0.0f;
            }
        }
    }

    return lightness;
}

RayTraceableModel::RayTraceableModel(Model* model) 
{
    this->model = model;
    this->type = RayTraceable::MODEL;
}

bool RayTraceableModel::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights) 
{
    res_lightness = 0;
    float min_t = std::numeric_limits<float>::max();
    Vec3f intersectPoint;
    int intersectFaceIdx = -1;
    Vec3f bc_coord;
    Vec3f render_bc_coord;
    Vec3f face_norm;
    if(this->GetIntersectFace(ray, intersectPoint, face_norm, intersectFaceIdx))
    {
        std::vector<int> face = model->face(intersectFaceIdx);
        std::vector<int> t_face = model->texture_face(intersectFaceIdx);
        Vec3f world_coords[3];
        Vec2f uv_coords[3];
        for(int j=0; j<3; j++)
        {
            Vec3f v = model->vert(face[j]);
            world_coords[j] = v;
            uv_coords[j] = model->uv(t_face[j]);
        }
        render_bc_coord = getBarycenterCoord(world_coords, intersectPoint);
        _intersectPoint = intersectPoint;
        float color_coord_x = uv_coords[0].x*render_bc_coord.raw[0] + uv_coords[1].x*render_bc_coord.raw[1] + uv_coords[2].x*render_bc_coord.raw[2];
        float color_coord_y = uv_coords[0].y*render_bc_coord.raw[0] + uv_coords[1].y*render_bc_coord.raw[1] + uv_coords[2].y*render_bc_coord.raw[2];
        res_color = model->getTextureColor(color_coord_x, color_coord_y);
        for(Model* light: lights)
        {
            res_lightness += this->lightnessAdding(this->father_tracer->scene, light, max_remain_lightness/lights.size(), intersectPoint, face_norm);
        }
        if(res_lightness > 1.0)
        {
            res_lightness = 1.0;
        }
        return true;
    }

    return false;
}

bool RayTraceableModel::SimpleIntersect(Ray ray, Vec3f& _intersectPoint, Vec3f& normal) 
{
    Vec3f bc_coord;
    Vec3f intersectPoint;
    float min_t = std::numeric_limits<float>::max();
    int intersectFaceIdx = -1;
    for(int i=0; i<model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f points[3];
        for(int j=0; j<3; j++)
        {
            points[j] = model->vert(face[j]);
        }
        if(ray.isIntersect(points, intersectPoint, bc_coord))
        {
            intersectFaceIdx = i;
            float t = (intersectPoint - ray.origin).norm();
            if(t < min_t)
            {
                min_t = t;
                _intersectPoint = intersectPoint;
                normal.x = bc_coord.raw[0]*model->norm(face[0]).x + bc_coord.raw[1]*model->norm(face[1]).x + bc_coord.raw[2]*model->norm(face[2]).x;
                normal.y = bc_coord.raw[0]*model->norm(face[0]).y + bc_coord.raw[1]*model->norm(face[1]).y + bc_coord.raw[2]*model->norm(face[2]).y;
                normal.z = bc_coord.raw[0]*model->norm(face[0]).z + bc_coord.raw[1]*model->norm(face[1]).z + bc_coord.raw[2]*model->norm(face[2]).z;
            }
        }
    }
    if(intersectFaceIdx == -1)
    {
        return false;
    }
    else 
    {
        return true;
    }
}

bool RayTraceableModel::GetIntersectFace(Ray ray, Vec3f& _intersectPoint, Vec3f& normal, int& face_index)
{
    Vec3f bc_coord;
    Vec3f intersectPoint;
    float min_t = std::numeric_limits<float>::max();
    int intersectFaceIdx = -1;
    for(int i=0; i<model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f points[3];
        for(int j=0; j<3; j++)
        {
            points[j] = model->vert(face[j]);
        }
        if(ray.isIntersect(points, intersectPoint, bc_coord))
        {
            intersectFaceIdx = i;
            float t = (intersectPoint - ray.origin).norm();
            if(t < min_t)
            {
                min_t = t;
                _intersectPoint = intersectPoint;
                face_index = intersectFaceIdx;
                normal.x = bc_coord.raw[0]*model->norm(face[0]).x + bc_coord.raw[1]*model->norm(face[1]).x + bc_coord.raw[2]*model->norm(face[2]).x;
                normal.y = bc_coord.raw[0]*model->norm(face[0]).y + bc_coord.raw[1]*model->norm(face[1]).y + bc_coord.raw[2]*model->norm(face[2]).y;
                normal.z = bc_coord.raw[0]*model->norm(face[0]).z + bc_coord.raw[1]*model->norm(face[1]).z + bc_coord.raw[2]*model->norm(face[2]).z;
            }
        }
    }
    if(intersectFaceIdx == -1)
    {
        return false;
    }
    else 
    {
        return true;
    }
}

RayTraceableSphere::RayTraceableSphere(Vec3f center, float r) 
{
    this->center = center;
    this->r = r;
    this->type = RayTraceable::SPHERE;
}

bool RayTraceableSphere::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights) 
{
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!this->SimpleIntersect(ray, intersectPoint, normal))
    {
        return false;
    }

    _intersectPoint = intersectPoint;
    
    for(Model* light: lights)
    {
        res_lightness += lightnessAdding((*this->father_tracer).scene, light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }

    return true;
}

bool RayTraceableSphere::SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal) 
{
    Vec3f origin = ray.origin;
    Vec3f direction = ray.direction;
    Vec3f oc = origin - center;
    float a = direction*direction;
    float b = oc*direction*2;
    float c = oc*oc - r*r;
    float delta = b*b - 4.0*a*c;
    if(delta < 0)
    {
        return false;
    }
    else
    {
        float t1 = (-b - sqrt(delta)) / (2.0*a);
        float t2 = (-b + sqrt(delta)) / (2.0*a);
        float t;
        if(t1 < 0 && t2 < 0)
        {
            return false;
        }
        else if( t1 > 0 && t2 > 0)
        {
            t = t1 < t2 ? t1 : t2;
        }
        else if( t1 > 0 )
        {
            t = t1;
        }
        else if(t2 > 0)
        {
            t = t2;
        }
        
        intersectPoint = origin + direction*t;
        normal = (intersectPoint - center).normalize();
        return true;
    }
}

RayTraceableInfPlane::RayTraceableInfPlane(float a, float b, float c, float d) 
{
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
    this->type = RayTraceable::INF_PLANE;
}

bool RayTraceableInfPlane::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!RayTraceableInfPlane::SimpleIntersect(ray, intersectPoint, normal))
    {
        return false;
    }
    
    _intersectPoint = intersectPoint;

    for(Model* light: lights)
    {
        res_lightness += lightnessAdding((*this->father_tracer).scene, light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }

    return true;
}

bool RayTraceableInfPlane::SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal) 
{
    Vec3f plane_norm = Vec3f(a, b, c);
    if(std::abs(a*ray.direction.x + b*ray.direction.y + c*ray.direction.z) <= 1e-10)
    {
        return false;
    }
    float t = -(a*ray.origin.x + b*ray.origin.y + c*ray.origin.z + d) / (a*ray.direction.x + b*ray.direction.y + c*ray.direction.z);
    if(t < 0)
    {
        return false;
    }
    intersectPoint = ray.origin + ray.direction*t;
    normal = plane_norm.normalize();
    return true;
}

RayTraceableBox::RayTraceableBox(std::vector<Vec3f> vertexs) 
{
    this->vertexes = vertexs;
    this->type = RayTraceable::BOX;
}

bool RayTraceableBox::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!RayTraceableBox::SimpleIntersect(ray, intersectPoint, normal))
    {
        return false;
    }
    
    _intersectPoint = intersectPoint;

    // res_lightness += 1;
    for(Model* light: lights)
    {
        res_lightness += lightnessAdding((*this->father_tracer).scene, light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }

    return true;
}

// 判断target是否在points构成的三角形内
bool isInface(Vec3f* points, Vec3f target)
{
    // Vec3f vec0 = points[0] - target;
    // Vec3f vec1 = points[0] - points[1];
    // Vec3f vec2 = points[0] - points[2];
    // Vec3f t0 = vec1^vec0;
    // Vec3f t1 = vec1^vec2;

    // vec0 = points[1] - target;
    // vec1 = points[1] - points[2];
    // vec2 = points[1] - points[0];
    // Vec3f t2 = vec1^vec0;
    // Vec3f t3 = vec1^vec2;

    // vec0 = points[2] - target;
    // vec1 = points[2] - points[1];
    // vec2 = points[2] - points[0];
    // Vec3f t4 = vec1^vec0;
    // Vec3f t5 = vec1^vec2;

    // if(t0*t1 > 0 
    // && t2*t3 > 0
    // && t4*t5 > 0)
    // {
    //     std::cout << "#" << target.x << " " << target.y << " " << target.z << std::endl;
    //     std::cout << "$0 " << points[0].x << " " << points[0].y << " " << points[0].z << std::endl;
    //     std::cout << "$1 " << points[1].x << " " << points[1].y << " " << points[1].z << std::endl;
    //     std::cout << "$2 " << points[2].x << " " << points[2].y << " " << points[2].z << std::endl;
    //     std::cout << "t0*t1: " << t0*t1 << std::endl;
    //     std::cout << "t2*t3: " << t2*t3 << std::endl;
    //     std::cout << "t4*t5: " << t4*t5 << std::endl;
    //     return true;
    // }
    // else 
    // {
    //     return false;
    // }
    Vec3f bc_coord = getBarycenterCoord(points, target);
    // std::cout << bc_coord.raw[0] << " " << bc_coord.raw[1] << " " << bc_coord.raw[2] << std::endl;
    if(bc_coord.raw[0] > -1e-6 && bc_coord.raw[0] < 1.0 + 1e-6 
    && bc_coord.raw[1] > -1e-6 && bc_coord.raw[1] < 1.0 + 1e-6
    && bc_coord.raw[2] > -1e-6 && bc_coord.raw[2] < 1.0 + 1e-6)
    {
        return true;
    }
    else 
    {
        return false;
    }
}

bool RayTraceableBox::SimpleIntersect(Ray ray, Vec3f& _intersectPoint, Vec3f& _normal) 
{
    Vec3f intersectPoint;
    Vec3f normal;
    Vec3f refer_norm;
    float min_t = std::numeric_limits<float>::max();
    bool hasIntersect = false;
    std::vector<Vec3f> points;
    points.resize(4);

    for(int i=0; i<6; i++)
    {
        if(i == 0)
        {
            points[0] = vertexes[0];
            points[1] = vertexes[2];
            points[2] = vertexes[3];
            points[3] = vertexes[1];
            refer_norm = vertexes[0] - vertexes[6];
        }
        else if(i == 1)
        {
            points[0] = vertexes[2];
            points[1] = vertexes[4];
            points[2] = vertexes[5];
            points[3] = vertexes[3];
            refer_norm = vertexes[2] - vertexes[0];
        }
        else if(i == 2)
        {
            points[0] = vertexes[4];
            points[1] = vertexes[5];
            points[2] = vertexes[7];
            points[3] = vertexes[6];
            refer_norm = vertexes[4] - vertexes[2];
        }
        else if(i == 3)
        {
            points[0] = vertexes[0];
            points[1] = vertexes[1];
            points[2] = vertexes[7];
            points[3] = vertexes[6];
            refer_norm = vertexes[0] - vertexes[2];
        }
        else if(i == 4)
        {
            points[0] = vertexes[0];
            points[1] = vertexes[2];
            points[2] = vertexes[4];
            points[3] = vertexes[6];
            refer_norm = vertexes[0] - vertexes[1];
        }
        else if(i == 5)
        {
            points[0] = vertexes[1];
            points[1] = vertexes[3];
            points[2] = vertexes[5];
            points[3] = vertexes[7];
            refer_norm = vertexes[1] - vertexes[0];
        }
        else 
        {}

        if(SimpleFaceIntersect(ray, points, refer_norm, intersectPoint, normal))
        {
            // std::cout << "yes " << std::endl;
            Vec3f vec_arr1[3] = {points[0], points[1], points[3]};
            Vec3f vec_arr2[3] = {points[2], points[3], points[1]};
            
            if(isInface(vec_arr1, intersectPoint) 
            || isInface(vec_arr2, intersectPoint))
            {
                float t = (intersectPoint - ray.origin).norm();
                if(t < min_t)
                {
                    min_t = t;
                    _intersectPoint = intersectPoint;
                    _normal = normal;
                    hasIntersect = true;
                }
            }
        }
    }

    return hasIntersect;
}

bool RayTraceableBox::SimpleFaceIntersect(Ray ray, std::vector<Vec3f> _points, Vec3f refer_norm, Vec3f& _intersectPoint, Vec3f& _normal)
{
    // Vec3f vec1 = _points[1] - _points[0];
    // Vec3f vec2 = _points[2] - _points[0];
    // Vec3f normal = vec1^vec2;
    // if(normal*refer_norm < 0)
    // {
    //     normal = Vec3f(0, 0, 0) - normal;
    // }
    // if(std::abs(normal*ray.direction) <= 1e-10)
    // {
    //     return false;
    // }

    // float t = (normal*(_points[0] - ray.origin)) / (normal*ray.direction);
    // if(t < 0)
    // {
    //     return false;
    // }
    // else 
    // {
    //     _intersectPoint = ray.origin + ray.direction*t;
    //     _normal = normal.normalize();
    //     return true;
    // }

    if(std::abs(refer_norm.x*ray.direction.x 
              + refer_norm.y*ray.direction.y 
              + refer_norm.z*ray.direction.z) <= 1e-10)
    {
        return false;
    }
    float d = -refer_norm.x*_points[0].x - refer_norm.y*_points[0].y - refer_norm.z*_points[0].z;
    float t = -(refer_norm.x*ray.origin.x + refer_norm.y*ray.origin.y + refer_norm.z*ray.origin.z + d) 
             / (refer_norm.x*ray.direction.x + refer_norm.y*ray.direction.y + refer_norm.z*ray.direction.z);
    if(t < 0)
    {
        return false;
    }
    _intersectPoint = ray.origin + ray.direction*t;
    _normal = refer_norm.normalize();
    return true;
}

void RayTraceableBox::rotate(float angle, Vec3f axis)
{
    for(int i=0; i<8; i++)
    {
        rotatePoint(this->vertexes[i], angle, axis);
    }
}

void RayTraceableBox::move(Vec3f offset)
{
    for(int i=0; i<8; i++)
    {
        this->vertexes[i] = this->vertexes[i] + offset;
    }
}

bool AbandonRayTraceable::getSphereIntersectPoint(Ray ray, Vec3f& intersectPoint, Vec3f& normal)
{
    Vec3f origin = ray.origin;
    Vec3f direction = ray.direction;
    Vec3f oc = origin - center;
    float a = direction*direction;
    float b = oc*direction*2;
    float c = oc*oc - r*r;
    float delta = b*b - 4.0*a*c;
    if(delta < 0)
    {
        // intersectPoint = Vec3f(0, 0, 0);
        return false;
    }
    else
    {
        float t1 = (-b - sqrt(delta)) / (2.0*a);
        float t2 = (-b + sqrt(delta)) / (2.0*a);
        float t;
        if(t1 < 0 && t2 < 0)
        {
            return false;
        }
        else if( t1 > 0 && t2 > 0)
        {
            t = t1 < t2 ? t1 : t2;
        }
        else if( t1 > 0 )
        {
            t = t1;
        }
        else if(t2 > 0)
        {
            t = t2;
        }
        
        intersectPoint = origin + direction*t;
        normal = (intersectPoint - center).normalize();
        return true;
    }
}

bool AbandonRayTraceable::getInfPlaneIntersectPoint(Ray ray, Vec3f& intersectPoint, Vec3f& normal)
{
    Vec3f plane_norm = Vec3f(a, b, c);
    if(std::abs(a*ray.direction.x + b*ray.direction.y + c*ray.direction.z) <= 1e-10)
    {
        return false;
    }
    // if(std::abs((a*ray.origin.x + b*ray.origin.y + c*ray.origin.z + d)) <= 1e-10)
    // {
    //     return false;
    // }
    // std::cout << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << std::endl;
    // std::cout << a << " " << b << " " << c << " " << d << std::endl;
    // std::cout << b*ray.origin.y << "#" << std::endl;
    // std::cout << b*ray.direction.y << std::endl;
    float t = -(a*ray.origin.x + b*ray.origin.y + c*ray.origin.z + d) / (a*ray.direction.x + b*ray.direction.y + c*ray.direction.z);
    // std::cout << t << std::endl;
    if(t < 0)
    {
        return false;
    }
    intersectPoint = ray.origin + ray.direction*t;
    normal = plane_norm.normalize();
    return true;
}

// ax+by+cz+d = 0
// k(d1, d2, d3) + (x0, y0, z0) 
// a(x0+kd1) + b(y0+kd2) + c(z0+kd3) + d = 0
// k = -(ax0+by0+cz0+d) / (ad1+bd2+cd3)

AbandonRayTraceable::AbandonRayTraceable(Vec3f center, float r)
{
    this->center = center;
    this->r = r;
    this->type = SPHERE;
}

AbandonRayTraceable::AbandonRayTraceable(Model* model)
{
    this->model = model;
    this->type = MODEL;
}

AbandonRayTraceable::AbandonRayTraceable(float a, float b, float c, float d)
{
    this->type = INF_PLANE;
    this->a = a;
    this->b = b;
    this->c = c;
    this->d = d;
}

float AbandonRayTraceable::lightnessAdding(Model* light, float max_lightness, Vec3f point)
{
    float lightness = 0.0f;
    float distance = (light->vert(0) - point).norm();
    lightness = max_lightness / (distance*distance + 1.0);
    return lightness;
}

// 计算亮度衰减的光源亮度单步叠加
float AbandonRayTraceable::lightnessAdding(Model* light, float max_lightness, Vec3f point, Vec3f normal)
{
    // float lightness = 1.0f;
    float lightness = 0.0f;
    float distance = (light->vert(0) - point).norm();
    Vec3f light_direction = (light->vert(0) - point).normalize();
    lightness = std::min(std::abs(max_lightness / (distance*distance)), max_lightness);
    lightness *= std::abs(normal*light_direction);
    if(lightness < 0)
    {
        lightness = 0;
    }
    return lightness;
}

float AbandonRayTraceable::lightnessAdding(std::vector<RayTraceable*> tracing_objects, Model* light, float max_lightness, Vec3f point, Vec3f normal)
{
    float lightness = 0.0f;
    float distance = (light->vert(0) - point).norm();
    Vec3f light_direction = (point - light->vert(0)).normalize();
    Ray light_ray(light->vert(0), light_direction);
    lightness = std::min(std::abs(max_lightness / (distance*distance)), max_lightness);
    lightness *= std::abs(normal*light_direction);
    for(RayTraceable* object: tracing_objects)
    {
        Vec3f intersectPoint;
        Vec3f objects_normal;
        if(object->SimpleIntersect(light_ray, intersectPoint, objects_normal))
        {
            float new_distance = (intersectPoint - light->vert(0)).norm();
            if(new_distance < distance - 1e-3)
            {
                return 0.0f;
            }
        }
    }

    return lightness;
}

void AbandonRayTraceable::ModelIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights)
{
    res_lightness = 0;
    float min_t = std::numeric_limits<float>::max();
    Vec3f intersectPoint;
    int intersectFaceIdx = -1;
    Vec3f bc_coord;
    Vec3f render_bc_coord;
    for(int i=0; i<model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f points[3];
        for(int j=0; j<3; j++)
        {
            points[j] = model->vert(face[j]);
        }
        if(ray.isIntersect(points, intersectPoint, bc_coord))
        {
            float t = (intersectPoint - ray.origin).norm();
            if(t < min_t)
            {
                min_t = t;
                intersectFaceIdx = i;
                render_bc_coord = getBarycenterCoord(points, intersectPoint);
            }
        }
    }
    if(intersectFaceIdx != -1)
    {
        std::vector<int> face = model->face(intersectFaceIdx);
        std::vector<int> t_face = model->texture_face(intersectFaceIdx);
        Vec3f world_coords[3];
        Vec2f uv_coords[3];
        for(int j=0; j<3; j++)
        {
            Vec3f v = model->vert(face[j]);
            world_coords[j] = v;
            uv_coords[j] = model->uv(t_face[j]);
        }
        float color_coord_x = uv_coords[0].x*render_bc_coord.raw[0] + uv_coords[1].x*render_bc_coord.raw[1] + uv_coords[2].x*render_bc_coord.raw[2];
        float color_coord_y = uv_coords[0].y*render_bc_coord.raw[0] + uv_coords[1].y*render_bc_coord.raw[1] + uv_coords[2].y*render_bc_coord.raw[2];
        res_color = model->getTextureColor(color_coord_x, color_coord_y);
        for(Model* light: lights)
        {
            res_lightness += lightnessAdding(light, max_remain_lightness/lights.size(), intersectPoint);
        }
        if(res_lightness > 1.0)
        {
            res_lightness = 1.0;
        }
    }
}

void AbandonRayTraceable::SphereIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights)
{
    res_color = TGAColor(255, 255, 255, 255);
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getSphereIntersectPoint(ray, intersectPoint, normal))
    {
        return;
    }
    
    for(Model* light: lights)
    {
        // res_lightness += lightnessAdding(light, max_remain_lightness/lights.size(), intersectPoint, normal);
        res_lightness += lightnessAdding(light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }
}

void AbandonRayTraceable::InfPlaneIntersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights)
{
    // res_color = TGAColor(255, 255, 255, 255);
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getInfPlaneIntersectPoint(ray, intersectPoint, normal))
    {
        return;
    }
    // std::cout << intersectPoint.x << " " << intersectPoint.y << "#" << std::endl;
    
    for(Model* light: lights)
    {
        // res_lightness += lightnessAdding(light, max_remain_lightness/lights.size(), intersectPoint, normal);
        res_lightness += lightnessAdding(light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }
}

bool AbandonRayTraceable::ModelIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    res_lightness = 0;
    float min_t = std::numeric_limits<float>::max();
    Vec3f intersectPoint;
    int intersectFaceIdx = -1;
    Vec3f bc_coord;
    Vec3f render_bc_coord;
    for(int i=0; i<model->nfaces(); i++)
    {
        std::vector<int> face = model->face(i);
        Vec3f points[3];
        for(int j=0; j<3; j++)
        {
            points[j] = model->vert(face[j]);
        }
        if(ray.isIntersect(points, intersectPoint, bc_coord))
        {
            float t = (intersectPoint - ray.origin).norm();
            if(t < min_t)
            {
                min_t = t;
                intersectFaceIdx = i;
                render_bc_coord = getBarycenterCoord(points, intersectPoint);
                _intersectPoint = intersectPoint;
            }
        }
    }
    if(intersectFaceIdx != -1)
    {
        std::vector<int> face = model->face(intersectFaceIdx);
        std::vector<int> t_face = model->texture_face(intersectFaceIdx);
        Vec3f world_coords[3];
        Vec2f uv_coords[3];
        for(int j=0; j<3; j++)
        {
            Vec3f v = model->vert(face[j]);
            world_coords[j] = v;
            uv_coords[j] = model->uv(t_face[j]);
        }
        float color_coord_x = uv_coords[0].x*render_bc_coord.raw[0] + uv_coords[1].x*render_bc_coord.raw[1] + uv_coords[2].x*render_bc_coord.raw[2];
        float color_coord_y = uv_coords[0].y*render_bc_coord.raw[0] + uv_coords[1].y*render_bc_coord.raw[1] + uv_coords[2].y*render_bc_coord.raw[2];
        res_color = model->getTextureColor(color_coord_x, color_coord_y);
        for(Model* light: lights)
        {
            res_lightness += lightnessAdding(light, max_remain_lightness/lights.size(), intersectPoint);
        }
        if(res_lightness > 1.0)
        {
            res_lightness = 1.0;
        }
        return true;
    }

    return false;
}

bool AbandonRayTraceable::SphereIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getSphereIntersectPoint(ray, intersectPoint, normal))
    {
        return false;
    }

    _intersectPoint = intersectPoint;
    
    for(Model* light: lights)
    {
        res_lightness += lightnessAdding((*this->father_tracer).scene, light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }

    return true;
}

bool AbandonRayTraceable::InfPlaneIntersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& _intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    res_color = this->color;
    res_lightness = 0;
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getInfPlaneIntersectPoint(ray, intersectPoint, normal))
    {
        return false;
    }
    
    _intersectPoint = intersectPoint;

    for(Model* light: lights)
    {
        res_lightness += lightnessAdding((*this->father_tracer).scene, light, max_remain_lightness, intersectPoint, normal);
    }

    if(res_lightness > 1.0)
    {
        res_lightness = 1.0;
    }

    return true;
}

void AbandonRayTraceable::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, float max_remain_lightness, std::vector<Model*> lights)
{
    if(this->type == RayTraceableType::MODEL)
    {
        ModelIntersect(ray, res_color, res_lightness, max_remain_lightness, lights);
    }
    else if(this->type == RayTraceableType::SPHERE)
    {
        SphereIntersect(ray, res_color, res_lightness, max_remain_lightness, lights);
    }
    else if(this->type == AbandonRayTraceable::INF_PLANE)
    {
        InfPlaneIntersect(ray, res_color, res_lightness, max_remain_lightness, lights);
    }
    else 
    {}
}

bool AbandonRayTraceable::Intersect(Ray ray, TGAColor& res_color, float& res_lightness, Vec3f& intersectPoint, float max_remain_lightness, std::vector<Model*> lights)
{
    if(this->type == RayTraceableType::MODEL)
    {
        return ModelIntersect(ray, res_color, res_lightness, intersectPoint, max_remain_lightness, lights);
    }
    else if(this->type == RayTraceableType::SPHERE)
    {
        res_color = this->color;
        return SphereIntersect(ray, res_color, res_lightness, intersectPoint, max_remain_lightness, lights);
    }
    else if(this->type == AbandonRayTraceable::INF_PLANE)
    {
        res_color = this->color;
        return InfPlaneIntersect(ray, res_color, res_lightness, intersectPoint, max_remain_lightness, lights);
    }
    else 
    {}
}

bool AbandonRayTraceable::SimpleIntersect(Ray ray, Vec3f& intersectPoint, Vec3f& normal)
{
    if(this->type == RayTraceableType::MODEL)
    {
        return SimpleModelIntersect(ray, intersectPoint);
    }
    else if(this->type == RayTraceableType::SPHERE)
    {
        return SimpleSphereIntersect(ray, intersectPoint);
    }
    else if(this->type == AbandonRayTraceable::INF_PLANE)
    {
        return SimpleInfPlaneIntersect(ray, intersectPoint);
    }
    else 
    {}
}

bool AbandonRayTraceable::SimpleModelIntersect(Ray ray, Vec3f& intersectPoint)
{
    return true;
}

bool AbandonRayTraceable::SimpleSphereIntersect(Ray ray, Vec3f& _intersectPoint)
{
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getSphereIntersectPoint(ray, intersectPoint, normal))
    {
        return false;
    }

    _intersectPoint = intersectPoint;
    
    return true;
}

bool AbandonRayTraceable::SimpleInfPlaneIntersect(Ray ray, Vec3f& _intersectPoint)
{
    Vec3f intersectPoint;
    Vec3f normal;
    if(!getInfPlaneIntersectPoint(ray, intersectPoint, normal))
    {
        return false;
    }
    
    _intersectPoint = intersectPoint;

    return true;
}



