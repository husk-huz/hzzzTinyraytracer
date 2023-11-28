#include "include/config.h"
#include "include/HzzzImage.h"
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include "include/RayTraceable.h"
#include "include/RayTracer.h"

TGAColor RandomColor()
{
    TGAColor color;
    color.r = rand() % 255;
    color.g = rand() % 255;
    color.b = rand() % 255;
    color.a = 255;
    return color;
}

HzzzImage::HzzzImage(int width, int height)
{
    this->width = width;
    this->height = height;
    this->image = new char[width * height * 4];
    this->zbuffer = new float[width * height];
    for(int i=0; i<width*height; i++)
    {
        this->zbuffer[i] = std::numeric_limits<float>::max();
    }
}

bool HzzzImage::set(int x, int y, const TGAColor &c)
{
    if (!this->image || x<0 || y<0 || x>=width || y>=height) 
    {
		return false;
	}
    int idx = x + y * width;
    image[idx*4    ] = c.r;
    image[idx*4 + 1] = c.g;
    image[idx*4 + 2] = c.b;
    image[idx*4 + 3] = c.a;

	return true;
}

bool HzzzImage::set(int x, int y, float z, const TGAColor &c)
{
    if (!this->image || x<0 || y<0 || x>=width || y>=height || c.val == -1) 
    {
		return false;
	}
    int idx = x + y * width;
    if (zbuffer[idx] > z) 
    {
        zbuffer[idx] = z;
        image[idx*4    ] = c.r;
        image[idx*4 + 1] = c.g;
        image[idx*4 + 2] = c.b;
        image[idx*4 + 3] = c.a;
    }
	return true;
}

bool HzzzImage::setCamera(Vec3f camera)
{
    this->camera = camera;
    return true;
}

bool HzzzImage::setViewDirection(Vec3f view_direction)
{
    this->view_direction = view_direction;
    return true;
}

Vec3f HzzzImage::getCamera()
{
    return this->camera;
}

Vec3f HzzzImage::getViewDirection()
{
    return this->view_direction;
}

void HzzzImage::Line(Vec2i p0, Vec2i p1, TGAColor color)
{
    bool steep = false;
    if (std::abs(p0.x-p1.x)<std::abs(p0.y-p1.y)) 
    {
        std::swap(p0.x, p0.y);
        std::swap(p1.x, p1.y);
        steep = true;
    }
    if (p0.x>p1.x) 
    {
        std::swap(p0.x, p1.x);
        std::swap(p0.y, p1.y);
    }

    for (int x=p0.x; x<=p1.x; x++) 
    {
        float t = (x-p0.x)/(float)(p1.x-p0.x);
        int y = p0.y*(1.-t) + p1.y*t;
        if (steep) 
        {
            this->set(y, x, color);
        } 
        else 
        {
            this->set(x, y, color);
        }
    }
}

void HzzzImage::Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, TGAImage& texture)
{
    // TGAColor color = RandomColor();
    Vec2i v3, v1, v2;
    v1 = this->toImageCoord(p1);
    v2 = this->toImageCoord(p2);
    v3 = this->toImageCoord(p3);

#ifdef DEBUG
    std::cout << uv[0].x << " " << uv[0].y << std::endl;
    std::cout << uv[1].x << " " << uv[1].y << std::endl;
    std::cout << uv[2].x << " " << uv[2].y << std::endl;
    std::cout << iuv0.x << " " << iuv0.y << std::endl;
    std::cout << iuv1.x << " " << iuv1.y << std::endl;
    std::cout << iuv2.x << " " << iuv2.y << std::endl;
#endif

    if (v1.y==v2.y && v1.y==v3.y) return; // I dont care about degenerate triangles 
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
    if (v1.y>v2.y)
    {   
        std::swap(v1, v2);
        std::swap(t1, t2);
        std::swap(p1, p2);
    } 
    if (v1.y>v3.y)
    { 
        std::swap(v1, v3);
        std::swap(t1, t3);
        std::swap(p1, p3);
    } 
    if (v2.y>v3.y)
    { 
        std::swap(v2, v3);
        std::swap(t2, t3);
        std::swap(p2, p3);
    } 
    int total_height = v3.y-v1.y; 
    Vec2i vertex[3] = {v1, v2, v3};
    // Vec2i iuv0, iuv1, iuv2;
    // iuv0.x = t1.x*texture.get_width();
    // iuv0.y = t1.y*texture.get_height();
    // iuv1.x = t2.x*texture.get_width();
    // iuv1.y = t2.y*texture.get_height();
    // iuv2.x = t3.x*texture.get_width();
    // iuv2.y = t3.y*texture.get_height();
    // Vec2i vec1 = iuv1 - iuv0;
    // Vec2i vec2 = iuv2 - iuv0;

    for (int i=0; i<total_height; i++) 
    { 
        bool second_half = i>v2.y-v1.y || v2.y==v1.y; 
        int segment_height = second_half ? v3.y-v2.y : v2.y-v1.y; 
        float alpha = (float)i/total_height; 
        float beta  = (float)(i-(second_half ? v2.y-v1.y : 0))/segment_height; // be careful: with above conditions no division by zero here 
        Vec2i A =               v1 + (v3-v1)*alpha; 
        Vec2i B = second_half ? v2 + (v3-v2)*beta : v1 + (v2-v1)*beta; 
        if (A.x>B.x)
        {
            std::swap(A, B);
        } 

        // Vec3 norm_vec = 

        for (int j=A.x; j<=B.x; j++) 
        { 
            // int uv_x;
            // int uv_y;
            // float denominator = (v3.x - v1.x)*(v2.y - v1.y) - (v2.x - v1.x)*(v3.y - v1.y);
            // float u, v;
            // u = ((j-v1.x)*(v2.y-v1.y)*1. - (i)*(v2.x-v1.x)*1.)/denominator;
            // v = ((v3.x-v1.x)*(i)*1. - (j-v1.x)*(v3.y-v1.y)*1.)/denominator;
            // uv_x = u*vec2.x + v*vec1.x + iuv0.x;
            // uv_y = u*vec2.y + v*vec1.y + iuv0.y;

            Vec3f bc_coord = getBarycenterCoord(vertex, Vec2i(j, v1.y+i));
            // if(bc_coord.raw[0] < 0 || bc_coord.raw[1] < 0 || bc_coord.raw[2] < 0)
            // {
            //     continue;
            // }
            // std::cout << "coord: " << bc_coord.raw[0] << " " << bc_coord.raw[1] << " " << bc_coord.raw[2] << std::endl;
            // std::cout << "this: " << j << " " << v1.y+i << std::endl;
            // std::cout << "coord sum: " << (bc_coord.raw[0]*v1.x + bc_coord.raw[1]*v2.x + bc_coord.raw[2]*v3.x) << "#" 
            //           << (bc_coord.raw[0]*v1.y + bc_coord.raw[1]*v2.y + bc_coord.raw[2]*v3.y) << std::endl;
            float z = bc_coord.raw[0]*p1.z + bc_coord.raw[1]*p2.z + bc_coord.raw[2]*p3.z;
            float uv_x = bc_coord.raw[0]*t1.x + bc_coord.raw[1]*t2.x + bc_coord.raw[2]*t3.x;
            float uv_y = bc_coord.raw[0]*t1.y + bc_coord.raw[1]*t2.y + bc_coord.raw[2]*t3.y;
            
            // std::cout << "uv_x: " << uv_x << " uv_y: " << uv_y << std::endl;
            // this->set(j, v1.y+i, z, texture.get(uv_y*texture.get_width(), uv_x*texture.get_height()));
            // this->set(j, v1.y+i, z, texture.get(uv_x, uv_y));
            this->set(j, v1.y+i, z, texture.get(uv_x, uv_y));
            // this->set(j, v1.y+i, z, color);
        } 
    } 

    return;
}

void HzzzImage::Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, TGAImage& texture, float depth)
{
    Vec2i v3, v1, v2;
    v1 = this->toImageCoord(p1);
    v2 = this->toImageCoord(p2);
    v3 = this->toImageCoord(p3);

    if (v1.y==v2.y && v1.y==v3.y) return; 
    if (v1.y>v2.y)
    {   
        std::swap(v1, v2);
        std::swap(t1, t2);
        std::swap(p1, p2);
    } 
    if (v1.y>v3.y)
    { 
        std::swap(v1, v3);
        std::swap(t1, t3);
        std::swap(p1, p3);
    } 
    if (v2.y>v3.y)
    { 
        std::swap(v2, v3);
        std::swap(t2, t3);
        std::swap(p2, p3);
    } 
    int total_height = v3.y-v1.y; 
    Vec2i vertex[3] = {v1, v2, v3};

    for (int i=0; i<total_height; i++) 
    { 
        bool second_half = i>v2.y-v1.y || v2.y==v1.y; 
        int segment_height = second_half ? v3.y-v2.y : v2.y-v1.y; 
        float alpha = (float)i/total_height; 
        float beta  = (float)(i-(second_half ? v2.y-v1.y : 0))/segment_height;  
        Vec2i A =               v1 + (v3-v1)*alpha; 
        Vec2i B = second_half ? v2 + (v3-v2)*beta : v1 + (v2-v1)*beta; 
        if (A.x>B.x)
        {
            std::swap(A, B);
        } 

        for (int j=A.x; j<=B.x; j++) 
        { 
            Vec3f bc_coord = getBarycenterCoord(vertex, Vec2i(j, v1.y+i));
            float z = bc_coord.raw[0]*p1.z + bc_coord.raw[1]*p2.z + bc_coord.raw[2]*p3.z;
            float uv_x = bc_coord.raw[0]*t1.x + bc_coord.raw[1]*t2.x + bc_coord.raw[2]*t3.x;
            float uv_y = bc_coord.raw[0]*t1.y + bc_coord.raw[1]*t2.y + bc_coord.raw[2]*t3.y;
            
            this->set(j, v1.y+i, z, texture.get(uv_x, uv_y) * depth);
        } 
    } 

    return;
}

void HzzzImage::Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, Model* model, float* depths)
{
    Vec2i v3, v1, v2;
    v1 = this->toImageCoord(p1);
    v2 = this->toImageCoord(p2);
    v3 = this->toImageCoord(p3);

    if (v1.y==v2.y && v1.y==v3.y) return; 
    if (v1.y>v2.y)
    {   
        std::swap(v1, v2);
        std::swap(t1, t2);
        std::swap(p1, p2);
        std::swap(depths[0], depths[1]);
    } 
    if (v1.y>v3.y)
    { 
        std::swap(v1, v3);
        std::swap(t1, t3);
        std::swap(p1, p3);
        std::swap(depths[0], depths[2]);
    } 
    if (v2.y>v3.y)
    { 
        std::swap(v2, v3);
        std::swap(t2, t3);
        std::swap(p2, p3);
        std::swap(depths[1], depths[2]);
    } 
    int total_height = v3.y-v1.y; 
    Vec2i vertex[3] = {v1, v2, v3};

    for (int i=0; i<total_height; i++) 
    { 
        bool second_half = i>v2.y-v1.y || v2.y==v1.y; 
        int segment_height = second_half ? v3.y-v2.y : v2.y-v1.y; 
        float alpha = (float)i/total_height; 
        float beta  = (float)(i-(second_half ? v2.y-v1.y : 0))/segment_height;  
        Vec2i A =               v1 + (v3-v1)*alpha; 
        Vec2i B = second_half ? v2 + (v3-v2)*beta : v1 + (v2-v1)*beta; 
        if (A.x>B.x)
        {
            std::swap(A, B);
        } 

        for (int j=A.x; j<=B.x; j++) 
        { 
            Vec3f bc_coord = getBarycenterCoord(vertex, Vec2i(j, v1.y+i));
            float z = bc_coord.raw[0]*p1.z + bc_coord.raw[1]*p2.z + bc_coord.raw[2]*p3.z;
            float uv_x = bc_coord.raw[0]*t1.x + bc_coord.raw[1]*t2.x + bc_coord.raw[2]*t3.x;
            float uv_y = bc_coord.raw[0]*t1.y + bc_coord.raw[1]*t2.y + bc_coord.raw[2]*t3.y;
            float depth = bc_coord.raw[0]*depths[0] + bc_coord.raw[1]*depths[1] + bc_coord.raw[2]*depths[2];
            this->set(j, v1.y+i, z, model->getTextureColor(uv_x, uv_y) * depth);
        } 
    } 

    return;
}

void HzzzImage::RenderModel(Model* model, TGAImage& texture)
{
    try
    {
        for(int i=0; i<model->nfaces(); i++)
        {
            std::vector<int> face = model->face(i);
            std::vector<int> t_face = model->texture_face(i);
            Vec3f world_coords[3];
            Vec2i screen_coords[3];
            Vec2f uv_coords[3];
            for(int j=0; j<3; j++)
            {
                Vec3f v = model->vert(face[j]);
                screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
                world_coords[j] = v;
                uv_coords[j] = model->uv(t_face[j]);
            }
            this->Triangle(world_coords[0], world_coords[1], world_coords[2], uv_coords[0], uv_coords[1], uv_coords[2], texture);
        }
        // std::cout << "Model rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

void HzzzImage::RenderModel(Model* model, TGAImage& texture, Vec3f light)
{
    try
    {
        for(int i=0; i<model->nfaces(); i++)
        {
            std::vector<int> face = model->face(i);
            std::vector<int> t_face = model->texture_face(i);
            Vec3f world_coords[3];
            Vec2i screen_coords[3];
            Vec2f uv_coords[3];
            Vec3f face_norm = model->face_norm(i);
            for(int j=0; j<3; j++)
            {
                Vec3f v = model->vert(face[j]);
                screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
                world_coords[j] = v;
                uv_coords[j] = model->uv(t_face[j]);
            }
            light.normalize();
            float depth = face_norm * light;
            depth = depth/2 + 0.5;
            this->Triangle(world_coords[0], world_coords[1], world_coords[2], 
                uv_coords[0], uv_coords[1], uv_coords[2], texture, depth);
        }
        // std::cout << "Model rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

void HzzzImage::RenderModel(Model* model, Vec3f light)
{
    try
    {
        for(int i=0; i<model->nfaces(); i++)
        {
            if(model->is_random_model)
            {
                this->Triangle(model->vert(0), model->vert(1), model->vert(2), 
                    Vec2f(0, 1), Vec2f(1, 0), Vec2f(1, 1), model, new float[3]{1, 1, 1});
                continue;
            }
            std::vector<int> face = model->face(i);
            std::vector<int> t_face = model->texture_face(i);
            Vec3f world_coords[3];
            Vec2i screen_coords[3];
            Vec2f uv_coords[3];
            // Vec3f face_norm = model->face_norm(i);
            float depths[3];
            for(int j=0; j<3; j++)
            {
                Vec3f v = model->vert(face[j]);
                screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
                world_coords[j] = v;
                uv_coords[j] = model->uv(t_face[j]);
                depths[j] = model->norm(face[j]) * light;
                // if(0.9 < depths[j] && depths[j] < 1.0)
                // {
                //     depths[j] = 2.0;
                // }
                depths[j] = depths[j]/2.0 + 0.5;
            }
            light.normalize();
            
            this->Triangle(world_coords[0], world_coords[1], world_coords[2], 
                uv_coords[0], uv_coords[1], uv_coords[2], model, depths);
        }
        // std::cout << "Model rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

void HzzzImage::RayTracingRender(Model* model, TGAImage& texture, float fov, std::vector<Model*> lights)
{
    try
    {
        // int max_threads = MAX_MECHINE_THREADS;
        // int pix_stp = this->width/max_threads;
        // TGAImage* texture_ptr = &texture;
        // for(int threads_num=0; threads_num<max_threads; threads_num++)
        // {
        //     std::thread t
        //     (
        //         [=]()
        //         {
        //             RayTraceable ray_traceable(model);
        //             for(int pix_x=threads_num*pix_stp; pix_x<(threads_num+1)*pix_stp; pix_x++)
        //             {
        //                 for(int pix_y=0; pix_y<this->height; pix_y++)
        //                 {
        //                     float x = screen_width *pix_x/width - 0.5*screen_width ;
        //                     float y = screen_height*pix_y/width - 0.5*screen_height;
        //                     Vec3f dir = Vec3f(x, y, 1);
        //                     dir.normalize();
        //                     Ray ray = Ray(Vec3f(0, 0, -3), dir);
        //                     TGAColor res_color;
        //                     float res_lightness;
        //                     ray_traceable.Intersect(ray, res_color, res_lightness, 1.0, lights);
        //                     // std::cout << res_lightness << " ";
        //                     this->set(pix_x, pix_y, res_color*res_lightness);
        //                 }
        //                 // std::cout << "One row rendered: " << pix_x << std::endl;
        //             }
        //         }
        //     );
        //     if(threads_num == max_threads - 1)
        //     {
        //         t.join();
        //     }
        //     else 
        //     {
        //         t.detach();
        //     }
            
        // }

        // for(int pix_x=0; pix_x<this->width; pix_x++)
        // // for(int pix_x=300; pix_x<700; pix_x++)
        // {
        //     for(int pix_y=0; pix_y<this->height; pix_y++)
        //     // for(int pix_y=300; pix_y<700; pix_y++)
        //     {
        //         // float x = (2 * (pix_x + 0.5) / (float)width - 1) * tan(fov / 2.0) * width / (float)height;
        //         // float y = -(2 * (pix_y + 0.5) / (float)height - 1) * tan(fov / 2.0);
        //         // float z = -1;
        //         float x = screen_width *pix_x/width - 0.5*screen_width ;
        //         float y = screen_height*pix_y/width - 0.5*screen_height;
        //         Vec3f dir = Vec3f(x, y, 1);
        //         dir.normalize();
        //         Ray ray = Ray(Vec3f(0, 0, -3), dir);
        //         float min_t = std::numeric_limits<float>::max();
        //         Vec3f intersectPoint;
        //         int intersectFaceIdx = -1;
        //         Vec3f bc_coord;
        //         Vec3f render_bc_coord;
        //         for(int i=0; i<model->nfaces(); i++)
        //         {
        //             std::vector<int> face = model->face(i);
        //             Vec3f points[3];
        //             for(int j=0; j<3; j++)
        //             {
        //                 points[j] = model->vert(face[j]);
        //             }
        //             if(ray.isIntersect(points, intersectPoint, bc_coord))
        //             {
        //                 float t = (intersectPoint - ray.origin).norm();
        //                 if(t < min_t)
        //                 {
        //                     min_t = t;
        //                     intersectFaceIdx = i;
        //                     render_bc_coord = getBarycenterCoord(points, intersectPoint);
        //                     // std::cout << "Intersect point: " << intersectPoint << std::endl;
        //                 }
        //             }
        //         }
        //         if(intersectFaceIdx != -1)
        //         {
        //             std::vector<int> face = model->face(intersectFaceIdx);
        //             std::vector<int> t_face = model->texture_face(intersectFaceIdx);
        //             Vec3f world_coords[3];
        //             // Vec2i screen_coords[3];
        //             Vec2f uv_coords[3];
        //             // Vec3f face_norm = model->face_norm(intersectFaceIdx);
        //             for(int j=0; j<3; j++)
        //             {
        //                 Vec3f v = model->vert(face[j]);
        //                 // screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);
        //                 world_coords[j] = v;
        //                 uv_coords[j] = model->uv(t_face[j]);
        //             }
                    
        //             // if(render_bc_coord.raw[0] > 1 || render_bc_coord.raw[0] < 0 || render_bc_coord.raw[1] > 1 || render_bc_coord.raw[1] < 0 || render_bc_coord.raw[2] > 1 || render_bc_coord.raw[2] < 0)
        //             // {
        //             //     continue;
        //             // }
        //             // std::cout << bc_coord.raw[0] << " " << bc_coord.raw[1] << " " << bc_coord.raw[2] << "#";
        //             // std::cout<< render_bc_coord.raw[0] << " " << render_bc_coord.raw[1] << " " << render_bc_coord.raw[2] << std::endl;
        //             // float color_coord_x = uv_coords[0].x*bc_coord.raw[0] + uv_coords[1].x*bc_coord.raw[1] + uv_coords[2].x*bc_coord.raw[2];
        //             // float color_coord_y = uv_coords[0].y*bc_coord.raw[0] + uv_coords[1].y*bc_coord.raw[1] + uv_coords[2].y*bc_coord.raw[2];

        //             float color_coord_x = uv_coords[0].x*render_bc_coord.raw[0] + uv_coords[1].x*render_bc_coord.raw[1] + uv_coords[2].x*render_bc_coord.raw[2];
        //             float color_coord_y = uv_coords[0].y*render_bc_coord.raw[0] + uv_coords[1].y*render_bc_coord.raw[1] + uv_coords[2].y*render_bc_coord.raw[2];
        //             this->set(pix_x, pix_y, intersectPoint.z, texture.get(color_coord_x, color_coord_y));
        //         }
        //     }
        //     std::cout << "One row rendered: " << pix_x << std::endl;
        // }
        // std::cout << "Model rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

void HzzzImage::RayTracingRender(RayTraceable* ray_traceable, float fov, std::vector<Model*> lights)
{
    try
    {
        int max_threads = MAX_MECHINE_THREADS;
        int pix_stp = this->width/max_threads;
        int* finish_threads = new int[1];
        *finish_threads = 0;

        for(int threads_num=0; threads_num<max_threads; threads_num++)
        {
            auto shared_operate = 
            [=]()
            {
                for(int pix_x=threads_num*pix_stp; pix_x<(threads_num+1)*pix_stp; pix_x++)
                {
                    for(int pix_y=0; pix_y<this->height; pix_y++)
                    {
                        float x = screen_width *pix_x/width - 0.5*screen_width ;
                        float y = screen_height*pix_y/width - 0.5*screen_height;
                        Vec3f dir = Vec3f(x, y, 1);
                        dir.normalize();
                        Ray ray = Ray(this->getCamera(), dir);
                        TGAColor res_color;
                        float res_lightness;
                        Vec3f intersectPoint;
                        if(ray_traceable->Intersect(ray, res_color, res_lightness, intersectPoint, 1.0, lights))
                        {
                            this->set(pix_x, pix_y, res_color*res_lightness);
                        }                        
                    }
                    std::cout << "One row rendered: " << pix_x << std::endl;
                }
            };

            if(threads_num == max_threads - 1)
            {
                std::thread t
                (
                    [=]()
                    {
                        // for(int pix_x=threads_num*pix_stp; pix_x<(threads_num+1)*pix_stp; pix_x++)
                        // {
                        //     for(int pix_y=0; pix_y<this->height; pix_y++)
                        //     {
                        //         float x = screen_width *pix_x/width - 0.5*screen_width ;
                        //         float y = screen_height*pix_y/width - 0.5*screen_height;
                        //         Vec3f dir = Vec3f(x, y, 1);
                        //         dir.normalize();
                        //         Ray ray = Ray(this->getCamera(), dir);
                        //         TGAColor res_color;
                        //         float res_lightness;
                        //         Vec3f intersectPoint;
                        //         ray_traceable->Intersect(ray, res_color, res_lightness, intersectPoint, 1.0, lights);
                        //         this->set(pix_x, pix_y, intersectPoint.z, res_color*res_lightness);
                        //         // this->set(pix_x, pix_y, res_color);
                        //     }
                        //     std::cout << "One row rendered: " << pix_x << std::endl;
                        // }
                        shared_operate();
                        finish_threads[0]++;
                        while(finish_threads[0] < max_threads)
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(100));
                        }
                    }
                );
                t.join();
            }
            else 
            {
                std::thread t
                (
                    [=]()
                    {
                        // for(int pix_x=threads_num*pix_stp; pix_x<(threads_num+1)*pix_stp; pix_x++)
                        // {
                        //     for(int pix_y=0; pix_y<this->height; pix_y++)
                        //     {
                        //         float x = screen_width *pix_x/width - 0.5*screen_width ;
                        //         float y = screen_height*pix_y/width - 0.5*screen_height;
                        //         Vec3f dir = Vec3f(x, y, 1);
                        //         dir.normalize();
                        //         Ray ray = Ray(this->getCamera(), dir);
                        //         TGAColor res_color;
                        //         float res_lightness;
                        //         Vec3f intersectPoint;
                        //         ray_traceable->Intersect(ray, res_color, res_lightness, intersectPoint, 1.0, lights);
                        //         this->set(pix_x, pix_y, intersectPoint.z, res_color*res_lightness);
                        //         // this->set(pix_x, pix_y, res_color);
                        //     }
                        //     std::cout << "One row rendered: " << pix_x << std::endl;
                        // }
                        shared_operate();
                        finish_threads[0]++;
                        while(finish_threads[0] < max_threads)
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(100));
                        }
                        finish_threads[0]++;
                    }
                );
                t.detach();
            }
            
        }

        // std::cout << "Traceable rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }

}

class Allocer 
{
public:
    int pix_x = 0;
    int max_x;
    int pix_y = 0;
    int max_y;
    bool is_finished = false;
    std::mutex mt;

    Allocer(int _max_x, int _max_y)
    {
        this->max_x = _max_x;
        this->max_y = _max_y;
    }

    bool isFinished(int& pix_x, int& pix_y)
    {
        mt.lock();
        bool is_finished = this->is_finished;
        if(is_finished)
        {
            mt.unlock();
            return is_finished;
        }
        pix_x = this->pix_x;
        pix_y = this->pix_y;
        this->pix_x += 1;
        if(this->pix_x == this->max_x)
        {
            this->pix_x = 0;
            this->pix_y += 1;
            if(this->pix_y == this->max_y)
            {
                this->is_finished = true;
            }
        }
        mt.unlock();
        return is_finished;
    }
};

void HzzzImage::RayTracingRender(std::vector<RayTraceable*> ray_traceables, std::vector<Model*> lights)
{
    try
    {
        int max_threads = MAX_MECHINE_THREADS;
        int pix_stp = this->width/max_threads;
        int* finish_threads = new int[1];
        *finish_threads = 0;
        Allocer* ac = new Allocer(this->width, this->height);

        for(int threads_num=0; threads_num<max_threads; threads_num++)
        {
            // auto shared_operate = 
            // [=]()
            // {
            //     for(int pix_x=threads_num*pix_stp; pix_x<(threads_num+1)*pix_stp; pix_x++)
            //     {
            //         for(int pix_y=0; pix_y<this->height; pix_y++)
            //         {
            //             // Vec3f rotate_vec = this->getViewDirection();
            //             // float theta = -std::atan2(rotate_vec.x, rotate_vec.y);
            //             // float phi = std::atan2(std::sqrt(rotate_vec.x*rotate_vec.x+rotate_vec.y*rotate_vec.y), rotate_vec.z);
            //             float x = screen_width *pix_x/width - 0.5*screen_width ;
            //             float y = screen_height*pix_y/width - 0.5*screen_height;
            //             Vec3f dir = Vec3f(x, y, std::tan(fov/2.0)*screen_height);
            //             dir.normalize();
            //             Ray ray = Ray(this->getCamera(), dir);
            //             TGAColor res_color;
            //             float res_lightness;
            //             Vec3f intersectPoint;
            //             float distance = std::numeric_limits<float>::max();
            //             for(RayTraceable* ray_traceable: ray_traceables)
            //             {
            //                 float temp_res_lightness;
            //                 TGAColor temp_res_color;
            //                 if(ray_traceable->Intersect(ray, temp_res_color, temp_res_lightness, intersectPoint, 1.0, lights))
            //                 {
            //                     float t = (intersectPoint - ray.origin).norm();
            //                     if(t < distance)
            //                     {
            //                         distance = t;
            //                         res_lightness = temp_res_lightness;
            //                         res_color = temp_res_color;
            //                     }
            //                 }
            //             }
            //             this->set(pix_x, pix_y, res_color*res_lightness);
            //         }
                    
            //         #ifndef SILENCE
            //         std::cout << "One row rendered: " << pix_x << std::endl;
            //         // std::cout << "#" << pix_x;
            //         #endif
            //     }
            // };
            
            auto shared_operate = 
            [=]()
            {
                int pix_x;
                int pix_y;
                while(!ac->isFinished(pix_x, pix_y))
                {
                    float x = screen_width *pix_x/width - 0.5*screen_width ;
                    float y = screen_height*pix_y/width - 0.5*screen_height;
                    Vec3f mid_axis_dir = this->getViewDirection().normalize();
                    Vec3f pix_x_dir = Vec3f(mid_axis_dir.z, 0, -mid_axis_dir.x);
                    pix_x_dir.normalize();
                    Vec3f pix_y_dir = Vec3f(0, 0, 0) - pix_x_dir^mid_axis_dir;
                    pix_y_dir.normalize();
                    Vec3f dir = mid_axis_dir*(std::tan(fov/2.0)*screen_height) 
                                + pix_x_dir*x + pix_y_dir*y;
                    // Vec3f pix_point = Vec3f(x, y, std::tan(fov/2.0)*screen_height);
                    // Vec3f dir = Vec3f(x, y, std::tan(fov/2.0)*screen_height);
                    dir.normalize();
                    Ray ray = Ray(this->getCamera(), dir);
                    TGAColor res_color;
                    float res_lightness;
                    Vec3f intersectPoint;
                    float distance = std::numeric_limits<float>::max();
                    for(RayTraceable* ray_traceable: ray_traceables)
                    {
                        float temp_res_lightness;
                        TGAColor temp_res_color;
                        if(ray_traceable->Intersect(ray, temp_res_color, temp_res_lightness, intersectPoint, 1.0, lights))
                        {
                            float t = (intersectPoint - ray.origin).norm();
                            if(t < distance)
                            {
                                distance = t;
                                res_lightness = temp_res_lightness;
                                res_color = temp_res_color;
                            }
                        }
                    }
                    this->set(pix_x, pix_y, res_color*res_lightness);
                }
            };

            if(threads_num == max_threads - 1)
            {
                std::thread t
                (
                    [=]()
                    {
                        shared_operate();
                        finish_threads[0]++;
                        while(finish_threads[0] < max_threads)
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(100));
                        }
                    }
                );
                t.join();
            }
            else 
            {
                std::thread t
                (
                    [=]()
                    {
                        shared_operate();
                        finish_threads[0]++;
                    }
                );
                t.detach();
            }
            
        }

        std::cout << "Traceable rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

void HzzzImage::RayTracingRender(RayTracer* rt)
{
    try
    {
        int max_threads = MAX_MECHINE_THREADS;
        int pix_stp = this->width/max_threads;
        int* finish_threads = new int[1];
        *finish_threads = 0;
        Allocer* ac = new Allocer(this->width, this->height);

        for(int threads_num=0; threads_num<max_threads; threads_num++)
        {
            auto shared_operate = 
            [=]()
            {
                int pix_x;
                int pix_y;
                while(!ac->isFinished(pix_x, pix_y))
                {
                    float x = screen_width *pix_x/width - 0.5*screen_width ;
                    float y = screen_height*pix_y/width - 0.5*screen_height;
                    Vec3f mid_axis_dir = this->getViewDirection().normalize();
                    Vec3f pix_x_dir = Vec3f(mid_axis_dir.z, 0, -mid_axis_dir.x);
                    pix_x_dir.normalize();
                    Vec3f pix_y_dir = Vec3f(0, 0, 0) - pix_x_dir^mid_axis_dir;
                    pix_y_dir.normalize();
                    Vec3f dir = mid_axis_dir*(std::tan(fov/2.0)*screen_height) 
                                + pix_x_dir*x + pix_y_dir*y;
                    // Vec3f pix_point = Vec3f(x, y, std::tan(fov/2.0)*screen_height);
                    // Vec3f dir = Vec3f(x, y, std::tan(fov/2.0)*screen_height);
                    dir.normalize();
                    Ray ray = Ray(this->getCamera(), dir);
                    TGAColor res_color;
                    float res_lightness;
                    Vec3f intersectPoint;
                    float distance = std::numeric_limits<float>::max();
                    for(RayTraceable* ray_traceable: rt->scene)
                    {
                        float temp_res_lightness;
                        TGAColor temp_res_color;
                        if(ray_traceable->Intersect(ray, temp_res_color, temp_res_lightness, intersectPoint, 1.0, rt->lights))
                        {
                            float t = (intersectPoint - ray.origin).norm();
                            if(t < distance)
                            {
                                distance = t;
                                res_lightness = temp_res_lightness;
                                res_color = temp_res_color;
                            }
                        }
                    }
                    this->set(pix_x, pix_y, res_color*res_lightness);
                }
            };

            if(threads_num == max_threads - 1)
            {
                std::thread t
                (
                    [=]()
                    {
                        shared_operate();
                        finish_threads[0]++;
                        while(finish_threads[0] < max_threads)
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(100));
                        }
                    }
                );
                t.join();
            }
            else 
            {
                std::thread t
                (
                    [=]()
                    {
                        shared_operate();
                        finish_threads[0]++;
                    }
                );
                t.detach();
            }
            
        }

        std::cout << "Traceable rendered." << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

TGAImage & HzzzImage::toTGAImage(TGAImage &image)
{
    image = TGAImage(this->width, this->height, TGAImage::RGB);
    for (int i = 0; i < this->width; i++)
    {
        for (int j = 0; j < this->height; j++)
        {
            int idx = i + j * this->width;
            image.set(i, j, TGAColor(this->image[idx*4], this->image[idx*4 + 1], this->image[idx*4 + 2]));
        }
    }
    return image;
}

TGAImage HzzzImage::toTGAImage()
{
    TGAImage image = TGAImage(this->width, this->height, TGAImage::RGB);
    for (int i = 0; i < this->width; i++)
    {
        for (int j = 0; j < this->height; j++)
        {
            int idx = i + j * this->width;
            image.set(i, j, TGAColor(this->image[idx*4], this->image[idx*4 + 1], this->image[idx*4 + 2]));
        }
    }
    return image;
}

Vec2i HzzzImage::toImageCoord(Vec3f v)
{
    Vec2i imageCoord;
    imageCoord.x = (v.x+1.)*width/2.;
    imageCoord.y = (v.y+1.)*height/2.;
    return imageCoord;
}

void HzzzImage::toFile(std::string filename)
{
    std::fstream fout;
    fout.open(filename, std::ios::out);

    fout << this->width << " " << this->height << std::endl;
    for(int i=0; i<this->width; i++)
    {
        for(int j=0; j<this->height; j++)
        {
            int idx = i + j * this->width;
            fout << int(static_cast<unsigned char>(this->image[idx*4])) << " " 
                 << int(static_cast<unsigned char>(this->image[idx*4 + 1])) << " " 
                 << int(static_cast<unsigned char>(this->image[idx*4 + 2])) << " ";
        }
        fout << std::endl;
    }
}

void HzzzImage::saveToTGAImage(std::string filename)
{
    TGAImage image = this->toTGAImage();
    image.flip_vertically();
    image.write_tga_file(filename.c_str());
}

HzzzImage::~HzzzImage()
{}



void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
    bool steep = false;
    if (std::abs(x0-x1)<std::abs(y0-y1)) {
        std::swap(x0, y0);
        std::swap(x1, y1);
        steep = true;
    }
    if (x0>x1) {
        std::swap(x0, x1);
        std::swap(y0, y1);
    }

    for (int x=x0; x<=x1; x++) {
        float t = (x-x0)/(float)(x1-x0);
        int y = y0*(1.-t) + y1*t;
        if (steep) {
            image.set(y, x, color);
        } else {
            image.set(x, y, color);
        }
    }
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) 
{ 
    if (t0.y==t1.y && t0.y==t2.y) return; // I dont care about degenerate triangles 
    // sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
    if (t0.y>t1.y) std::swap(t0, t1); 
    if (t0.y>t2.y) std::swap(t0, t2); 
    if (t1.y>t2.y) std::swap(t1, t2); 
    int total_height = t2.y-t0.y; 
    for (int i=0; i<total_height; i++) 
    { 
        bool second_half = i>t1.y-t0.y || t1.y==t0.y; 
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y; 
        float alpha = (float)i/total_height; 
        float beta  = (float)(i-(second_half ? t1.y-t0.y : 0))/segment_height; // be careful: with above conditions no division by zero here 
        Vec2i A =               t0 + (t2-t0)*alpha; 
        Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta; 
        if (A.x>B.x) std::swap(A, B); 
        for (int j=A.x; j<=B.x; j++) 
        { 
            image.set(j, t0.y+i, color); // attention, due to int casts t0.y+i != A.y 
        } 
    } 
}

void mytriangle(Vec3f v1, Vec3f v2, Vec3f v3, int height, int width, TGAImage &image, TGAColor color)
{
    Vec2i t0, t1, t2;
    t0.x = (v1.x+1.)*width/2.;
    t0.y = (v1.y+1.)*height/2.;
    t1.x = (v2.x+1.)*width/2.;
    t1.y = (v2.y+1.)*height/2.;
    t2.x = (v3.x+1.)*width/2.;
    t2.y = (v3.y+1.)*height/2.;
    triangle(t0, t1, t2, image, color);

    return;
}

void rotateVec(Vec3f& v, float theta, float phi)
{
    v.x = v.x*cos(theta) + v.y*sin(theta);
    v.y = -v.x*sin(theta) + v.y*cos(theta);
    v.z = v.z*cos(phi) + v.y*sin(phi);

    return;
}

