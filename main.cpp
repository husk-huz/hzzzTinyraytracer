#include <vector>
#include <cmath>
#include <iostream> 
#include <string>
#include "include/tgaimage.h"
#include "include/model.h"
#include "include/geometry.h"
#include "include/HzzzImage.h"
#include "include/RayTraceable.h"
#include "include/RayTracer.h"
#include "include/config.h"

#define DEBUG

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;

#ifdef HIGH_RES
    const int width  = 1000;
    const int height = 1000;
#else 
    const int width  = 100;
    const int height = 100;
#endif



// float zbuffer[width][height];

HzzzImage RenderImage(Model* model, TGAImage& texture)
{
    try
    {
        HzzzImage image(width, height);
        for(int i=0; i<model->nfaces(); i++)
        // for(int i=0; i<; i++)
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
            image.Triangle(world_coords[0], world_coords[1], world_coords[2], uv_coords[0], uv_coords[1], uv_coords[2], texture);
            // image.Line(screen_coords[0], screen_coords[1], white);
            // image.Line(screen_coords[1], screen_coords[2], white);
            // image.Line(screen_coords[2], screen_coords[0], white);
        }

        return image;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "Fail to render." << std::endl;
    }
}

RayTracer GenCornellBox()
{
    RayTracer tracer;
    tracer.scene = std::vector<RayTraceable*>();
    RayTraceable* object = new RayTraceableInfPlane(0, 0, 1, -3);
    object->color = white;
    tracer.addSon(object);
    object = new RayTraceableInfPlane(1, 0, 0, 1.5);
    object->color = red;
    tracer.addSon(object);
    object = new RayTraceableInfPlane(1, 0, 0, -1.5);
    object->color = green;
    tracer.addSon(object);
    object = new RayTraceableInfPlane(0, 1, 0, 3);
    object->color = white;
    tracer.addSon(object);
    object = new RayTraceableInfPlane(0, 1, 0, 0);
    object->color = white;
    tracer.addSon(object);
    object = new RayTraceableInfPlane(0, 1, 0, -3);
    object->color = white;
    tracer.addSon(object);

    float a = 0.5;
    object = new RayTraceableBox(std::vector<Vec3f>{Vec3f(-a, a, a), Vec3f(a, a, a), 
                                                    Vec3f(-a, a, -a), Vec3f(a, a, -a), 
                                                    Vec3f(-a, -a, -a), Vec3f(a, -a, -a), 
                                                    Vec3f(-a, -a, a), Vec3f(a, -a, a)});
    object->color = white;
    static_cast<RayTraceableBox*>(object)->rotate(0.2, Vec3f(0, 1, 0));
    static_cast<RayTraceableBox*>(object)->move(Vec3f(0.3, a, 0));
    tracer.addSon(object);

    object = new RayTraceableBox(std::vector<Vec3f>{Vec3f(-a, 2*a, a), Vec3f(a, 2*a, a), 
                                                    Vec3f(-a, 2*a, -a), Vec3f(a, 2*a, -a), 
                                                    Vec3f(-a, -2*a, -a), Vec3f(a, -2*a, -a), 
                                                    Vec3f(-a, -2*a, a), Vec3f(a, -2*a, a)});
    object->color = white;
    static_cast<RayTraceableBox*>(object)->rotate(-0.2, Vec3f(0, 1, 0));
    static_cast<RayTraceableBox*>(object)->move(Vec3f(-0.3, 2*a, 2));
    tracer.addSon(object);

    return tracer;
}

int main()
{
    RayTracer rt = GenCornellBox();
    HzzzImage image(width, height);
    image.setCamera(Vec3f(0, 1, -10));
    image.setViewDirection(Vec3f(0, 0, 1));
    std::vector<Model*> light_models;
    for(int i=0; i<4; i++)
    {
        light_models.push_back(new Model());
    }
    image.RayTracingRender(rt.scene, light_models);
    image.saveToTGAImage("output.tga");

    return 0;
}

int not_main(int argc, char** argv) 
{
    // TGAImage frame(200, 200, TGAImage::RGB); 
    // Vec2i pts[3] = {Vec2i(10,10), Vec2i(100, 30), Vec2i(190, 160)}; 
    // triangle(pts[0], pts[1], pts[2], frame, TGAColor(255, 0, 0)); 
    // frame.flip_vertically(); // to place the origin in the bottom left corner of the image 
    // frame.write_tga_file("framebuffer.tga");
    // return 0; 

    // for(int i=0; i<width; i++)
    // {
    //     for(int j=0; j<height; j++)
    //     {
    //         zbuffer[i][j] = -std::numeric_limits<float>::max();
    //     }
    // }
    // Model* model_rock = new Model("obj/models/rock/rock.obj");
    // if (2==argc) 
    // {
    //     model = new Model(argv[1]);
    // } else 
    // {
    //     // model = new Model("/obj/spot_triangulated.obj");
    //     model = new Model("obj/new.obj");
    //     // model->calNorm();
    //     // std::cout << model->toString() << std::endl;
    //     // model = new Model("obj/models/rock/rock.obj");
        
    //     // model = new Model("obj/african_head.obj");
    //     std::cout << model->toString() << std::endl;
    //     std::cout << "model loaded" << std::endl;
    // }

    model = new Model("obj/new.obj");
    TGAImage texture;
    // TGAImage texture_rock;
    // texture.read_tga_file("obj/african_head_diffuse.tga"); 
    // texture.read_tga_file("obj/models/rock/rock.tga"); 
    texture.read_tga_file("obj/spot_texture.tga"); 
    texture.flip_vertically();
    model->setTexture(&texture);
    // texture_rock.read_tga_file("obj/models/rock/rock.tga");
    // texture.flip_vertically();
    // texture_rock.flip_vertically();
    // model->setTexture(&texture);
    // model_rock->setTexture(&texture_rock);

    float d = M_PI/30;
    std::vector<Vec3f> light_vector;
    std::vector<Model*> light_models;
    
    // model_rock->scale(0.1);
    // RayTraceable rt(Vec3f(0, 1, 0), 1);
    // RayTraceable sphere(Vec3f(1, 0.2, 0), 0.2);
    // RayTraceable plane(0, 1, 0, 0);
    // RayTraceable plane2(1, 1, 1, -3);
    RayTraceable** rts = new RayTraceable*[10];
    RayTracer tracer;
    rts[0] = new RayTraceableInfPlane(0, 1, 0, 1);
    rts[0]->color = TGAColor(200, 120, 200, 255);
    rts[0]->father_tracer = &tracer;
    rts[1] = new RayTraceableInfPlane(0, 0, 1, -1.5);
    rts[1]->color = TGAColor(200, 200, 120, 255);
    rts[1]->father_tracer = &tracer;
    rts[2] = new RayTraceableSphere(Vec3f(1, 0.2, 0), 0.2);
    rts[2]->color = TGAColor(120, 200, 200, 255);
    rts[2]->father_tracer = &tracer;
    rts[3] = new RayTraceableSphere(Vec3f(-2, 1, 0), 1);
    rts[3]->color = TGAColor(255, 255, 255, 255);
    rts[3]->father_tracer = &tracer;
    // rts[4] = new RayTraceableModel(model);
    float a = 0.3;
    rts[4] = new RayTraceableBox(std::vector<Vec3f>{Vec3f(-a, a, a), Vec3f(a, a, a), 
                                                    Vec3f(-a, a, -a), Vec3f(a, a, -a), 
                                                    Vec3f(-a, -a, -a), Vec3f(a, -a, -a), 
                                                    Vec3f(-a, -a, a), Vec3f(a, -a, a)});
    static_cast<RayTraceableBox*>(rts[4])->rotate(0.5, Vec3f(0, 1, 0));
    static_cast<RayTraceableBox*>(rts[4])->rotate(0.5, Vec3f(1, 0, 0));
    static_cast<RayTraceableBox*>(rts[4])->move(Vec3f(0, -0.3, 0));
    rts[4]->color = TGAColor(255, 255, 255, 255);
    rts[4]->father_tracer = &tracer;

    tracer.scene = std::vector<RayTraceable*>();
    // tracer.lights = std::vector<Model*>();

    for(int i=0; i<5; i++)
    {
        tracer.scene.push_back(rts[i]);
    }

    HzzzImage hi(width, height);
    hi.setCamera(Vec3f(0, 1, -10));
    hi.setViewDirection(Vec3f(0, 0, 1));
    
    for(int i=0; i<4; i++)
    {
        light_models.push_back(new Model());
    }

    // hi.RayTracingRender(tracer.scene, light_models);

    // HzzzImage hi(width, height);
    // hi.setCamera(Vec3f(0, 1, -10));
    // hi.setViewDirection(Vec3f(0, 0, 1));
    // // std::vector<RayTraceable*> ray_traceables;
    // // ray_traceables = std::vector<RayTraceable*>();
    // rt.color = TGAColor(255, 255, 255, 255);
    // plane.color = TGAColor(200, 120, 200, 255);
    // plane2.color = TGAColor(200, 200, 120, 255);
    // sphere.color = TGAColor(120, 200, 200, 255);
    // tracer.scene = std::vector<RayTraceable*>();
    // tracer.scene.push_back(&rt);
    // rt.father_tracer = &tracer;
    // tracer.scene.push_back(&sphere);
    // sphere.father_tracer = &tracer;
    // tracer.scene.push_back(&plane);
    // plane.father_tracer = &tracer;
    // tracer.scene.push_back(&plane2);
    // plane2.father_tracer = &tracer;
    // // ray_traceables.push_back(&rt);
    // // ray_traceables.push_back(&plane);
    // // ray_traceables.push_back(&plane2);
    // // ray_traceables.push_back(&sphere);
    // hi.RayTracingRender(tracer.scene, light_models);
    for(int i=0; i<1; i++)
    {
        HzzzImage hz(width, height);
        hz.setCamera(Vec3f(0, 1, -10));
        hz.setViewDirection(Vec3f(-0.1, -0.1, 1));
        hz.RayTracingRender(tracer.scene, light_models);
        // hz.toFile("outfiles/" + std::to_string(i) + std::string(".hzzz"));
        TGAImage tga_image = hz.toTGAImage();
        tga_image.flip_vertically();
        tga_image.write_tga_file("output.tga");
        static_cast<RayTraceableBox*>(rts[4])->rotate(d, Vec3f(0, 1, 0));
        for(int j=0; j<4; j++)
        {
            light_models[j]->move(Vec3f(cos(d*i)/4, 0, 0));
        }
        std::cout << "#" << i << std::endl;
    }
    // hi.RayTracingRender(&rt, 90, light_models);
    // hi.RayTracingRender(&plane, 90, light_models);
    // TGAImage tga_image = hi.toTGAImage();
    // tga_image.flip_vertically();
    // tga_image.write_tga_file("output.tga");
    std::cout << "C++ ALL DONE" << std::endl;

    // image.RayTracingRender(&rt, 90, light_models);
    // for(int i=0; i<60; i++)
    // {
    //     // image.RayTracingRender(model, texture, 90, light_models);
    //     HzzzImage image(width, height);
    //     // image.RayTracingRender(&rt, 90, light_models);
    //     image.setCamera(Vec3f(0, 1, -4));
    //     image.setViewDirection(Vec3f(0, 0, 1));
    //     image.RayTracingRender(ray_traceables, light_models);
    //     image.toFile("outfiles/" + std::to_string(i) + std::string(".hzzz"));
    //     std::cout << "#" ;
    //     for(int j=0; j<5; j++)
    //     {
    //         light_models[j]->move(Vec3f(-0.05, 0, 0));
    //     }
    //     // for(int j=0; j<5; j++)
    //     // {
    //     //     light_models[j]->scale(1.05);
    //     // }
    // }
    // image.RayTracingRender(model, texture, 90, light_models);
    // image.RenderModel(model, texture, Vec3f(cos(d*2), sin(d*2), -1));
    // image.RenderModel(model_rock, texture_rock, Vec3f(cos(d*2), sin(d*2), -1));
    // image.RenderModel(model, Vec3f(cos(d*2), sin(d*2), -1));
    // for(int i=0; i<100; i++)
    // {
    //     image.RenderModel(new Model(), Vec3f(cos(d*i), sin(d*i), -1));
    // }
    // image.RayTracingRender(model, texture, 90, light_vector);
    // TGAImage tga_image = image.toTGAImage();
    // tga_image.flip_vertically();
    // tga_image.write_tga_file("output.tga");

    // model->rotate(20*d);
    
    // int video = 1;
    // if(video)
    // {
    //     for(int k=0; k<60; k++)
    //     {
    //         // HzzzImage image = RenderImage(model, texture);
    //         // model->rotate(d*k);
    //         HzzzImage image(width, height);
    //         image.RenderModel(model, texture, Vec3f(cos(d*k*2), sin(d*k*2), -1));
    //         image.toFile("outfiles/" + std::to_string(k) + std::string(".hzzz"));
    //         // TGAImage tgaimage = image.toTGAImage();
    //         // tgaimage.flip_vertically();
    //         // tgaimage.write_tga_file(("outputs/" + std::to_string(k) + std::string(".tga")).data());
    //         model->rotate(d);
    //         std::cout << "#" ;
    //     }
    // }
    // else 
    // {
    //     int k = 2;
    //     HzzzImage image(width, height);
    //     image.RenderModel(model, texture, Vec3f(cos(d*k*2), sin(d*k*2), -1));
    //     TGAImage tgaimage = image.toTGAImage();
    //     tgaimage.flip_vertically();
    //     // tgaimage.write_tga_file(("outputs/" + std::to_string(k) + std::string(".tga")).data());
    //     tgaimage.write_tga_file("output.tga");
    // }
    
    

    
    // TGAImage tgaimage = image.toTGAImage();
    // tgaimage.flip_vertically();
    // tgaimage.write_tga_file("output.tga");


    // TGAImage image(width, height, TGAImage::RGB);
    // TGAImage uvmap;
    // uvmap.read_tga_file("obj/spot_texture.tga");
    
    // uvmap.read_tga_file("obj/african_head_diffuse.tga");
    // uvmap.read_tga_file("obj/models/rock/rock.tga");

    // std::cout << "read tga file" << std::endl;
    // uvmap.flip_vertically();
    // float l = 1.0;

    // Vec3f light_dir(0,0,-1); // define light_dir

    // std::cout << "nfaces: " << model->nfaces() << std::endl;
    // std::cout << "rendering" << std::endl;
    // int k=0;
    // for (int i=0; i<model->nfaces(); i++) 
    // {
    //     std::cout << "#" ;
    //     std::vector<int> face = model->face(i);
    //     std::vector<int> full_face = model->full_face(i);
    //     Vec2i screen_coords[3]; 
    //     Vec3f world_coords[3]; 
    //     Vec2f uv_coords[3];
    //     for (int j=0; j<3; j++) 
    //     { 
    //         // std::cout << face[j] << " " ;
    //         Vec3f v = model->vert(face[j]); 
    //         screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.); 
    //         world_coords[j]  = v; 
    //         uv_coords[j] = model->uv(full_face[3*j+1]); 
    //         // std::cout << full_face[j] << std::endl;
    //         // std::cout << uv_coords[j].x << " " << uv_coords[j].y << std::endl;
    //     } 
    //     // Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]); 
    //     // n.normalize(); 
    //     // float intensity = n*light_dir; 
    //     // if (intensity>0) 
    //     if(1)
    //     { 
    //         // triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255)); 
    //         mapTriangle(full_face, world_coords, uv_coords, height, width, image, uvmap);
            
    //         // k++;
    //         // if(k == 100)
    //         // {
    //         //     break;
    //         // }
    //     } 

    //     // mytriangle(model->vert(face[0]), model->vert(face[1]), model->vert(face[2]), height, width, image, TGAColor(rand()%255, rand()%255, rand()%255, 255));
    //     // for (int j=0; j<3; j++) {
    //     //     Vec3f v0 = model->vert(face[j]);
    //     //     Vec3f v1 = model->vert(face[(j+1)%3]);
    //     //     int x0 = (v0.x+1.)*width/2.;
    //     //     int y0 = (v0.y+1.)*height/2.;
    //     //     int x1 = (v1.x+1.)*width/2.;
    //     //     int y1 = (v1.y+1.)*height/2.;
    //     //     // int x0 = (v0.x*(l/(l+std::abs(v0.z+1.)))+1.)*width/2.;
    //     //     // int y0 = (v0.y*(l/(l+std::abs(v0.z+1.)))+1.)*height/2.;
    //     //     // int x1 = (v1.x*(l/(l+std::abs(v1.z+1.)))+1.)*width/2.;
    //     //     // int y1 = (v1.y*(l/(l+std::abs(v1.z+1.)))+1.)*height/2.;
            
    //     //     line(x0, y0, x1, y1, image, white);
            
    //     // }
    // }

    // image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    // std::cout << uvmap.get_width() << " " << uvmap.get_height() << std::endl;
    // image.write_tga_file("output.tga");
    delete model;
    return 0;
}

