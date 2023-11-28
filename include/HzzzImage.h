#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "Ray.h"
#include "RayTraceable.h"

class HzzzImage
{
private:
    char* image;
    float* zbuffer;
    int width;
    int height;
    float screen_width  = 1.0f;
    float screen_height = 1.0f;
    Vec3f camera;
    Vec3f view_direction;

public:
    float fov = 90;
    HzzzImage(int width, int height);
    bool set(int x, int y, const TGAColor &color);
    bool set(int x, int y, float z, const TGAColor &color);
    bool setCamera(Vec3f camera);
    bool setViewDirection(Vec3f view_direction);
    Vec3f getCamera();
    Vec3f getViewDirection();
    void Line(Vec2i p0, Vec2i p1, TGAColor color);
    void Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, TGAImage& texture);
    // depth: 0~1
    void Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, TGAImage& texture, float depth);
    void Triangle(Vec3f p1, Vec3f p2, Vec3f p3, Vec2f t1, Vec2f t2, Vec2f t3, Model* model, float* depths);
    void RenderModel(Model* model, TGAImage& texture);
    void RenderModel(Model* model, TGAImage& texture, Vec3f light);
    void RenderModel(Model* model, Vec3f light);
    void RayTracingRender(Model* model, TGAImage& texture, float fov, std::vector<Model*> light);
    void RayTracingRender(RayTraceable* ray_traceable, float fov, std::vector<Model*> light);
    void RayTracingRender(std::vector<RayTraceable*> ray_traceables, std::vector<Model*> lights);
    void RayTracingRender(RayTracer* rt);
    TGAImage & toTGAImage(TGAImage &image);
    TGAImage toTGAImage();
    Vec2i toImageCoord(Vec3f v);
    void toFile(std::string filename);
    void saveToTGAImage(std::string filename);
    ~HzzzImage();
};

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color);

void rotateVec(Vec3f& v, float theta, float phi);
