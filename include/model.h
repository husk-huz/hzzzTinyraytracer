#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include <string>
#include "geometry.h"
#include "tgaimage.h"

class Model 
{
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> norms_;
	std::vector<Vec2f> uv_;
	std::vector<std::vector<int> > faces_;
	std::vector<std::vector<int> > texture_faces_;	
	std::vector<std::vector<int> > norm_vertexs_;
	std::vector<Vec3f> face_norm_;
	TGAImage* texture;

public:
	bool is_random_model;
	bool has_texture;
	Model(const char *filename);
	Model(float distance);
	// generate a model with a only one face randomly
	Model();
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	Vec2f uv(int i);
	Vec3f norm(int i);
	Vec3f face_norm(int i);
	std::vector<int> face(int idx);
	std::vector<int> texture_face(int idx);
	std::vector<int> norm_vertex(int idx);
	std::string toString();
	void calNorm();
	void rotate(float angle);
	void scale(float _scale);
	void move(Vec3f v);
	void setTexture(TGAImage* _texture);
	TGAImage* getTexture();
	TGAColor getTextureColor(float x, float y);
};

#endif //__MODEL_H__
