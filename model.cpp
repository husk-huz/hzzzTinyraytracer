#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "include/model.h"

Model::Model(const char *filename) : verts_(), faces_() 
{
    this->is_random_model = false;
    this->has_texture = false;
    std::cout << "Load model: " << filename << std::endl;
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) 
    {   
        std::cout << "fail to open" << std::endl;
        return;
    }
    std::string line;
    try
    {
        int count = 0;
        while (!in.eof()) 
        {
            std::getline(in, line);
            std::istringstream iss(line.c_str());
            char trash;
            if (!line.compare(0, 2, "v ")) 
            {
                iss >> trash;
                Vec3f v;
                for (int i=0;i<3;i++) iss >> v.raw[i];
                verts_.push_back(v);
            }
            else if (!line.compare(0, 2, "vt"))
            {
                iss >> trash >> trash;
                Vec2f uv;
                // for (int i=0;i<2;i++) iss >> uv.raw[i];
                iss >> uv.x >> uv.y;
                uv_.push_back(uv);
                // std:: cout << uv.x << " " << uv.y << std::endl;
            } 
            else if(!line.compare(0, 2, "vn"))
            {
                iss >> trash >> trash;
                Vec3f n;
                for (int i=0;i<3;i++) iss >> n.raw[i];
                norms_.push_back(n);
            }
            else if (!line.compare(0, 2, "f ")) 
            {
                std::vector<int> f;
                int t;
                iss >> trash;
                char c;

                // x1/x2 x3/x4 x5/x6          9
                // 0  2  3  5  6  8
                // x1/x2/x3 x4/x5/x6 x7/x8/x9 15
                // 0  2  4  5  7  9  10 12 14

                count = 0;
                while(iss >> c)
                {
                    if(isdigit(c))
                    {
                        iss.unget();
                        iss >> t;
                        f.push_back(t);
                        count++;
                    }
                }

                std::vector<int> f_v;
                std::vector<int> f_t;
                std::vector<int> f_n;

                if(count == 9)
                // with norm given
                {
                    f_v.push_back(f[0]-1);
                    f_v.push_back(f[3]-1);
                    f_v.push_back(f[6]-1);

                    f_t.push_back(f[1]-1);
                    f_t.push_back(f[4]-1);
                    f_t.push_back(f[7]-1);

                    f_n.push_back(f[2]-1);
                    f_n.push_back(f[5]-1);
                    f_n.push_back(f[8]-1);
                }
                else if(count == 6)
                // no norm given
                {
                    f_v.push_back(f[0]-1);
                    f_v.push_back(f[2]-1);
                    f_v.push_back(f[4]-1);

                    f_t.push_back(f[1]-1);
                    f_t.push_back(f[3]-1);
                    f_t.push_back(f[5]-1);
                }
                else 
                {
                    std::cout << "Face with count : " << count << std::endl;
                }

                faces_.push_back(f_v);
                texture_faces_.push_back(f_t);
                norm_vertexs_.push_back(f_n);
            }
        }
        std::cout << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;
        if(count == 6)
        {
            this->calNorm();
        }
        else if(count == 9)
        {
            norm_vertexs_.resize(faces_.size());
            face_norm_.resize(faces_.size());
            for(int i=0; i<norm_vertexs_.size(); i++)
            {
                norm_vertexs_[i].resize(3);
            }
            for(int i=0; i<faces_.size(); i++)
            {
                std::vector<int> face = faces_[i];
                for(int j=0; j<face.size(); j++)
                {
                    norm_vertexs_[i][j] = face[j];
                }

                face_norm_[i] = norms_[face[0]] + norms_[face[1]] + norms_[face[2]];
                face_norm_[i].normalize();
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        std::cout << "fail to load model" << std::endl;
    }

    std::cout << "load model success" << std::endl;
}

Model::Model()
{
    this->is_random_model = true;
    this->has_texture = false;
    Vec3f center;
    center.x = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.y = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.z = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.normalize();
    center = center;
    center.x += 0.0;
    center.y += 2.0;
    center.z += 1.0;

    Vec3f dir_1;
    dir_1.x = rand() % 1000 / 50000.0 - 0.01;
    dir_1.y = rand() % 1000 / 50000.0 - 0.01;
    dir_1.z = rand() % 1000 / 50000.0 - 0.01;

    Vec3f dir_2;
    dir_2.x = rand() % 1000 / 50000.0 - 0.01;
    dir_2.y = rand() % 1000 / 50000.0 - 0.01;
    dir_2.z = rand() % 1000 / 50000.0 - 0.01;

    Vec3f dir_3;
    dir_3.x = rand() % 1000 / 50000.0 - 0.01;
    dir_3.y = rand() % 1000 / 50000.0 - 0.01;
    dir_3.z = rand() % 1000 / 50000.0 - 0.01;

    verts_.push_back(center + dir_1);
    verts_.push_back(center + dir_2);
    verts_.push_back(center + dir_3);

    std::vector<int> face;
    face.push_back(0);
    face.push_back(1);
    face.push_back(2);

    faces_.push_back(face);
}

Model::Model(float distance)
{
    this->is_random_model = true;
    this->has_texture = false;
    Vec3f center;
    center.x = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.y = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.z = 1.0*(rand() % 1000 / 500.0 - 1.0);
    center.normalize();
    center = center * distance;

    Vec3f dir_1;
    dir_1.x = rand() % 1000 / 50000.0 - 0.01;
    dir_1.y = rand() % 1000 / 50000.0 - 0.01;
    dir_1.z = rand() % 1000 / 50000.0 - 0.01;

    Vec3f dir_2;
    dir_2.x = rand() % 1000 / 50000.0 - 0.01;
    dir_2.y = rand() % 1000 / 50000.0 - 0.01;
    dir_2.z = rand() % 1000 / 50000.0 - 0.01;

    Vec3f dir_3;
    dir_3.x = rand() % 1000 / 50000.0 - 0.01;
    dir_3.y = rand() % 1000 / 50000.0 - 0.01;
    dir_3.z = rand() % 1000 / 50000.0 - 0.01;

    verts_.push_back(center + dir_1);
    verts_.push_back(center + dir_2);
    verts_.push_back(center + dir_3);

    std::vector<int> face;
    face.push_back(0);
    face.push_back(1);
    face.push_back(2);

    faces_.push_back(face);
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() 
{
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

std::vector<int> Model::texture_face(int idx)
{
    return texture_faces_[idx];
}

std::vector<int> Model::norm_vertex(int idx)
{
    return norm_vertexs_[idx];
}

Vec2f Model::uv(int idx)
{
    return uv_[idx];
}

Vec3f Model::norm(int i) 
{
    return norms_[i];
}

Vec3f Model::vert(int i) 
{
    return verts_[i];
}

Vec3f Model::face_norm(int i)
{
    return face_norm_[i];
}

std::string Model::toString()
{
    std::string str = "model: \n";
    for(int i = 0; i < std::min(10, int(this->verts_.size())); i++)
    {
        str += "v ";
        str += std::to_string(verts_[i].x) + " ";
        str += std::to_string(verts_[i].y) + " ";
        str += std::to_string(verts_[i].z) + "\n";
    }
    for(int i = 0; i < std::min(10, int(faces_.size())); i++)
    {
        str += "f ";
        for(int j = 0; j < faces_[i].size(); j++)
        {
            str += std::to_string(faces_[i][j]) + " ";
        }
        str += "\n";
    }
    for(int i = 0; i < std::min(20, int(norms_.size())); i++)
    {
        str += "vn ";
        str += std::to_string(norms_[i].x) + " " +
               std::to_string(norms_[i].y) + " " +
               std::to_string(norms_[i].z);
        str += "\n";
    }
    return str;
}

void Model::calNorm()
{
    norm_vertexs_.resize(faces_.size());
    norms_.resize(verts_.size());
    face_norm_.resize(faces_.size());
    for(int i=0; i<norm_vertexs_.size(); i++)
    {
        norm_vertexs_[i].resize(3);
    }

    std::vector<std::vector<int> > vertex_faces;
    for(int i=0; i<faces_.size(); i++)
    {
        std::vector<int> face = faces_[i];
        for(int j=0; j<face.size(); j++)
        {
            int v = face[j]+1;
            if(vertex_faces.size() < v)
            {
                vertex_faces.resize(v);
            }
            vertex_faces[v-1].push_back(i);
        }
    }

    // std::cout << "vertex_faces size: " << vertex_faces.size() << std::endl;

    for (int i = 0; i < verts_.size(); i++)
    {
        for (int j = 0; j < vertex_faces[i].size(); j++)
        {
            std::vector<int> face = faces_[vertex_faces[i][j]];
            int x1 = face[0];
            int x2 = face[1];
            int x3 = face[2];
            // std::cout << x1 << " " << x2 << " " << x3 << "#";
            //计算顶点连接面的法向量
            float vx1 = verts_[x1].x - verts_[x2].x;
            float vy1 = verts_[x1].y - verts_[x2].y;
            float vz1 = verts_[x1].z - verts_[x2].z;

            float vx2 = verts_[x1].x - verts_[x3].x;
            float vy2 = verts_[x1].y - verts_[x3].y;
            float vz2 = verts_[x1].z - verts_[x3].z;

            float nx = vy1 * vz2 - vz1 * vy2;
            float ny = vz1 * vx2 - vx1 * vz2;
            float nz = vx1 * vy2 - vy1* vx2;

            //将连接面的法向量相加
            norms_[i].x += nx*10000.0;
            norms_[i].y += ny*10000.0;
            norms_[i].z += nz*10000.0;
            // std::cout << nx << " " << ny << " " << nz << "#";
        }
        
        norms_[i].normalize();
    }

    for(int i=0; i<faces_.size(); i++)
    {
        std::vector<int> face = faces_[i];
        for(int j=0; j<face.size(); j++)
        {
            norm_vertexs_[i][j] = face[j];
        }

        face_norm_[i] = norms_[face[0]] + norms_[face[1]] + norms_[face[2]];
        face_norm_[i].normalize();
    }

    std::cout << "Finish cal norm for vertexs" << std::endl;
}

void Model::rotate(float angle)
{
    for(int i = 0; i < verts_.size(); i++)
    {
        // verts_[i].rotate(angle);
        float x = verts_[i].x;
        float y = verts_[i].z;

        verts_[i].x = x * cos(angle) - y * sin(angle);
        verts_[i].z = x * sin(angle) + y * cos(angle);
    }
}

void Model::scale(float _scale)
{
    for(int i = 0; i < verts_.size(); i++)
    {
        verts_[i].x *= _scale;
        verts_[i].y *= _scale;
        verts_[i].z *= _scale;
    }
}

void Model::move(Vec3f dir)
{
    for(int i = 0; i < verts_.size(); i++)
    {
        verts_[i] = verts_[i] + dir;
    }
}

void Model::setTexture(TGAImage* _texture)
{
    this->texture = _texture;
    this->has_texture = true;
} 

TGAImage* Model::getTexture()
{
    return this->texture;
}

TGAColor Model::getTextureColor(float x, float y)
{
    if(this->has_texture)
    {
        return this->texture->get(x, y);
    }
    else 
    {
        return TGAColor(255, 255, 255, 255);
    }
}







