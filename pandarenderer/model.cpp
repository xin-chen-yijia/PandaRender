#include "model.h"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

using namespace panda;

Model::Model(const char* filename) : verts_(), faces_() {
    std::ifstream in;
    in.open(filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3 v;
            //for (int i = 0;i < 3;i++) iss >> v.raw[i];
            iss >> v.x;
            iss >> v.y;
            iss >> v.z;
            verts_.push_back(v);
        }
        else if (!line.compare(0, 2, "f ")) {
            std::vector<vec3i> f;
            int itrash;
            vec3i tmp;
            iss >> trash;
            while (iss >> tmp[0] >> trash >> tmp[1] >> trash >> tmp[2]) {
                //idx--; // in wavefront obj all indices start at 1, not zero
                for (int j = 0;j < 3;++j)
                {
                    tmp[j]--;
                }
                f.push_back(tmp);
            }
            faces_.push_back(f);
        }
        else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            vec2 uv;
            for (int i = 0;i < 2;i++) iss >> uv[i];
            uv_.push_back(uv);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# " << faces_.size() << std::endl;
}

int Model::nverts()
{
    return (int)verts_.size();
}

int Model::nfaces()
{
    return (int)faces_.size();
}

vec3 Model::vert(int i)
{
    assert(i >= 0 && i < verts_.size());
    return verts_[i];
}

std::vector<int> Model::face(int idx)
{
    assert(idx >= 0 && idx < faces_.size());
    std::vector<int> tmp;
    for (int i = 0;i < (int)faces_[idx].size();++i)
    {
        tmp.push_back(faces_[idx][i][0]);
    }
    return tmp;
}

vec2 Model::uv(int iface, int nvert) {
    int idx = faces_[iface][nvert][1];
    return uv_[idx];
}

