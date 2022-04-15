#pragma once

#include <vector>
#include <string>

#include "TGAImage.h"
#include "geomerty.h"

namespace panda
{
	class Model
	{
	private:
		std::vector<vec3f> verts_;
		std::vector<std::vector<vec3i> > faces_;
		std::vector<vec2f> uv_;
		std::vector<vec3f> norms_;

		TGAImage diffusemap_;
		TGAImage normalmap_;
		TGAImage specularmap_;

		void load_texture(std::string filename, const char* suffix, TGAImage& img);

	public:
		Model(const char* filename);
		int nverts();
		int nfaces();
		vec3f vert(int i);
		std::vector<int> face(int idx);
		vec2f uv(int iface, int nvert);
		vec3f normal(int iface, int nthvert);
		vec3f normal(vec2f uv);


	};
}
