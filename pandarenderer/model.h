#pragma once

#include <vector>
#include "geomerty.h"

namespace panda
{
	class Model
	{
	private:
	private:
		std::vector<vec3> verts_;
		std::vector<std::vector<vec3i> > faces_;
		std::vector<vec2> uv_;

	public:
		Model(const char* filename);
		int nverts();
		int nfaces();
		vec3 vert(int i);
		std::vector<int> face(int idx);
		vec2 uv(int iface, int nvert);
	};
}
