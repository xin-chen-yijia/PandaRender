#pragma once

#include "geomerty.h"
#include "TGAImage.h"

namespace panda
{
	class IShader
	{
	public:
		virtual ~IShader() {}
		virtual vec4f vertex(int iface, int nthvert) = 0;
		virtual bool fragment(vec3f bar, TGAColor& color) = 0;
	};

	void line(int x0, int y0, int x1, int y1, panda::TGAImage& image, const panda::TGAColor& color);

	void triangle(vec2i p0, vec2i p1, vec2i p2, panda::TGAImage& image, const panda::TGAColor& color);
	void triangle(vec3f* pts, float* zbuffer, panda::TGAImage& image, const panda::TGAColor& color);
	void triangle(vec3f* pts, vec2i* uvs, panda::TGAImage& inImage, float* zbuffer, panda::TGAImage& image);
	void triangle(vec4f* clip_verts, float* zbuffer, const mat4x4& viewPort, IShader* shader, panda::TGAImage& image);

	vec3f barycentric(const vec3f& A, const vec3f& B, const vec3f& C, const vec3f& P);
	vec3f barycentric(const vec2f tri[3], const vec2f& P);

	mat4x4 Ortho(float bottom, float top, float left, float right, float near, float far);
	mat4x4 Frustum(float bottom, float top, float left, float right, float near, float far);
	mat4x4 Perspective(float fovy, float aspect, float zNear, float zFar);
	mat4x4 Lookat(const vec3f eye, const vec3f center, const vec3f up);
	mat4x4 Viewport(int x, int y, int w, int h);
}// namespace panda
