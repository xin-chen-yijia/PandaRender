#include "panda_gl.h"

void line(int x0, int y0, int x1, int y1, panda::TGAImage& image, const panda::TGAColor& color)
{
	if (abs(y1 - y0) > abs(x1 - x0))
	{
		if (y0 > y1)
		{
			std::swap(y0, y1);
			std::swap(x0, x1);
		}
		int dx = x1 - x0;
		int dy = y1 - y0;
		int delta = std::abs(dx) * 2;
		int accDelta = 0;
		int t = x0;
		for (int i = y0;i <= y1;++i)
		{
			image.set(t, i, color);
			
			accDelta += delta;
			if (accDelta > dy)
			{
				accDelta -= dy * 2;
				t += (x0 < x1 ? 1 : -1);
			}
		}
	}
	else
	{
		if (x0 > x1)
		{
			std::swap(x0, x1);
			std::swap(y0, y1);
		}
		int dx = x1 - x0;
		int dy = y1 - y0;
		int delta = std::abs(dy) * 2;
		int accDelta = 0;
		int t = y0;
		for (int i = x0;i <= x1;++i)
		{
			image.set(i, t, color);

			accDelta += delta;
			if (accDelta > dx)
			{
				accDelta -= dx * 2;
				t += (y0 < y1 ? 1 : -1);
			}
		}
	}

}

void triangle(vec2i p0, vec2i p1, vec2i p2, panda::TGAImage& image, const panda::TGAColor& color)
{
	if (p0.y == p1.y && p0.y == p2.y) return;
	if (p0.y > p1.y) std::swap(p0, p1);
	if (p0.y > p2.y) std::swap(p0, p2);
	if (p1.y > p2.y) std::swap(p1, p2);

	int total_height = p2.y - p0.y; //p0 < p1 < p2
	for (int y = 0;y < total_height;++y)
	{
		bool is_second_half = y > (p1.y - p0.y) || p0.y == p1.y;

		float k1 = y / (float)total_height;
		float k2 = is_second_half ? (y - (p1.y-p0.y)) / float(p2.y - p1.y) : (y) / (float)(p1.y - p0.y);
		vec2i left = p0 +(p2 - p0) * k1;
		vec2i right = is_second_half ? (p1 + (p2 - p1) * k2) : (p0 + (p1 - p0) * k2);

		if (left.x > right.x) std::swap(left, right);
		for (int x = left.x;x < right.x;++x)
		{
			image.set(x, p0.y + y, color);
		}
	}


}

void triangle(vec3* pts, float* zbuffer, panda::TGAImage& image, const panda::TGAColor& color)
{
	vec2 boxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	vec2 boxmax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
	vec2 clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0;i < 3;++i)
	{
		boxmin.x = std::max(0.0f, std::min(boxmin.x, pts[i].x));
		boxmin.y = std::max(0.0f, std::min(boxmin.y, pts[i].y));

		boxmax.x = std::min(clamp.x, std::max(boxmax.x, pts[i].x));
		boxmax.y = std::min(clamp.y, std::max(boxmax.y, pts[i].y));
	}


	vec3 P;
	for (P.x = boxmin.x;P.x <= boxmax.x;++P.x)
	{
		for (P.y = boxmin.y;P.y <= boxmax.y;++P.y)
		{
			vec3 coor = barycentric(pts[0], pts[1], pts[2], P);
			if (coor.x < 0 || coor.y < 0 || coor.z < 0) continue;

			P.z = 0;
			for (int i = 0;i < 3;++i)
			{
				P.z += pts[i][2] * coor[i];
			}
			int idx = int(P.x + P.y * image.get_width());
			if (zbuffer[idx] < P.z)
			{
				zbuffer[idx] = P.z;
				image.set(P.x, P.y, color);
			}
		}
	}
}

void triangle(vec3* pts, vec2i* uvs, panda::TGAImage& inImage, float* zbuffer, panda::TGAImage& image, const panda::TGAColor& color)
{
	vec2 boxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	vec2 boxmax(std::numeric_limits<float>::min(), std::numeric_limits<float>::min());
	vec2 clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0;i < 3;++i)
	{
		boxmin.x = std::max(0.0f, std::min(boxmin.x, pts[i].x));
		boxmin.y = std::max(0.0f, std::min(boxmin.y, pts[i].y));

		boxmax.x = std::min(clamp.x, std::max(boxmax.x, pts[i].x));
		boxmax.y = std::min(clamp.y, std::max(boxmax.y, pts[i].y));
	}


	vec3 P;
	for (P.x = boxmin.x;P.x <= boxmax.x;++P.x)
	{
		for (P.y = boxmin.y;P.y <= boxmax.y;++P.y)
		{
			vec3 coor = barycentric(pts[0], pts[1], pts[2], P);
			if (coor.x < 0 || coor.y < 0 || coor.z < 0) continue;

			P.z = 0;
			for (int i = 0;i < 3;++i)
			{
				P.z += pts[i][2] * coor[i];
			}

			//uv
			vec2i uv(0, 0);
			for (int i = 0;i < 3;++i)
			{
				uv.x += uvs[i][0] * coor[i];
				uv.y += uvs[i][1] * coor[i];
			}

			int idx = int(P.x + P.y * image.get_width());
			if (zbuffer[idx] < P.z)
			{
				zbuffer[idx] = P.z;
				image.set(P.x, P.y, inImage.get(uv.x,uv.y));
			}
		}
	}
}

vec3 barycentric(const vec2 tri[3], const vec2 P) {
	mat<3, 3> ABC = { {embed<Float,3>(tri[0]), embed<Float,3>(tri[1]), embed<Float,3>(tri[2])} };
	if (ABC.det() < 1e-3) return vec3(-1, 1, 1); // for a degenerate triangle generate negative coordinates, it will be thrown away by the rasterizator
	return ABC.invert_transpose() * embed<Float,3>(P);
}

vec3 barycentric(const vec3& A, const vec3& B, const vec3& C, const vec3& P) {
	vec3 s[2];
	for (int i = 2; i--; ) {
		s[i][0] = C[i] - A[i];
		s[i][1] = B[i] - A[i];
		s[i][2] = A[i] - P[i];
	}
	vec3 u = cross(s[0], s[1]);
	if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return vec3(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return vec3(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

mat4x4 Ortho(float bottom, float top, float left, float right, float near, float far)
{
	mat4x4 tmp = { {{2/(right-left),0,0,-(right+left)/(right-left)},   {0,2/(top-bottom),0,-(top+bottom)/(top-bottom)},   {0,0,-2/(far-near),-(far+near)/(far-near)},   {0,0,0,1.0f}} };
	return tmp;
}

mat4x4 Frustum(float bottom, float top, float left, float right, float near, float far)
{
	mat4x4 tmp = { {{2*near / (right - left),0,(right + left) / (right - left), 0},   {0,2*near / (top - bottom),(top + bottom) / (top - bottom),0},   {0,0,-(far + near) / (far - near),-2*far*near / (far - near)},   {0,0,-1,0}} };
	return tmp;
}

mat4x4 Perspective(float fovy, float aspect, float zNear, float zFar)
{
	float half = fovy * 0.5f;
	float h = tan(half / 180.0f * 3.141592653) * zNear;
	float w = h * aspect;

	return Frustum(-h, h, -w, w, zNear, zFar);
}

mat4x4 lookat(const vec3 eye, const vec3 center, const vec3 up) { // check https://github.com/ssloy/tinyrenderer/wiki/Lesson-5-Moving-the-camera
	vec3 z = (center - eye).normalize();
	vec3 x = cross(up, z).normalize();
	vec3 y = cross(z, x).normalize();
	mat4x4 Minv = { {{x.x,x.y,x.z,0.0f},   {y.x,y.y,y.z,0.0f},   {z.x,z.y,z.z,0.0f},   {0,0,0,1.0f}} };
	mat4x4 Tr = { {{1,0,0,-eye.x}, {0,1,0,-eye.y}, {0,0,1,-eye.z}, {0,0,0,1}} };
	return  Minv * Tr;
}
