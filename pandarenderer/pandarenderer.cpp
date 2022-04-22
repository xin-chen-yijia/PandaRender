// pandarenderer.cpp: 定义应用程序的入口点。
//

#include "pandarenderer.h"
#include "TGAImage.h"
#include "panda_gl.h"
#include "model.h"
#include <memory>
#include <algorithm>

using namespace std;
using namespace panda;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const int width = 500;
const int height = 500;

shared_ptr<Model> headModel;
vec3f light_dir(0, 0, 1); // define light_dir
float zbuffer[width * height] = { 0 };

class GouraudShader : public IShader {
public:
	mat4x4 proj;
	mat4x4 view;

	TGAImage diffuse;
	vec3f varying_intensity;
	vec2f varying_uvs[3];

	vec4f vertex(int iface, int nthvert) override
	{
		varying_intensity[nthvert] = headModel->normal(iface, nthvert)* light_dir; // get diffuse lighting intensity
		varying_uvs[nthvert] = headModel->uv(iface, nthvert);
		auto v = headModel->vert(headModel->face(iface)[nthvert]);
		vec4f tmp = { v.x,v.y,v.z,1.0f };
		return  proj * view * tmp;
	}

	bool fragment(vec3f bar, TGAColor& color) override
	{
		vec2f uv=(varying_uvs[0] * bar[0] + varying_uvs[1] * bar[1] + varying_uvs[2] * bar[2]);
		uv.y = 1.0f - uv.y;

		float intensity = varying_intensity * bar;
		int dx = int(uv.x * diffuse.get_width());
		int dy = (int)(uv.y* diffuse.get_height());
		color = diffuse.get(dx, dy);//*intensity;
		//color = TGAColor(255, 0, 0, 255);
		return false;
	}
};

int main()
{
	for (int i = 0;i < width * height;++i)
	{
		zbuffer[i] = std::numeric_limits<float>::max();
	}

	TGAImage image(width, height, TGAImage::RGB);

	TGAImage diffuse;
	if (!diffuse.read_tga_file("Resource/african_head_diffuse.tga"))
		//if (!diffuse.read_tga_file("Resource/floor_diffuse.tga"))
	{
		printf("read diffuse fail...\n");
		return 1;
	}

	mat4x4 viewport = Viewport(width / 8.0f, height / 8.0f, width * 0.75f, height * 0.75f);
	mat4x4 view = Lookat({ 0.0f,0.0f,2.0f }, { 0,0.0f,0 }, { 0,1,0 });
	mat4x4 proj = Ortho(-1.0f, 1.0f, -1.0f, 1.0f, 0.0f, 100.0f);
	//mat4x4 proj = Frustum(-1.0f, 1.0f, -1.0f, 1.0f, 0.3f, 100.0f);
	//mat4x4 proj = Perspective(60.0,1.0f, 0.3f, 20.0f);

	GouraudShader gaudShader;
	gaudShader.proj = proj;
	gaudShader.view = view;
	gaudShader.diffuse = diffuse;

	headModel = make_shared<Model>("Resource/african_head.obj");
	//headModel = make_shared<Model>("Resource/floor.obj");
	for (int i = 0;i < headModel->nfaces();++i)
	{
		vec4f vexPoses[3];
		for (int j = 0;j < 3;++j)
		{
			vexPoses [j]= gaudShader.vertex(i, j);
		}

		triangle(vexPoses, zbuffer, viewport, &gaudShader, image);
	}
	image.write_tga_file("d:/Temp/a.tga",false,false);

	return 0;
}
