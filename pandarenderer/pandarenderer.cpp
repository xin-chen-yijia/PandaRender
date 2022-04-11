// pandarenderer.cpp: 定义应用程序的入口点。
//

#include "pandarenderer.h"
#include "TGAImage.h"
#include "panda_gl.h"
#include "model.h"
#include <memory>

using namespace std;
using namespace panda;

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

class A {
public:
	A() { printf("A cons\n"); }
};

int main()
{
	const int width = 500;
	const int height = 500;
	TGAImage image(width, height, TGAImage::RGB);
	//line(100, 100, 200, 450, image, white);
	//line(500, 500, 100, 100, image, white);

	//triangle(vec2i(100, 100), vec2i(100, 300), vec2i(300, 200), image, red);
	//triangle(vec2i(300, 300), vec2i(450, 450), vec2i(400, 250), image, red);
	//triangle(vec2i(100, 300), vec2i(300, 300), vec2i(250, 400), image, white);

	TGAImage diffuse;
	if (!diffuse.read_tga_file("Resource/african_head_diffuse.tga"))
	{
		printf("read diffuse fail...\n");
		return 1;
	}

	vec3 light_dir(0, 0, -1); // define light_dir
	float zbuffer[width * height] = { 0 };
	//light_dir.normalize();

	auto worldToScreenPos = [width,height](vec3 a) {
		return vec3(int((a.x + 1.0f) * width * 0.5f + 0.5f), int((a.y + 1.0f) * height * 0.5f + 0.5f), a.z);
	};

	Model head("Resource/african_head.obj");
	for (int i = 0;i < head.nfaces();++i)
	{
		//for (int j = 0;j < 3;++j)
		//{
		//	auto v0 = head.vert(head.face(i)[j]);
		//	auto v1 = head.vert(head.face(i)[(j+1)%3]);

		//	int x0 = (v0.x + 1.0) * width / 2;
		//	int y0 = (v0.y + 1.0) * height/ 2;
		//	int x1 = (v1.x + 1.0) * width / 2;
		//	int y1 = (v1.y + 1.0) * height/ 2;

		//	line(x0, y0, x1, y1, image, white);
		//}

		vec2i points[3];
		vec3 worldPoses[3];
		for (int j = 0;j < 3;++j)
		{
			auto v0 = head.vert(head.face(i)[j]);
			int x0 = (v0.x + 1.0) * width / 2.0;
			int y0 = (v0.y + 1.0) * height/ 2.0;
			points[j] = vec2i(x0, y0);
			worldPoses[j] = worldToScreenPos(v0);
		}

		vec3 A = head.vert(head.face(i)[1]) - head.vert(head.face(i)[0]);
		vec3 B = head.vert(head.face(i)[2]) - head.vert(head.face(i)[0]);
		
		vec3 normal = cross(B, A);
		normal.normalize();
		float intensity = normal * light_dir;
		if (intensity < 0)
		{
			intensity = 0.0f;
		}

		//uv

		vec2i uvs[3];
		for (int j = 0;j < 3;++j)
		{
			auto tmp = head.uv(i, j);
			uvs[j].x = tmp.x * diffuse.get_width();
			uvs[j].y = (1.0-tmp.y) * diffuse.get_height();
		}


		triangle(worldPoses, uvs,diffuse, zbuffer, image, TGAColor(intensity * 255, intensity* 255, intensity * 255,255));
	}

	image.write_tga_file("d:/Temp/a.tga",true,false);

	return 0;
}
