#pragma once

#include "geomerty.h"
#include "TGAImage.h"

void line(int x0, int y0, int x1, int y1, panda::TGAImage& image, const panda::TGAColor& color);

void triangle(vec2i p0, vec2i p1, vec2i p2, panda::TGAImage& image, const panda::TGAColor& color);
void triangle(vec3* pts, float* zbuffer, panda::TGAImage& image, const panda::TGAColor& color);
void triangle(vec3* pts, vec2i* uvs, panda::TGAImage& inImage, float* zbuffer, panda::TGAImage& image, const panda::TGAColor& color);

vec3 barycentric(const vec3& A, const vec3& B, const vec3& C, const vec3& P);
vec3 barycentric(const vec2 tri[3], const vec2 P);
