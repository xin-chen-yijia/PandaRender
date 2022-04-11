# CMakeList.txt: pandarenderer 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 将源代码添加到此项目的可执行文件。
add_executable (pandarenderer "pandarenderer.cpp" "pandarenderer.h" "TGAImage.h" "TGAImage.cpp" "geomerty.h" "geometry.cpp" "panda_gl.h" "panda_gl.cpp" "model.h" "model.cpp")

# TODO: 如有需要，请添加测试并安装目标。
