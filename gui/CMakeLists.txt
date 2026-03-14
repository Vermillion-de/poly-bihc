cmake_minimum_required(VERSION 3.20)
project(bih LANGUAGES CUDA CXX)

set(CMAKE_BUILD_TYPE release)

# CUDA
enable_language("CUDA")
find_package(CUDA REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(ZLIB REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")
# Easy3d
add_subdirectory("./deps/Easy3D")
# OpenMP
find_package(OpenMP REQUIRED)

include_directories("./include/")
include_directories("./deps/eigen/")
include_directories("./deps/include/")

add_library(lib_viewer "./include/imgui_viewer.hpp" "./src/imgui_viewer.cpp")
target_link_libraries(lib_viewer easy3d::viewer easy3d::gui 3rd_imgui 3rd_glfw OpenMP::OpenMP_CXX)
target_compile_features(lib_viewer PUBLIC cxx_std_17)

add_executable(main "./main.cpp")
target_link_libraries(main easy3d::viewer lib_viewer ZLIB::ZLIB)