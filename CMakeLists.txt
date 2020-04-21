cmake_minimum_required(VERSION 3.1)

project(tetlib)

set(CGAL_DIR=~/src/cgal)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"    ON)
option(LIBIGL_WITH_PNG               "Use PNG"      ON)


find_package(CGAL REQUIRED)
find_package(LIBIGL REQUIRED QUIET)
#find_package(Eigen3 REQUIRED)


file(GLOB LIB_SOURCES src/*.cpp )
list(REMOVE_ITEM LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_generation.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/tet_experiment.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh_viewer.cpp")

message(${LIB_SOURCES})

file( GLOB LIB_HEADERS src/*.h )
add_library( TetLib ${LIB_SOURCES} ${LIB_HEADERS} )
target_link_libraries( TetLib igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png CGAL::CGAL)

add_executable(mesh_generation "src/mesh_generation.cpp")
target_link_libraries(mesh_generation TetLib igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png CGAL::CGAL)

add_executable(tet_experiment "src/tet_experiment.cpp")
target_link_libraries(tet_experiment TetLib igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png CGAL::CGAL)

add_executable(mesh_viewer "src/mesh_viewer.cpp")
target_link_libraries(mesh_viewer TetLib igl::core igl::opengl_glfw igl::opengl_glfw_imgui igl::png CGAL::CGAL)
