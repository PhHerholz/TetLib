#include <igl/readOFF.h>
#include <igl/writeOFF.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/colormap.h>
#include <imgui/imgui.h>
#include <iostream>

#include "lodepng.h"
#include "CGALTriangulation.hpp"
#include "CGALPolyhedron.hpp"
#include "TetgenMeshPolyhedron.hpp"
#include "SolveConstrained.hpp"

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;

void solveDirichletProblem(CGALTriangulation<Kernel>& tri, Eigen::MatrixXd& x)
{
    const int cntr = tri.centerVertex();
    auto constr = tri.surfaceVertices();
    
    const int n = tri.mesh.number_of_vertices();
    
    Eigen::MatrixXd b(n, 1);
    b.setZero();
    
    Eigen::MatrixXd C(constr.size() + 1, 1);
    C.setZero();
    C(constr.size()) = 1.;
    
    constr.push_back(cntr);

    Eigen::SparseMatrix<double> A, L, M;
    tri.DECLaplacian(L, &M);
    const double t = tri.meanEdgeLengthSquared();
    A = M - t * L;

    //Eigen::SparseMatrix<double> LTL = L.transpose() * L;
    //Eigen::MatrixXd LTb = L.transpose() * b;

    
    solveConstrainedSymmetric(A, b, constr, C, x);
    
   


    Eigen::MatrixXd pts;
    tri.getPoints(pts);
    
    const double div = sqrt(- 4 * t * log(1e-20));
    
    
    std::ofstream file("../x");
    
    for(int i = 0; i < x.size(); ++i)
    {
         file <<  (pts.row(i) - pts.row(cntr)).norm() << " " << x(i) << "\n";
    //    x(i) = (pts.row(i) - pts.row(cntr)).norm();
        x(i) = sqrt(- 4 * t * log((1. - 1e-20) * x(i) + 1e-20)) / div;
    
       
    }
    
    
    file.close();
    
    
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    tri.marchingTets(x, V, F, 0.5);

    igl::writeOFF("../iso.off", V, F);
}


void setTexture(const std::string filename, igl::opengl::glfw::Viewer& viewer)
{
    static Eigen::Matrix<unsigned char, -1, -1> R, G, B;
    static std::string currentFile;
    
    if(currentFile.compare(filename))
    {
        currentFile = filename;
        
        unsigned w, h;
        std::vector<unsigned char> image; //the raw pixels
        
        if( lodepng::decode(image, w, h, filename) ) return;
        
        R.resize(w, h);
        G.resize(w, h);
        B.resize(w, h);
        
        for(int i = 0; i < w * h; ++i)
        {
            R.data()[i] = image[4 * i];
            G.data()[i] = image[4 * i + 1];
            B.data()[i] = image[4 * i + 2];
        }
    }
    
    viewer.data().set_texture(R, G, B);
}

int main(int argc, char *argv[])
{
    CGALPolyhedron<Kernel> p;
    CGALTriangulation<Kernel> tri;
    
    p.load("../data/sphere.off");
    tetgenMeshPolyhedron(p, tri, 0.0001);
    
    Eigen::MatrixXd x;
    solveDirichletProblem(tri, x);
    
    std::cout << "done" << std::endl;
       
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    std::array<double, 4> plane{1,0,0,100};
    tri.cutMesh(plane, V, F);
    
    // Style
    Eigen::Vector3d ambient(.1,.1,.1);
    Eigen::Vector3d diffuse(.7,.7,.7);
    Eigen::Vector3d specular(.9,.9,.9);
    
    // Init the viewer
    igl::opengl::glfw::Viewer viewer;
        
    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
    
    // Customize the menu
    double offset = 0.0;
    int dir = 0;
    
    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        //   menu.draw_viewer_menu();
        
        // Add new group
        if (ImGui::CollapsingHeader("Cut View", ImGuiTreeNodeFlags_DefaultOpen))
        {
            double mi = -1.;
            double ma = 1.;
            double oldoffset = offset;
            ImGui::DragScalar("double", ImGuiDataType_Double, &offset, 0.1, &mi, &ma, "%.4f");
            
            int oldDir = dir;
            ImGui::Combo("Direction", (int *)(&dir), "X\0Y\0Z\0\0");
            
            if(oldoffset != offset || dir != oldDir)
            {
                std::array<double, 4> plane{0,0,0,offset};
                plane[dir] = 1;
                
                Eigen::MatrixXd V;
                Eigen::MatrixXi F;
                    
                auto ids = tri.cutMesh(plane, V, F);
             
                if(ids.size())
                {
                    viewer.data().clear();
                    viewer.data().set_mesh(V, F);
                    viewer.data().uniform_colors(ambient, diffuse, specular);
                    viewer.data().show_texture = 1;
                    
                    Eigen::VectorXd x2(ids.size());
                    
                    for(int i = 0; i < ids.size(); ++i)
                        x2(i) = x(ids[i]);
                    
                  //  Eigen::MatrixXd color;
                  //  igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
                    //viewer.data().set_colors(color);
               
                    setTexture("../data/tex.png", viewer);
                    Eigen::MatrixXd UV(x2.rows(), 2);
                    UV.col(0) = UV.col(1) = x2;
                    viewer.data().set_uv(UV);
                }
            }
        }
    };
    
    
    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.core().background_color.setOnes();
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = 1;
    setTexture("../data/tex.png", viewer);
    Eigen::MatrixXd UV(V.rows(), 2);
    viewer.data().set_uv(UV.setZero());
    
    viewer.launch();
}
