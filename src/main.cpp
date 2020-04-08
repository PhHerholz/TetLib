#include <igl/readOFF.h>
#include <igl/writeOFF.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/colormap.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

#include <imgui/imgui.h>
#include <iostream>

#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/ArpackSupport>
#include "lodepng.h"
#include "CGALTriangulation.hpp"
#include "CGALPolyhedron.hpp"
#include "TetgenMeshPolyhedron.hpp"
#include "CGALMeshPolyhedron.hpp"
#include "SolveConstrained.hpp"


#include "tetgen.h"
#include "IndexedTetMesh.hpp"

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> Kernel;


#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    

/*
  We need a functor that can pretend it's const,
  but to be a good random number generator 
  it needs mutable state.
*/
namespace Eigen {
namespace internal {
template<typename Scalar> 
struct scalar_normal_dist_op 
{
  static boost::mt19937 rng;    // The uniform pseudo-random algorithm
  mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

  EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

  template<typename Index>
  inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;

template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };
} // end namespace internal
} // end namespace Eigen

/*
  Draw nn samples from a size-dimensional normal distribution
  with a specified mean and covariance
*/
Eigen::MatrixXd randomPoints(int nn) 
{
  int size = 3; // Dimensionality (rows)
  Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
  Eigen::internal::scalar_normal_dist_op<double>::rng.seed(1); // Seed the rng

  // Define mean and covariance of the distribution
  Eigen::VectorXd mean(size);       
  Eigen::MatrixXd covar(size,size);

  mean  <<  0,  0, 0;
  covar <<  1, 0, 0,
            0, 1, 0,
            0, 0, 1;

  Eigen::MatrixXd normTransform(size,size);

  Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

  // We can only use the cholesky decomposition if 
  // the covariance matrix is symmetric, pos-definite.
  // But a covariance matrix might be pos-semi-definite.
  // In that case, we'll go to an EigenSolver
  if (cholSolver.info()==Eigen::Success) {
    // Use cholesky solver
    normTransform = cholSolver.matrixL();
  } else {
    // Use eigen solver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    normTransform = eigenSolver.eigenvectors() 
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::MatrixXd samples = (normTransform 
                           * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise() 
                           + mean;
  return samples;
}

void setTexture(const std::string filename, igl::opengl::ViewerData& viewerData)
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
    
    viewerData.set_texture(R, G, B);
}

void
tetgenMeshSphere(CGALTriangulation<Kernel>& tri, int num_samples=500)
{
    using namespace std;
    using namespace std::chrono;

	const double area = .1;

	std::cout << "Generate random samples" << std::endl;
	Eigen::MatrixXd randomsamples = randomPoints(num_samples);
    
    tetgenio in, out, out0;
    
    in.firstnumber = 0;

    in.numberofpoints = num_samples + 1; 
    in.pointlist = new REAL[in.numberofpoints * 3];

	// add origin
	in.pointlist[0]=0;
	in.pointlist[1]=0;
	in.pointlist[2]=0;

	// add normalized samples (on 3d shpere)
	for(int i=0; i< randomsamples.cols(); i++){
		for(int j=0; j<3; j++){
			in.pointlist[3+3*i+j] = randomsamples(j,i) / randomsamples.col(i).norm();
		}
	}

	// generate simple mesh
	std::cout << "Step 1: generate starmesh" << std::endl;
    tetgenbehavior settings;
    string opts = string(""); //string("q1.414a") + to_string(area);
    settings.parse_commandline((char*)opts.c_str());
    settings.quiet = 1;
    tetrahedralize(&settings, &in, &out0);
	std::cout << "...done" << std::endl;


	std::cout << "Step 2: optimize" << std::endl;
    opts = string("rq1.414"); //string("q1.414a") + to_string(area);
    settings.parse_commandline((char*)opts.c_str());
    settings.quiet = 1;
    tetrahedralize(&settings, &out0, &out);
    
    
    IndexedTetMesh mesh;
    
    mesh.vertices.resize(out.numberofpoints);
    copy_n(out.pointlist, 3 * out.numberofpoints, (double*)mesh.vertices.data() );
    
    mesh.tets.resize(out.numberoftetrahedra);
    copy_n(out.tetrahedronlist, 4 * out.numberoftetrahedra, (int*)mesh.tets.data() );
    
    mesh.convert(tri);
}


int main(int argc, char *argv[])
{
    CGALPolyhedron<Kernel> p;
    CGALTriangulation<Kernel> tri;
    
    //p.load("../data/bunny.off");
   // tetgenMeshPolyhedron(p, tri, 0.0001);
    //meshPolyhedron(p, tri, 0.01);
   // tutte3d(tri);
    //  return 0;
	
	std::cout << "START" << std::endl;
	tetgenMeshSphere(tri);
	std::cout << "Finished tetgen" << std::endl;

    Eigen::MatrixXd x;
    x.resize(tri.mesh.number_of_vertices(), 1);
    x.setZero();
  //  solveDirichletProblem(tri, x);
    
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
    int cutMeshId = -1;
    int isoMeshId = -1;
    float iso = 0.5f;
    bool orientation = false;
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
            ImGui::DragScalar("cut offset", ImGuiDataType_Double, &offset, 0.1, &mi, &ma, "%.4f");
            
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
                    viewer.data(cutMeshId).clear();
                    viewer.data(cutMeshId).set_mesh(V, F);
                    viewer.data(cutMeshId).uniform_colors(ambient, diffuse, specular);
                    viewer.data(cutMeshId).show_texture = 1;
                    
                    Eigen::VectorXd x2(ids.size());
                    
                    for(int i = 0; i < ids.size(); ++i)
                        x2(i) = x(ids[i]);
                    
                    //  Eigen::MatrixXd color;
                    //  igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
                    //viewer.data().set_colors(color);
                    
                    setTexture("../data/tex.png", viewer.data(cutMeshId));
                    Eigen::MatrixXd UV(x2.rows(), 2);
                    UV.col(0) = UV.col(1) = x2;
                    viewer.data(cutMeshId).set_uv(UV);
                }
            }
        }
        
        
        if (ImGui::CollapsingHeader("Iso surface", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if(ImGui::SliderFloat("iso value", &iso, 0.f, 1.f, "%.4f") || ImGui::Checkbox("orientation", &orientation) )
            {
                Eigen::MatrixXd Vi;
                Eigen::MatrixXi Fi;
                
                tri.marchingTets(x, Vi, Fi, iso);
                
                if(orientation)
                {
                    for(int i = 0; i < Fi.rows(); ++i)
                        std::swap(Fi(i, 1), Fi(i, 2));
                }
                
                viewer.data(isoMeshId).clear();
                viewer.data(isoMeshId).set_mesh(Vi, Fi);
                viewer.data(isoMeshId).uniform_colors(ambient, diffuse, specular);
            }
            
            if(ImGui::Button("clear iso surface"))
            {
                viewer.data(isoMeshId).clear();
            }
        }
    };
    
    
    // Plot the mesh
    viewer.data().set_mesh(V, F);
    viewer.core().background_color.setOnes();
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_texture = 1;
    setTexture("../data/tex.png", viewer.data());
    Eigen::MatrixXd UV(V.rows(), 2);
    viewer.data().set_uv(UV.setZero());
    cutMeshId = viewer.selected_data_index;
    
    isoMeshId = viewer.append_mesh();
    
    viewer.launch();
}
