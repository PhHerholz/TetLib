#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/colormap.h>

#include <igl/png/writePNG.h>
#include <igl/png/readPNG.h>

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
typedef Kernel::Point_3 Point;


#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    

void solveDirichletProblem(CGALTriangulation<Kernel>& tri, Eigen::MatrixXd& x)
{
    const int cntr = tri.centerVertex();
    auto constr = tri.surfaceVertices();
    
    const int n = tri.mesh.number_of_vertices();
    
    constr.push_back(cntr);
    
    Eigen::SparseMatrix<double> A, L, L2, M;
    tri.massMatrix(M);
    tri.FEMLaplacian(L);
    tri.DECLaplacian(L2);
    
    {
        constr.pop_back();
        
        std::ofstream file("../constr");
        for(int i : constr) file << i << "\n";
        file.close();
        
        Eigen::saveMarket(L, "../LFEM.mtx");
        Eigen::saveMarket(L2, "../LDEC.mtx");
    }
    
    
    const double t = tri.meanEdgeLengthSquared();
    A = M + t * L;
    
    //Eigen::SparseMatrix<double> LTL = L.transpose() * L;
    //Eigen::MatrixXd LTb = L.transpose() * b;
    
    Eigen::VectorXd bb(A.cols());
    bb.setZero();
    bb(cntr) = 1.;
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol;
    chol.analyzePattern(A);
    chol.factorize(A);
    Eigen::VectorXd h = chol.solve(bb);
    
    {
        double mi = h.minCoeff();
        double ma = h.maxCoeff();
        
        for(int i = 0; i < h.size(); ++i)
            h(i) = (h(i) - mi) / (ma - mi);
        
        x = h;
    }
    
    //solveConstrainedSymmetric(A, b, constr, C, x);
    
    Eigen::MatrixXd pts;
    tri.getPoints(pts);
    
    const double div = sqrt(- 4 * t * log(1e-20));
    
    std::ofstream file("../x");
    
    for(int i = 0; i < x.size(); ++i)
    {
        //    x(i) = (pts.row(i) - pts.row(cntr)).norm();
        double y = sqrt(- 4 * t * log((1. - 1e-20) * x(i) + 1e-20)) / div;
        file <<  (pts.row(i) - pts.row(cntr)).norm() << " " << x(i) << " " << y << "\n";
        
        x(i) = y;
    }
    
    
    file.close();
}


void screenshot(igl::opengl::glfw::Viewer& viewer, std::string filename) {

    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

    // Draw the scene in the buffers
    viewer.core().draw_buffer(
      viewer.data(),false,R,G,B,A);

    // Save it to a PNG
    igl::png::writePNG(R,G,B,A, "out/" + filename +  "out.png");

}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
	/*
	  std::cout << FILENAME_base << std::endl;
	  screenshot(viewer, FILENAME);
	  if (metric_shown == minangle) {
		  metric_shown = amips;
		  FILENAME = FILENAME_base + "amips";
	  } else if (metric_shown == amips) {
		  metric_shown = volume;
		  FILENAME = FILENAME_base + "volume";
	  } else {
		  metric_shown = minangle;
		  FILENAME = FILENAME_base + "minangle";
	  }

	  facecolors.resize(faceids.size(), 3);
	  for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
	  viewer.data().set_colors(facecolors);

	  if(metric_shown == minangle) {
		viewer.launch_shut(); 
	  } 
	  */
  }
  
}

int 
loadMeshWithOrbitpoints(CGALTriangulation<Kernel> &tri, std::vector<int> &orbitinds, std::string filepath)
{
	// load Triangulation from file and read orbitpoints
	tri.read(filepath);
	std::cout << tri.mesh.number_of_vertices() << std::endl;

	typedef CGALTriangulation<Kernel>::Vertex_handle Vertex_handle;
	typedef CGALTriangulation<Kernel>::Point Point;

	Point origin(0.,0.,0.);
	int originind = -1;
	
	for (Vertex_handle vh: tri.mesh.finite_vertex_handles()) {
		double dist = sqrt(CGAL::squared_distance(vh->point(), origin));		
		std::cout << std::endl << vh->point() << std::endl;
		std::cout << dist << std::endl;

		if (dist < (1e-10)) {
		// origin
			if (originind < 0) {
				originind = vh->info();
			} else {
				std::cout << "ERROR: File containts two origins" << std::endl;	
				return 0;
			}
		}
		if (fabs(dist - 0.5) < (1e-6)) {
		// orbitpoints
			orbitinds.push_back(vh->info());
		}
	}

	std::cout << "Loaded file with " << orbitinds.size() << " orbitpoints" << std::endl;
	orbitinds.push_back(originind);
	return 1;
}


int main(int argc, char *argv[])
{
	CGALTriangulation<Kernel> tri;
	int originind;
	std::vector<int> orbitinds;
	std::string filepath = "out/run_00/plots/Sphere_0.5_2_0_0_0_0_.meshfile";

	if(loadMeshWithOrbitpoints(tri, orbitinds, filepath)) {
		originind = orbitinds.back();
		orbitinds.pop_back();
		std::cout << "Oritin ind: " << originind << std::endl;
		std::cout << "Orbitinds: " << std::endl;
		for(int i: orbitinds) std::cout << i << " " << std::endl;
	
	} else {
		std::cout << "Something went wrong loading the mesh" << std::endl;	
	}


	// #########################################
	std::cout << "METRICS"  << std::endl;
	// #########################################
	enum Metric {minangle=0, amips, volume};
	Eigen::MatrixXd facecolors;
	std::map<Metric, Eigen::MatrixXd> cellcolors;
	std::vector<int>  faceids; 
	Metric metric_shown;


	std::cout << "Calc metrics" << std::endl;
	metric_shown = minangle;
	std::map<Metric, Eigen::VectorXd> cell_metrics;
	std::map<Metric, std::string> metric_names;
	metric_names[minangle] = "minangle";
	metric_names[amips] = "amips";
	metric_names[volume] = "volume";

	Eigen::VectorXd Vol, Minang, Amips; 
	tri.calcVolumeAllCells(Vol);
	tri.calcMinAngleAllCells(Minang);
	tri.calcAMIPSAllCells(Amips);
	cell_metrics[volume]=Vol;
	cell_metrics[minangle]=Minang;
	cell_metrics[amips]=Amips;

	bool normalize=false;
	double amips_max = 100;
	if (!normalize) {
		//normalize minangle by 70.5 deg
		for (int i=0; i < cell_metrics[minangle].size(); i++) cell_metrics[minangle][i] = cell_metrics[minangle][i] / 70.5;
		// normalize amips using the heuristically chosen max val amips_max and the min value 3
		for (int i=0; i < cell_metrics[amips].size(); i++) cell_metrics[amips][i] = (cell_metrics[amips][i] - 3) / amips_max;
	} 


	/* ################## SHOW  MESH ####################*/


	Eigen::MatrixXd cellcolors_volume; 
	igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[volume], true, cellcolors_volume);
	cellcolors[volume] = cellcolors_volume;
	Eigen::MatrixXd cellcolors_minangle; 
	igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[minangle], normalize, cellcolors_minangle);
	cellcolors[minangle] = cellcolors_minangle;
	Eigen::MatrixXd cellcolors_amips; 
	igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[amips], normalize, cellcolors_amips);
	cellcolors[amips] = cellcolors_amips;

    Point origin = Point(0., 0., 0.);

    Eigen::MatrixXd x;
    x.resize(tri.mesh.number_of_vertices(), 1);
    x.setZero();
  //  solveDirichletProblem(tri, x);
    
    std::cout << "done" << std::endl;


	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	double offset = 0.1;

	std::array<double, 4> plane{1,0,0,offset};
	std::vector<std::vector<int>> cmreturn = tri.cutMesh(plane, V, F);
	std::vector<int>  ids     = cmreturn[0];
	faceids = cmreturn[1];
	facecolors.resize(faceids.size(), 3);
	for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);

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
	float iso = 0.5f;
	bool orientation = false;
	int dir = 0;

	Metric metric = minangle;
	static bool showValues = false;

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		
		if (ImGui::CollapsingHeader("Presentation", ImGuiTreeNodeFlags_DefaultOpen))
		{
			Metric oldmetric = metric_shown;
			ImGui::Combo("Metric", (int *)(&metric_shown), "MinAngle\0AMIPS\0Volume\0\0");

			/*
			if(ImGui::Checkbox("Show Metric Values", &showValues)) {
				for(int i = 0; i < F.rows(); ++i) {
					const Eigen::Vector3d FaceCenter( (V(F(i,0), 0)+ V(F(i,1), 0)+ V(F(i,2), 0))/ 3.,  (V(F(i,0), 1)+ V(F(i,1), 1)+ V(F(i,2), 1))/ 3.,  (V(F(i,0), 2)+ V(F(i,1), 2)+ V(F(i,2), 2))/ 3.);
					viewer.data().add_label(FaceCenter,std::to_string(cell_metrics[metric][faceids[i]]));
				}
			} 
			*/
			//else {
			//	viewer.data().clear_labels();	
			//}

			if (oldmetric != metric_shown){
				facecolors.resize(faceids.size(), 3);
				for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
				viewer.data().set_colors(facecolors);
			}
		}
		
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
				
				cmreturn = tri.cutMesh(plane, V, F);
				ids = cmreturn[0];
				faceids = cmreturn[1];
				
				if(ids.size())
				{
					viewer.data().clear();
					viewer.data().set_mesh(V, F);

					facecolors.resize(faceids.size(), 3);
					for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
					viewer.data().set_colors(facecolors);
				}
			}
		}
	};

	// Plot the mesh
	viewer.data().set_mesh(V, F);
	viewer.data().set_colors(facecolors);
	viewer.core().background_color.setOnes();
	viewer.callback_key_down = &key_down;
		
	viewer.launch();
}
