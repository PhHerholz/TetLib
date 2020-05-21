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
#include <random>

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


void solveHeatProblem(CGALTriangulation<Kernel>& tri, Eigen::MatrixXd& h_fem, Eigen::MatrixXd& h_dec)
{
	const int cntr = tri.centerVertex();
	const int n = tri.mesh.number_of_vertices();

	// Construct A 
	Eigen::SparseMatrix<double> A_fem, A_dec, L_fem, L_dec, M;
	//tri.massMatrix(M);
	tri.FEMLaplacian(L_fem);

	std::cout << "Using Mixed DEC" << std::endl;
	tri.DECLaplacianMixed(L_dec, &M);
	// tri.DECLaplacian(L_dec, &M);
	const double t = tri.meanEdgeLengthSquared();
	A_fem = M + t * L_fem;
	A_dec = M - t * L_dec; 

	// solve the constrained problems
	std::vector<int> boundary_indices = tri.surfaceVertices(); 
	boundary_indices.push_back(cntr);

	Eigen::MatrixXd B(n, 1); B.setZero();
	Eigen::MatrixXd constrValues(boundary_indices.size(), 1);
	constrValues.setZero();
	constrValues(boundary_indices.size()-1, 0) = 1;

	std::cout << "SHAPES " <<std::endl;
	std::cout << constrValues.size() << std::endl;
	std::cout << B.rows() << ", " << B.cols()  << std::endl;
	std::cout << "FEM: " << std::endl;
	std::cout << A_fem.rows() << ", " << A_fem.cols()  << std::endl;
	std::cout << "DEC: " << std::endl;
	std::cout << A_dec.rows() << ", " << A_dec.cols()  << std::endl;

	solveConstrainedSymmetric(A_fem, B, boundary_indices, constrValues, h_fem);
	std::cout << h_fem.rows() << ", " << h_fem.cols() << std::endl;

	solveConstrainedSymmetric(A_dec, B, boundary_indices, constrValues, h_dec);
	std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;
}


void solveDirichletProblem(CGALTriangulation<Kernel>& tri, Eigen::MatrixXd& h_fem, Eigen::MatrixXd& h_dec)
{
    const int cntr = tri.centerVertex();
    const int n = tri.mesh.number_of_vertices();

    // Construct A 
    Eigen::SparseMatrix<double> A_fem, A_dec, L_fem, L_dec, M;
    //tri.massMatrix(M);
    tri.FEMLaplacian(L_fem);

	std::cout << "Using Mixed DEC" << std::endl;
    tri.DECLaplacianMixed(L_dec, &M);
    // tri.DECLaplacian(L_dec, &M);
    A_fem = L_fem;
	A_dec =-L_dec; 

    // solve the constrained problems
    std::vector<int> boundary_indices = tri.surfaceVertices(); 
    boundary_indices.push_back(cntr);

    Eigen::MatrixXd B(n, 1); B.setZero();
    Eigen::MatrixXd constrValues(boundary_indices.size(), 1);
    constrValues.setZero();
    constrValues(boundary_indices.size()-1, 0) = 1;

	/*
	std::cout << "SHAPES " <<std::endl;
	std::cout << constrValues.size() << std::endl;
	std::cout << B.rows() << ", " << B.cols()  << std::endl;
	std::cout << "FEM: " << std::endl;
	std::cout << A_fem.rows() << ", " << A_fem.cols()  << std::endl;
	std::cout << "DEC: " << std::endl;
	std::cout << A_dec.rows() << ", " << A_dec.cols()  << std::endl;
	*/

    solveConstrainedSymmetric(A_fem, B, boundary_indices, constrValues, h_fem);
	std::cout << h_fem.rows() << ", " << h_fem.cols() << std::endl;

    solveConstrainedSymmetric(A_dec, B, boundary_indices, constrValues, h_dec);
	std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;
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
  if (key == '2')
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

int 
loadMeshWithOrbitpoints(CGALTriangulation<Kernel> &tri, std::vector<int> &orbitinds, std::string filepath)
{
	// load Triangulation from file and read orbitpoints
	tri.read(filepath);
	std::cout << tri.mesh.number_of_vertices() << std::endl;

	typedef CGALTriangulation<Kernel>::Vertex_handle Vertex_handle;
	typedef CGALTriangulation<Kernel>::Point Point;

	orbitinds.clear();

	Point origin(0.,0.,0.);
	int originind = -1;
	
	for (Vertex_handle vh: tri.mesh.finite_vertex_handles()) {
		double dist = sqrt(CGAL::squared_distance(vh->point(), origin));		
		//std::cout << std::endl << vh->point() << std::endl;
		//std::cout << dist << std::endl;

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

Eigen::MatrixXd  normalizeHeatValues(Eigen::MatrixXd h) {

	Eigen::MatrixXd h_normed;
	h_normed.resize(h.rows(), 1);
	double mi = h.minCoeff();
	double ma = h.maxCoeff();
	for(int i = 0; i < h.size(); ++i)
		h_normed(i) = (h(i) - mi) / (ma - mi);
	return h_normed;
}

int main(int argc, char *argv[])
{

	if (argc < 2) {
		std::cout << "usage: argv[0] run_folder variance r_postfix (1 for gui output)" << std::endl;
		return 0;
	}

	double regnoise;
	std::string run_postfix = "";
	// no gui output
	bool silent = true;

	// no calculations, only show mesh
	bool viewer_only=false;
	std::string run_folder = argv[1];

	if (argc >= 3) {
		regnoise = std::stod(argv[2]);
		if (regnoise > 1e-10) {
			std::cout << "Regnoise specified" << std::endl;
			run_postfix += std::to_string(regnoise) + std::string("_"); 
		}
	}

	if (argc >= 4){
	run_postfix += argv[3] + std::string("_"); 
	//run_postfix += "_";	
	}

	if (argc >= 5) {
		if (atoi(argv[4])) silent = false;
	}
	bool meshwrite_only = false;
	if (argc >= 6) {
		if (atoi(argv[5])) meshwrite_only = true;
	}

	CGALTriangulation<Kernel> tri;
	int originind;
	std::vector<int> orbitinds;

	std::string meshNamesFile= run_folder + "meshes.txt";
	std::vector<std::string> meshNames;

	std::ifstream mNfile;
	mNfile.open(meshNamesFile); 
	std::string line;
	while(std::getline(mNfile, line)) meshNames.push_back(line);
	mNfile.close();

	int run_idx = 0;
	for(std::string run_name : meshNames) {
		std::cout << "processing run " << run_idx+1 << "/" << meshNames.size() << ": "  << run_name << std::endl;

		std::string filepath = run_folder + run_name + ".meshfile";

		if(loadMeshWithOrbitpoints(tri, orbitinds, filepath)) {
			originind = orbitinds.back();
			orbitinds.pop_back();
			std::cout << "Origin ind: " << originind << std::endl;
			std::cout << "Orbitinds: " << std::endl;
			for(int i: orbitinds) std::cout << i << " ";
			std::cout << std::endl;
		
		} else {
			std::cout << "Something went wrong loading the mesh" << std::endl;	
		}

		if (regnoise > 0) {
			// #################################
			// Replace by regular triangulation
			// #################################
			
			// set noise ifo mels
			double mels = tri.meanEdgeLengthSquared();
			double variance = mels * regnoise;  // (mels*mels*regnoise) * (mels*mels*regnoise);

			std::cout << "Old orbitinds size: " << orbitinds.size() << std::endl;
			std::cout << "Old orbitinds: ";
			for (int i : orbitinds ) std::cout << i << " ";
			std::cout << std::endl;
			tri.replaceMeshByRegular(variance, orbitinds, originind);
			std::cout << "New orbitinds size: " << orbitinds.size() << std::endl;
			std::cout << "New orbitinds: ";
			for (int i : orbitinds ) std::cout << i << " ";
			std::cout << std::endl;

			std::cout << "-------------------------------------" << std::endl;
			std::cout << run_name << std::endl;
			std::cout << "MELS: " << tri.meanEdgeLengthSquared() << std::endl;
			std::cout << "-------------------------------------" << std::endl;

			if (meshwrite_only) {
				std::string outfile = run_folder + run_name  + run_postfix + ".meshfile";
				tri.write(outfile);
				continue;
			}
		}

		// #########################################
		std::cout << "METRICS"  << std::endl;
		// #########################################
		enum Metric {minangle=0, amips, volume};
		Eigen::MatrixXd facecolors, x;
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

		std::cout << "clcd" << std::endl;

		bool normalize=false;
		double amips_max = 100;
		if (!normalize) {
			//normalize minangle by 70.5 deg
			for (int i=0; i < cell_metrics[minangle].size(); i++) cell_metrics[minangle][i] = cell_metrics[minangle][i] / 70.5;
			// normalize amips using the heuristically chosen max val amips_max and the min value 3
			for (int i=0; i < cell_metrics[amips].size(); i++) cell_metrics[amips][i] = (cell_metrics[amips][i] - 3) / amips_max;
		} 

		Eigen::MatrixXd cellcolors_volume; 
		igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[volume], true, cellcolors_volume);
		cellcolors[volume] = cellcolors_volume;
		Eigen::MatrixXd cellcolors_minangle; 
		igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[minangle], normalize, cellcolors_minangle);
		cellcolors[minangle] = cellcolors_minangle;
		Eigen::MatrixXd cellcolors_amips; 
		igl::colormap(igl::COLOR_MAP_TYPE_VIRIDIS, cell_metrics[amips], normalize, cellcolors_amips);
		cellcolors[amips] = cellcolors_amips;

		std::cout << "Viewer stuff " << std::endl;

		/* ################## HEAT DIFFUSION ################ */

		if (orbitinds.size() < 1) {
			std::cout << "No orbit points in file, abort" << std::endl;	
			return 0;
		}

		Eigen::MatrixXd h_fem, h_dec, h;
		solveHeatProblem(tri, h_fem, h_dec);
		//solveDirichletProblem(tri, h_fem, h_dec);

		std::cout << "...write heat vals to file... " << std::endl;

		std::string res_out_path = run_folder + run_name + run_postfix + "heatvals.csv";
		std::ofstream feil;
		feil.open(res_out_path);
		feil << "origin index"  << std::endl;
		feil << originind       << std::endl;
		feil << "orbit indices" << std::endl;
		for(int i=0; i < orbitinds.size() - 1; i++) feil << orbitinds[i] << ", ";
		feil << orbitinds[orbitinds.size()-1] << std::endl;
		feil << "h_fem" << ", " << "h_dec" << std::endl;

		for (int r = 0; r < h_fem.rows(); ++r) {
			feil << h_fem(r) << ", " << h_dec(r) << std::endl;
		}
		feil.close();

		std::cout << "Finished the feil" << std::endl;

		// normalize the heat values:
		x = normalizeHeatValues(h_dec);

		if (!silent) {
			/* ################## SHOW  MESH ####################*/


			Eigen::MatrixXd V;
			Eigen::MatrixXi F;

			double offset = 0.1;

			std::array<double, 4> plane{1,0,0,offset};
			std::vector<std::vector<int>> cmreturn = tri.cutMesh(plane, V, F);
			std::vector<int>  ids     = cmreturn[0];
			faceids = cmreturn[1];
			facecolors.resize(faceids.size(), 3);
			for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);

			/*
			Eigen::VectorXd x2();
			for(int i = 0; i < ids.size(); ++i)
				x2(i) = x(ids[i]);


			Eigen::MatrixXd color;
			igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
			viewer.data().set_colors(color);
			*/


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
			static bool showHeat   = false;

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

							if (!showHeat) {
								facecolors.resize(faceids.size(), 3);
								for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
								viewer.data().set_colors(facecolors);
							} else {

								Eigen::VectorXd x2(ids.size());
								for(int i = 0; i < ids.size(); ++i)
									x2(i) = x(ids[i]);


								Eigen::MatrixXd color;
								igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
								viewer.data().set_colors(color);
							
							}

							/*
							viewer.data().clear();
							viewer.data().set_mesh(V, F);
							viewer.data().uniform_colors(ambient, diffuse, specular);
							viewer.data().show_texture = 1;
							
							Eigen::VectorXd x2(ids.size());
							for(int i = 0; i < ids.size(); ++i)
								x2(i) = x(ids[i]);
							
							
							setTexture("../data/tex.png", viewer.data(cutMeshId));
							Eigen::MatrixXd UV(x2.rows(), 2);
							UV.col(0) = UV.col(1) = x2;
							viewer.data(cutMeshId).set_uv(UV);
							*/


						}
					}
				}
			};

			// Plot the mesh
			viewer.data().set_mesh(V, F);
			viewer.data().set_colors(facecolors);

			viewer.core().background_color.setOnes();
			viewer.callback_key_down = &key_down;

			/*
			if (!showHeat) {
				viewer.data().set_colors(facecolors);
			} else {

			Eigen::VectorXd x2(ids.size());
			for(int i = 0; i < ids.size(); ++i)
				x2(i) = x(ids[i]);


			Eigen::MatrixXd color;
			igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
			viewer.data().set_colors(color);
			
			}
			*/
			viewer.launch();
		}
	run_idx++;
	}
}


