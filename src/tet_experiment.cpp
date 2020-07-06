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
#include <cmath>

#include <unsupported/Eigen/SparseExtra>
#include <unsupported/Eigen/ArpackSupport>
#include "lodepng.h"
#include "CGALTriangulation.hpp"
#include "CGALPolyhedron.hpp"
#include "TetgenMeshPolyhedron.hpp"
#include "CGALMeshPolyhedron.hpp"
#include "SolveConstrained.hpp"

#include <unsupported/Eigen/SparseExtra>

#include "tetgen.h"
#include "IndexedTetMesh.hpp"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double> Kernel;
//typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel; //EPECKernel;
typedef Kernel::Point_3 Point;


#include <Eigen/Dense>

void retrieveShellIndices(CGALTriangulation<Kernel> tri, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell)
{
	/// retrieve the indices of all points lying on the three spherical shells by their distance to the origin

	typedef CGALTriangulation<Kernel>::Point Point;
	innerShell.clear();
	middleShell.clear();
	outerShell.clear();
	double d_i = 1.0;
	double d_m = 1.5;
	double d_o = 2.0;
	
	// this is the smallest threshold that seems to find all points
	double eps=1e-5;
	
	for (auto vh: tri.mesh.finite_vertex_handles()) {
		double dist = sqrt(CGAL::squared_distance(vh->point(), Point(CGAL::ORIGIN)));		
		if (dist < d_i + eps) {
			// inner shell, includes all points in the hole (if filled by our reg tri insertion)
			innerShell.push_back(vh->info());
		}
		if (fabs(dist - d_m) < eps) {
			middleShell.push_back(vh->info());	
		}
		if (fabs(dist - d_o) < eps) {
			outerShell.push_back(vh->info());	
		}
	}

}

void 
loadMeshWithShellIndices(CGALTriangulation<Kernel> &tri, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell, std::string filepath)
{
	// load Triangulation from file and read orbitpoints
	tri.read(filepath);
	retrieveShellIndices(tri, innerShell, middleShell, outerShell);
}

void solveHeatProblem(CGALTriangulation<Kernel>& tri, CGALTriangulation<Kernel>::Regular* reg, std::vector<int> innerShell, std::vector<int> outerShell, Eigen::MatrixXd& h_fem, Eigen::MatrixXd& h_dec, Eigen::MatrixXd& h_decreg)
{
	//const int cntr = tri.centerVertex();
	const int n = tri.mesh.number_of_vertices();

	// Construct A 
	Eigen::SparseMatrix<double> A_fem, A_dec, A_decreg, L_fem, L_dec, L_r, M, M_r;

	if (reg) {
		std::cout << "Regular found, calc the declaplaceregular" << std::endl;	
		tri.DECLaplacianRegular(*reg, L_r, &M_r);
	}

	//tri.massMatrix(M);
	tri.FEMLaplacian(L_fem);

	// std::cout << "Using Mixed DEC" << std::endl;
	// tri.DECLaplacianMixed(L_dec, &M);
	tri.DECLaplacian(L_dec, &M);
	const double t = tri.meanEdgeLengthSquared();
	A_fem = M + t * L_fem;
	A_dec = M - t * L_dec; 

	// solve the constrained problems
	std::vector<int> constrIndices(outerShell);
	constrIndices.insert(constrIndices.end(), innerShell.begin(), innerShell.end());

	//constrIndices.push_back(outerShell);
	//constrIndices.push_back(innerShell);

	Eigen::MatrixXd B(n, 1); B.setZero();
	Eigen::MatrixXd constrValues(constrIndices.size(), 1);
	constrValues.setZero();
	// set inner shell values to 1
	constrValues.block(outerShell.size(), 0, innerShell.size(), 1) = Eigen::MatrixXd::Ones(innerShell.size(),1);

	std::cout << "SHAPES " <<std::endl;
	std::cout << constrValues.size() << std::endl;
	std::cout << B.rows() << ", " << B.cols()  << std::endl;
	std::cout << "FEM: " << std::endl;
	std::cout << A_fem.rows() << ", " << A_fem.cols()  << std::endl;
	std::cout << "DEC: " << std::endl;
	std::cout << A_dec.rows() << ", " << A_dec.cols()  << std::endl;

	solveConstrainedSymmetric(A_fem, B, constrIndices, constrValues, h_fem);
	std::cout << h_fem.rows() << ", " << h_fem.cols() << std::endl;

	solveConstrainedSymmetric(A_dec, B, constrIndices, constrValues, h_dec);
	std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;

	if (reg) {	
		
		std::cout << "------------------- REG L CHECK ------------------" << std::endl;
		std::cout << "(L_dec - L_r).squaredNorm(): ";
		std::cout << (L_dec - L_r).squaredNorm() << std::endl;
		std::cout << "(L_dec - L_fem).squaredNorm(): ";
		std::cout << (L_dec - L_fem).squaredNorm() << std::endl;
		std::cout << "------------------- /REG L CHECK -----------------" << std::endl;


		A_decreg = M_r - t * L_r;
		std::cout << "DECreg: " << std::endl;
		std::cout << A_decreg.rows() << ", " << A_decreg.cols()  << std::endl;
		solveConstrainedSymmetric(A_decreg, B, constrIndices, constrValues, h_decreg);
		std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;
	}
}

void checkLP(CGALTriangulation<Kernel>& tri, Eigen::SparseMatrix<double> L, std::vector<int> innerShell) {
    const int nv = tri.mesh.number_of_vertices();
	// check linear precision
	Eigen::MatrixXd V(nv, 3);
	for(auto h : tri.mesh.finite_vertex_handles())
	{
		V(h->info(), 0) = h->point().x();
		V(h->info(), 1) = h->point().y();
		V(h->info(), 2) = h->point().z();
	}
	
	Eigen::MatrixXd LV = L * V;
	for(int i : tri.surfaceVerticesSlow())
	{
		LV.row(i).setZero();
	}
	/*
	for(int i : innerShell) 
	{
		LV.row(i).setZero();
	}
	*/

	std::cout << "linear precision " << LV.norm() << std::endl;
}

void solveDirichletProblem(CGALTriangulation<Kernel>& tri, CGALTriangulation<Kernel>::Regular* reg, std::vector<int> innerShell, std::vector<int> middleShell, std::vector<int> outerShell, Eigen::MatrixXd& h_fem, Eigen::MatrixXd& h_dec, Eigen::MatrixXd& h_opt, Eigen::MatrixXd& h_decreg, std::string laplaceSavePath)
{
	//const int cntr = tri.centerVertex();
	const int n = tri.mesh.number_of_vertices();

	// Construct A 
	Eigen::SparseMatrix<double> A_fem, A_dec, A_opt, A_decreg, L_fem, L_dec, L_opt, L_r, M, M_r;

	if (reg) {
		std::cout << "Regular found, calc the declaplaceregular" << std::endl;	
		tri.DECLaplacianRegular(*reg, L_r, &M_r);
	}

	//tri.massMatrix(M);
	tri.FEMLaplacian(L_fem);

	// std::cout << "Using Mixed DEC" << std::endl;
	// tri.DECLaplacianMixed(L_dec, &M);
	tri.DECLaplacian(L_dec, &M);

	// optimize laplacian
	L_opt = L_dec;
	float stepsize = 1;
	int maxits     = 10;
	int targetstyle = 1;
	//std::vector<int> emptylist;
	tri.DECLaplacianOptimized(L_opt, stepsize, maxits, targetstyle, innerShell);

	std::cout << "(L_dec - L_opt).norm() = " << (L_dec - L_opt).norm() << std::endl;

	if (!laplaceSavePath.empty()) {
		Eigen::saveMarket( L_dec, laplaceSavePath + "_Ldec.mtx");
		Eigen::saveMarket( L_fem, laplaceSavePath + "_Lfem.mtx");
		Eigen::saveMarket( L_opt, laplaceSavePath + "_Lopt.mtx");
	}


	const double t = tri.meanEdgeLengthSquared();
	A_fem =   L_fem;
	A_dec = - L_dec; 
	A_opt = - L_opt; 

	// solve the constrained problems
	std::vector<int> constrIndices(outerShell);
	constrIndices.insert(constrIndices.end(), innerShell.begin(), innerShell.end());

	//constrIndices.push_back(outerShell);
	//constrIndices.push_back(innerShell);

	Eigen::MatrixXd B(n, 1); B.setZero();
	Eigen::MatrixXd constrValues(constrIndices.size(), 1);
	constrValues.setZero();
	// set inner shell values to 1
	constrValues.block(outerShell.size(), 0, innerShell.size(), 1) = Eigen::MatrixXd::Ones(innerShell.size(),1);

	std::cout << "SHAPES " <<std::endl;
	std::cout << constrValues.size() << std::endl;
	std::cout << B.rows() << ", " << B.cols()  << std::endl;
	std::cout << "FEM: " << std::endl;
	std::cout << A_fem.rows() << ", " << A_fem.cols()  << std::endl;
	std::cout << "DEC: " << std::endl;
	std::cout << A_dec.rows() << ", " << A_dec.cols()  << std::endl;

	solveConstrainedSymmetric(A_fem, B, constrIndices, constrValues, h_fem);
	std::cout << h_fem.rows() << ", " << h_fem.cols() << std::endl;

	solveConstrainedSymmetric(A_dec, B, constrIndices, constrValues, h_dec);
	std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;

	solveConstrainedSymmetric(A_opt, B, constrIndices, constrValues, h_opt);
	std::cout << h_opt.size() << ", " << h_opt.cols() << std::endl;

	bool dbg = true;
	if (dbg) {
		std::cout << "--------- LP TEST ---------" << std::endl;
		std::vector<int> ignoreIndices;
		std::cout << " DEC: " << std::endl;
		checkLP(tri, L_dec, ignoreIndices);
		std::cout << " OPT: " << std::endl;
		checkLP(tri, L_opt, ignoreIndices);
		std::cout << "--------- /LP TEST ---------" << std::endl;
	}


	if (reg) {	
		
		std::cout << "------------------- REG L CHECK ------------------" << std::endl;
		std::cout << "(L_dec - L_r).squaredNorm(): ";
		std::cout << (L_dec - L_r).squaredNorm() << std::endl;
		std::cout << "(L_dec - L_fem).squaredNorm(): ";
		std::cout << (L_dec - L_fem).squaredNorm() << std::endl;
		/*
		std::cout << "Samples: L_dec" << std::endl;
		std::cout << L_dec.block(0,0,100,100) << std::endl;
		std::cout << "Samples: L_r" << std::endl;
		std::cout << L_r.block(0,0,100,100) << std::endl;
		*/

		std::cout << "------------------- /REG L CHECK -----------------" << std::endl;


		A_decreg = - L_r;
		std::cout << "DECreg: " << std::endl;
		std::cout << A_decreg.rows() << ", " << A_decreg.cols()  << std::endl;
		solveConstrainedSymmetric(A_decreg, B, constrIndices, constrValues, h_decreg);
		std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;
	}
}

void testOptimizedLaplace(CGALTriangulation<Kernel>& tri, double stepsize, int maxits, std::vector<int> ignoreIndices){
	std::cout << "Test opt Laplace" << std::endl;
	Eigen::SparseMatrix<double> L_dec, L_optimized, M;
	tri.DECLaplacian(L_dec, &M);
	L_optimized = L_dec;
	std::cout << "run optlaplace" << std::endl;
	tri.DECLaplacianOptimized(L_optimized, stepsize, maxits, 0, ignoreIndices);
}

void testSurfaceVertices(CGALTriangulation<Kernel> tri){
	std::vector<int> surfA, surfB;
	for(auto h : tri.mesh.finite_vertex_handles()) {
		std::vector<CGALTriangulation<Kernel>::Vertex_handle> adj;
		tri.mesh.adjacent_vertices(h, std::back_inserter(adj));
		bool adjtoinf=false;
		for (auto h: adj) {
			if (h->info() == -1) {
				adjtoinf=true;
			}	
		}
		if (adjtoinf) {
			surfA.push_back(h->info());
		}
	}

	surfB = tri.surfaceVertices();
	std::cout << "A.size " << surfA.size() << std::endl;
	std::cout << "B.size " << surfB.size() << std::endl;
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


Eigen::MatrixXd  normalizeHeatValues(Eigen::MatrixXd h) {

	Eigen::MatrixXd h_normed;
	h_normed.resize(h.rows(), 1);

	//lognorm=false;
	
	double mi = h.minCoeff();
	double ma = h.maxCoeff();
	//ma /= 10.;
	for(int i = 0; i < h.size(); ++i)
		h_normed(i) = (h(i) - mi) / (ma - mi);

	/*
	mi = h.minCoeff();
	ma = h.maxCoeff();
	//ma /= 10.;
	for(int i = 0; i < h.size(); ++i)
		h_normed(i) = (h_normed(i) - mi) / (ma - mi);
	*/

	return h_normed;
}

int main(int argc, char *argv[])
{

	// -------- ARG HANDLING ---------
	if (argc < 2) {
		std::cout << "usage: argv[0] run_folder variance r_postfix (1 for gui output)" << std::endl;
		return 0;
	}

	double regnoise;
	std::string run_postfix = "";
	// no gui output
	bool silent = false;

	// no calculations, only show mesh
	bool viewer_only=false;
	std::string run_folder = argv[1];

	//double stepsize			   = std::stod(argv[2]);
	//int    maxits              = atoi(argv[3]);

	bool meshwrite = false;
	bool saveLaplacians = false;
	bool output_mel = false;
	bool meshwrite_only = false;
	regnoise = -1;
	if (argc >= 3) {
		regnoise = std::stod(argv[2]);
		if (regnoise >= 0) {
			std::cout << "Regnoise specified" << std::endl;
		}
		run_postfix += std::to_string(regnoise) + std::string("_"); 
	}

	if (argc >= 4){
	run_postfix += argv[3] + std::string("_"); 
	//run_postfix += "_";	
	}

	if (argc >= 5) {
		if (atoi(argv[4])) silent = false;
	}
	if (argc >= 6) {
		if (atoi(argv[5])) meshwrite = true;
	}
	if (argc >= 7) {
		if (atoi(argv[6])) saveLaplacians = true;
	}
	if (argc >= 8) {
		if (atoi(argv[7])) output_mel = true;
	}


	/*
	if (argc >= 6) {
		if (atoi(argv[5])) meshwrite_only = true;
	}
	*/
	
	// --------  END ARG HANDLING ---------

	CGALTriangulation<Kernel> tri;

	// regular triangulation used for the DECLaplacianRegular and stored externally. 
	CGALTriangulation<Kernel>::Regular *reg=nullptr;
	CGALTriangulation<Kernel>::Regular regtri;

	std::vector<int> innerShell;
	std::vector<int> middleShell;
	std::vector<int> outerShell;

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

		loadMeshWithShellIndices(tri, innerShell, middleShell, outerShell, filepath);

		double mels = tri.meanEdgeLengthSquared();
		if (regnoise >= 0) {
			// #################################
			// Replace by regular triangulation
			// #################################
			
			// set noise ifo mels
			double variance = mels * regnoise;  // (mels*mels*regnoise) * (mels*mels*regnoise);

			std::cout << "Number of Cells pre: " << tri.mesh.number_of_finite_cells() << std::endl;


			// generate the regular triangulation
			regtri = tri.generateRandomRegular(variance);
			reg    = &regtri;

			std::cout << "Number of Cells reg: " << reg->number_of_finite_cells() << std::endl;

			// list with shellIndices that will be updated
			/*
			std::vector<std::vector<int>> shellIndices;
			shellIndices.push_back(innerShell);
			shellIndices.push_back(middleShell);
			shellIndices.push_back(outerShell);
			*/

			// replace the triangulation by the regular one
			tri.replaceMeshByRegular(*reg, innerShell, middleShell, outerShell);

			std::cout << "Middle Shell pre: " << std::endl;
			for (int i =0; i < 100; i++) {
			std::cout << middleShell[i] << ", ";
			}
			std::cout << std::endl;

			// update the shellIndices
			/*
			innerShell  = shellIndices[0];
			middleShell = shellIndices[1];
			outerShell  = shellIndices[2];
			*/

			std::cout << "Middle Shell post: " << std::endl;
			for (int i =0; i < 100; i++) {
			std::cout << middleShell[i] << ", ";
			}
			std::cout << std::endl;


			std::cout << "Number of Cells post: " << tri.mesh.number_of_finite_cells() << std::endl;

			mels = tri.meanEdgeLengthSquared(); 

			std::cout << "-------------------------------------" << std::endl;
			std::cout << run_name << std::endl;
			std::cout << "MELS: " << mels << std::endl;
			std::cout << "-------------------------------------" << std::endl;

			if (meshwrite) {
				std::string outfile = run_folder + run_name  + run_postfix + ".meshfile";
				tri.write(outfile);
			}
			/*
			if (meshwrite_only) {
				std::string outfile = run_folder + run_name  + run_postfix + ".meshfile";
				tri.write(outfile);
				continue;
			}
			*/
		}

		if (output_mel) {
			std::ofstream melFile;
			std::string melFileName =run_folder + "melvalues.txt";
			melFile.open(melFileName, std::ios_base::app); // append instead of overwrite
			std::string prfx = run_name + run_postfix;
			melFile << prfx << ", " << sqrt(mels)  << std::endl; 
			melFile.close();
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

		if (regnoise >= 0) {
			// write values to file again (reg might have changed them)
			std::ofstream feil;
			std::string metrics_out_path = run_folder + run_name + run_postfix + "metrics.csv";
			feil.open(metrics_out_path);
			feil << metric_names[minangle] << "," << metric_names[amips] << "," << metric_names[volume] << std::endl;
			for(int i; i < cell_metrics[volume].size(); i++) {
				feil << cell_metrics[minangle](i) << "," << cell_metrics[amips](i) << "," << cell_metrics[volume](i) << std::endl;
			}
			feil.close();
		}

		std::cout << "Viewer stuff " << std::endl;

		/* ################## HEAT DIFFUSION ################ */

		if (innerShell.size() < 1 || middleShell.size() < 1 || outerShell.size() < 1) {
			std::cout << "Error with the shell points, abort." << std::endl;	
			std::cout << "innerShell.size: " << innerShell.size() << std::endl;	
			std::cout << "middleShell.size: " << middleShell.size() << std::endl;	
			std::cout << "outerShell.size: " << outerShell.size() << std::endl;	
			return 0;
		}

		testSurfaceVertices(tri);

		Eigen::MatrixXd h_fem, h_dec, h_opt, h_decreg;
		//solveHeatProblem(tri, reg, innerShell, outerShell, h_fem, h_dec, h_decreg);
		std::string laplacesavepath = "";
		if (saveLaplacians) {
			laplacesavepath = run_folder + run_name  + run_postfix;
		}

		// TEST OPTIMIZED LAPLACE
		//testOptimizedLaplace(tri, stepsize, maxits);
		//return 0;

		solveDirichletProblem(tri, reg, innerShell, middleShell, outerShell, h_fem, h_dec, h_opt, h_decreg, laplacesavepath);

		std::cout << "...write heat vals to file... " << std::endl;
		std::string res_out_path = run_folder + run_name + run_postfix + "heatvals.csv";
		std::ofstream feil;
		feil.open(res_out_path);
		feil << "middle shell indices" << std::endl;
		for(int i=0; i < middleShell.size() - 1; i++) feil << middleShell[i] << ", ";
		feil << middleShell[middleShell.size()-1] << std::endl;
		
		if (reg) {
			// include values of the reg laplacian if caluclated
			feil << "h_fem" << "," << "h_dec" << "," << "h_opt"<< "," << "h_decreg" << std::endl;
			for (int r = 0; r < h_fem.rows(); ++r) {
				feil << h_fem(r) << "," << h_dec(r) << "," << h_opt(r) << "," << h_decreg(r) << std::endl;
			}
		} else {
			// include values of the reg laplacian if caluclated
			feil << "h_fem" << "," << "h_dec" << "," << "h_opt"<< std::endl;
			for (int r = 0; r < h_fem.rows(); ++r) {
				feil << h_fem(r) << "," << h_dec(r) << "," << h_opt(r)<< std::endl;
			}
		
		}

		feil.close();
		std::cout << "Finished the feil" << std::endl;

		// normalize the heat values:
		//x = normalizeHeatValues(h_dec);
		x = normalizeHeatValues(h_fem);

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
			static bool showHeat   = true;

			// Add content to the default menu window
			menu.callback_draw_viewer_menu = [&]()
			{
				// Draw parent menu content
				
				if (ImGui::CollapsingHeader("Presentation", ImGuiTreeNodeFlags_DefaultOpen))
				{
					Metric oldmetric = metric_shown;
					ImGui::Combo("Metric", (int *)(&metric_shown), "MinAngle\0AMIPS\0Volume\0\0");

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


