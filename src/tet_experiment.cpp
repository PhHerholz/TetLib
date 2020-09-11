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

void retrieveShellIndices(CGALTriangulation<Kernel> tri, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell, int &originind, double radius=2., double eps=1e-4)
{
	/// retrieve the indices of all points lying on the three spherical shells by their distance to the origin

	typedef CGALTriangulation<Kernel>::Point Point;
	innerShell.clear();
	middleShell.clear();
	outerShell.clear();
	/*
	double d_i = 1.0;
	double d_m = 1.5;
	double d_o = 2.0;
	*/
	double d_i = radius/2.;
	double d_m = radius*3/4;
	double d_o = radius;
	
	// this is the smallest threshold that seems to find all points
	// (for doublesphere it seemed to be 1e-5, for single 1e-4
	
	for (auto vh: tri.mesh.finite_vertex_handles()) {
		double dist = sqrt(CGAL::squared_distance(vh->point(), Point(CGAL::ORIGIN)));		

		if (dist < eps) {
			originind = vh->info();	
			std::cout << "NEW ORIGIN: " << vh->point() << std::endl;
			std::cout << "Origin ind: " << originind << std::endl;
		}

		if (fabs(dist - d_i) < eps) {
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
loadMeshWithShellIndices(CGALTriangulation<Kernel> &tri, std::vector<int> &innerShell, std::vector<int> &middleShell, std::vector<int> &outerShell, std::string filepath, int &originind, double outerrad, double eps)
{
	// load Triangulation from file and read orbitpoints
	tri.read(filepath);

	/*
	double maxdist=0.;
	for (auto vh: tri.mesh.finite_vertex_handles()) {
		double dist = sqrt(CGAL::squared_distance(vh->point(), Point(CGAL::ORIGIN)));		
		if (dist > maxdist) maxdist = dist;
	}
	std::cout << "MAXDIST: " << maxdist << std::endl;
	*/
	retrieveShellIndices(tri, innerShell, middleShell, outerShell, originind, outerrad, eps);
}


void 
calcAndWriteHeatGradsAllVertices(CGALTriangulation<Kernel> tri, Eigen::MatrixXd h_fem, Eigen::MatrixXd h_dec, Eigen::MatrixXd h_reg, bool includereg, std::string outpath)
{
	// calc heat gradients
	Eigen::MatrixXd heatGradField_fem, heatGradField_dec, heatGradField_reg;
	tri.calcHeatGradientAllVertices(h_fem, heatGradField_fem);
	tri.calcHeatGradientAllVertices(h_dec, heatGradField_dec);
	if (includereg) {
		tri.calcHeatGradientAllVertices(h_reg, heatGradField_reg);
	}

	int nv = tri.mesh.number_of_vertices();
	Eigen::MatrixXd V(nv, 3);
	for(auto h : tri.mesh.finite_vertex_handles())
	{
		V(h->info(), 0) = h->point().x();
		V(h->info(), 1) = h->point().y();
		V(h->info(), 2) = h->point().z();
	}


	std::ofstream feil;
	feil.open(outpath);

	// vertex header
	feil << "vx,vy,vz";
	// hg headers
	feil << "," << "hg_fem_x" << "," << "hg_fem_y" << "," << "hg_fem_z";
	feil << "," << "hg_dec_x" << "," << "hg_dec_y" << "," << "hg_dec_z";
	if (includereg) feil << "," << "hg_reg_x" << "," << "hg_reg_y" << "," << "hg_reg_z";
	feil << std::endl;

	for (int i=0; i<heatGradField_fem.rows();++i) {

		feil << V(i, 0) << "," << V(i,1) << "," << V(i,2);
		// hg values
		feil << "," << heatGradField_fem(i, 0) <<  "," << heatGradField_fem(i, 1) << "," << heatGradField_fem(i, 2);
		feil << "," << heatGradField_dec(i, 0) <<  "," << heatGradField_dec(i, 1) << "," << heatGradField_dec(i, 2);
		if (includereg) feil << "," << heatGradField_reg(i, 0) <<  "," << heatGradField_reg(i, 1) << "," << heatGradField_reg(i, 2);

		feil << std::endl;
	}
}

void 
calcAndWriteHeatGradsAndCentroidsAllCells(CGALTriangulation<Kernel> tri, Eigen::MatrixXd h_fem, Eigen::MatrixXd h_dec, Eigen::MatrixXd h_reg, bool includereg, std::string outpath)
{
	//calc cell centroids
	Eigen::MatrixXd heatGradField, cellCentroids;
	tri.calcCentroidAllCells(cellCentroids);

	// calc heat gradients
	Eigen::MatrixXd heatGradField_fem, heatGradField_dec, heatGradField_reg;
	tri.calcHeatGradientAllCells(h_fem, heatGradField_fem);
	tri.calcHeatGradientAllCells(h_dec, heatGradField_dec);
	if (includereg) {
		tri.calcHeatGradientAllCells(h_reg, heatGradField_reg);
	}

	std::ofstream feil;
	feil.open(outpath);
	// centroid header
	feil << "centroid_x" << "," << "centroid_y" << "," << "centroid_z";
	// hg headers
	feil << "," << "hg_fem_x" << "," << "hg_fem_y" << "," << "hg_fem_z";
	feil << "," << "hg_dec_x" << "," << "hg_dec_y" << "," << "hg_dec_z";
	if (includereg) feil << "," << "hg_reg_x" << "," << "hg_reg_y" << "," << "hg_reg_z";
	feil << std::endl;

	for (int i=0; i<cellCentroids.rows();++i) {
		// centroid values
		feil << cellCentroids(i, 0) <<  "," << cellCentroids(i, 1) << "," << cellCentroids(i, 2);

		// hg values
		feil << "," << heatGradField_fem(i, 0) <<  "," << heatGradField_fem(i, 1) << "," << heatGradField_fem(i, 2);
		feil << "," << heatGradField_dec(i, 0) <<  "," << heatGradField_dec(i, 1) << "," << heatGradField_dec(i, 2);
		if (includereg) feil << "," << heatGradField_reg(i, 0) <<  "," << heatGradField_reg(i, 1) << "," << heatGradField_reg(i, 2);

		feil << std::endl;
	}
}

void solveHeatProblem(CGALTriangulation<Kernel>& tri, CGALTriangulation<Kernel>::Regular* reg, std::vector<int> innerShell, std::vector<int> outerShell, Eigen::MatrixXd& h_fem, Eigen::MatrixXd& h_dec, Eigen::MatrixXd& h_decreg, Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& M_dec, Eigen::SparseMatrix<double>& M_decreg, std::string reglaplaceloadpath, double t=-1)
{
	//const int cntr = tri.centerVertex();
	const int n = tri.mesh.number_of_vertices();

	// Construct A 
	Eigen::SparseMatrix<double> A_fem, A_dec, A_decreg, L_fem, L_dec, L_r, Id(n,n);
	Id.setIdentity();


	bool useregular = (reg || !reglaplaceloadpath.empty());

	if (useregular) {
		if (!reglaplaceloadpath.empty()) {
			std::cout << "Regpath given, load decreg L and M from file" << std::endl;
			std::cout << "path: " << reglaplaceloadpath << std::endl;
			Eigen::loadMarket(L_r,      reglaplaceloadpath + "decregLoptimized.mtx");
			Eigen::loadMarket(M_decreg, reglaplaceloadpath + "decregMoptimized.mtx");

			//std::cout << "L_r.shape = " << L_r.rows() << "," << L_r.cols() << std::endl;
			//std::cout << "L_r.shape = " << M_decreg.rows() << "," << M_decreg.cols() << std::endl;
		} else {
			std::cout << "Regular found, calc the declaplaceregular" << std::endl;	
			tri.DECLaplacianRegular(*reg, L_r, &M_decreg);
		}
	}

	tri.massMatrix(M);
	tri.FEMLaplacian(L_fem);
	tri.DECLaplacian(L_dec, &M_dec);

	if (t<0) {
		std::cout << "t < 0 -> set to mean edge length squared" << std::endl;
		t = tri.meanEdgeLengthSquared();
	}
	std::cout << "t = " << t << std::endl;

	A_fem = M     + t * L_fem;
	A_dec = M_dec - t * L_dec; 

	Eigen::MatrixXd b_base(n, 1); b_base.setZero();
	for (int i: innerShell) {
		b_base(i) = 1.;
	}

	std::cout << "innerShell.size(): " << innerShell.size() << std::endl;
	std::cout << "start solve" << std::endl;

	bool enforce_border = false;
	// Solve FEM
	Eigen::MatrixXd B_fem = b_base; //M * b_base;
	if (enforce_border) {
		Eigen::MatrixXd constrValuesFEM(outerShell.size(), 1);
		constrValuesFEM.setZero();
		solveConstrainedSymmetric(A_fem, B_fem, outerShell, constrValuesFEM, h_fem);
	} else {
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol_fem;
		chol_fem.analyzePattern(A_fem);
		chol_fem.factorize(A_fem);
		h_fem = chol_fem.solve(B_fem);
	}

	// Solve DEC
	Eigen::MatrixXd B_dec = b_base; // M_dec * b_base;
	if (enforce_border) {
		std::cout << "...fixing border" << std::endl;
		std::cout << "outerShelll.size = " << outerShell.size() << std::endl;
		Eigen::MatrixXd constrValuesDEC(outerShell.size(), 1);
		constrValuesDEC.setZero();
		solveConstrainedSymmetric(A_dec, B_dec, outerShell, constrValuesDEC, h_dec);
	} else {
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol_dec;
		chol_dec.analyzePattern(A_dec);
		chol_dec.factorize(A_dec);
		h_dec = chol_dec.solve(B_dec);
	}

	if (useregular) {	
		A_decreg = M_decreg - t * L_r;
		Eigen::MatrixXd B_reg = b_base; // M_dec * b_base;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol_reg;
		chol_reg.analyzePattern(A_decreg);
		chol_reg.factorize(A_decreg);
		h_decreg = chol_reg.solve(B_reg);
	}
}

bool checkLP(CGALTriangulation<Kernel>& tri, Eigen::SparseMatrix<double> L, std::vector<int> ignoreIndices) {
    const int nv = tri.mesh.number_of_vertices();
	// check linear precision
	Eigen::MatrixXd V(nv, 3);
	for(auto h : tri.mesh.finite_vertex_handles())
	{
		V(h->info(), 0) = h->point().x();
		V(h->info(), 1) = h->point().y();
		V(h->info(), 2) = h->point().z();
	}

	if (ignoreIndices.size() == 0) {
		ignoreIndices = tri.surfaceVerticesSlow();
	}
	
	Eigen::MatrixXd LV = L * V;
	for(int i : ignoreIndices)
	{
		LV.row(i).setZero();
	}

	std::cout << "LV.sum() = " << LV.sum() << std::endl;
	if (LV.norm() > 1e-7) {
		std::cout << "Linear Precision Error !! " << std::endl;
		std::cout << "LV.norm= " << LV.norm() << std::endl;
		std::cout << "Problemrows: ";
		int cnt = 0;
		for (int i=0; i<LV.rows(); ++i) {
			double rowsum = LV.row(i).sum(); 
			if (rowsum> 1e-7) {
				std::cout << i << ": " << rowsum << ", ";	
				cnt++;
			}
		}
		std::cout << std::endl;
		std::cout << "Total of " << cnt << "problematic rows" << std::endl;
		//return false;
	}
	return true;
}

bool solveDirichletProblem(CGALTriangulation<Kernel>& tri, CGALTriangulation<Kernel>::Regular* reg, std::vector<int> innerShell, std::vector<int> middleShell, std::vector<int> outerShell, Eigen::MatrixXd& h_fem, Eigen::SparseMatrix<double>& M, Eigen::MatrixXd& h_dec,Eigen::SparseMatrix<double>& M_dec, Eigen::MatrixXd& h_opt, Eigen::MatrixXd& h_bndropt, Eigen::MatrixXd& h_decreg, Eigen::SparseMatrix<double>& M_decreg, std::string laplaceSavePath, int maxits, bool ignoreInner, bool addOpt, bool addBoundaryOpt, std::string traininglogfilebase, bool useMassMatrices, std::string reglaplaceloadpath)
{
	//const int cntr = tri.centerVertex();
	const int n = tri.mesh.number_of_vertices();

	// Construct A 
	Eigen::SparseMatrix<double> A_fem, A_dec, A_opt, A_bndropt, A_decreg, L_fem, L_dec, L_opt, L_bndropt, L_r; 

	bool useregular = (reg || !reglaplaceloadpath.empty());

	if (useregular) {
		if (!reglaplaceloadpath.empty()) {
			std::cout << "Regpath given, load decreg L and M from file" << std::endl;
			std::cout << "path: " << reglaplaceloadpath << std::endl;
			Eigen::loadMarket(L_r,      reglaplaceloadpath + "decregLoptimized.mtx");
			Eigen::loadMarket(M_decreg, reglaplaceloadpath + "decregMoptimized.mtx");

			//std::cout << "L_r.shape = " << L_r.rows() << "," << L_r.cols() << std::endl;
			//std::cout << "L_r.shape = " << M_decreg.rows() << "," << M_decreg.cols() << std::endl;
		} else {
			std::cout << "Regular found, calc the declaplaceregular" << std::endl;	
			tri.DECLaplacianRegular(*reg, L_r, &M_decreg);
		}
	}

	tri.massMatrix(M);
	tri.FEMLaplacian(L_fem);

	// std::cout << "Using Mixed DEC" << std::endl;
	// tri.DECLaplacianMixed(L_dec, &M);
	tri.DECLaplacian(L_dec, &M_dec);

	// optimize laplacian
	float stepsize = 0.1; // 0.0001;
	int targetstyle = 1; // sum of negative off-diagonal entries

	/*
	std::vector<int> ignoreIndices;
	ignoreIndices.insert(ignoreIndices.end(), outerShell.begin(), outerShell.end());
	if (ignoreInner) {
		ignoreIndices.insert(ignoreIndices.end(), innerShell.begin(), innerShell.end());
	} 
	*/
	std::vector<int> ignoreIndices = tri.surfaceVerticesSlow();

	if (addOpt) {
		L_opt = L_dec;
		tri.DECLaplacianOptimized(L_opt, stepsize, maxits, targetstyle, ignoreIndices, true, traininglogfilebase + "opt_trainlogs.csv");
	}

	if (addBoundaryOpt) {
		L_bndropt = L_dec;
		tri.DECLaplacianOptimized(L_bndropt, stepsize, maxits, targetstyle, ignoreIndices, false, traininglogfilebase +  "bndropt_trainlogs.csv");
	}

	if (!laplaceSavePath.empty()) {
		Eigen::saveMarket( L_dec, laplaceSavePath + "_Ldec.mtx");
		Eigen::saveMarket( L_fem, laplaceSavePath + "_Lfem.mtx");
		if (addOpt) {
			Eigen::saveMarket( L_opt, laplaceSavePath + "_Lopt.mtx");
		}
		if (useregular) {
			Eigen::saveMarket(L_r, laplaceSavePath + "_Lreg.mtx");
		}
		if (addBoundaryOpt) {
			Eigen::saveMarket(L_bndropt, laplaceSavePath + "_Lbndropt.mtx");
		}

		std::string res_out_path = laplaceSavePath + "_shellIndices.csv";
		std::ofstream feil;
		feil.open(res_out_path);

		feil << "inner shell indices" << std::endl;
		for(int i=0; i < innerShell.size() - 1; i++) feil << innerShell[i] << ", ";
		feil << innerShell[innerShell.size()-1] << std::endl;

		feil << "middle shell indices" << std::endl;
		for(int i=0; i < middleShell.size() - 1; i++) feil << middleShell[i] << ", ";
		feil << middleShell[middleShell.size()-1] << std::endl;

		feil << "outer shell indices" << std::endl;
		for(int i=0; i < outerShell.size() - 1; i++) feil << outerShell[i] << ", ";
		feil << outerShell[outerShell.size()-1] << std::endl;
		feil.close();
	}


	const double t = tri.meanEdgeLengthSquared();
	A_fem		=   L_fem;
	A_dec		= - L_dec; 
	A_opt		= - L_opt; 
	A_bndropt	= - L_bndropt; 

	// solve the constrained problems
	std::vector<int> constrIndices(outerShell);
	constrIndices.insert(constrIndices.end(), innerShell.begin(), innerShell.end());

	// CONSTR VALUES
	Eigen::MatrixXd B(n, 1); B.setZero();
	Eigen::MatrixXd constrValues(constrIndices.size(), 1);
	constrValues.setZero();
	// set inner shell values to 1
	constrValues.block(outerShell.size(), 0, innerShell.size(), 1) = Eigen::MatrixXd::Ones(innerShell.size(),1);

	double massscaling= 1.; //1e-6;//M.coeffs().minCoeff();
	Eigen::MatrixXd constrValuesFEM = constrValues;
	Eigen::MatrixXd constrValuesDEC = constrValues;
	if (useMassMatrices) {
		for (int i=0; i < innerShell.size(); ++i) {
			constrValuesFEM(outerShell.size() + i, 0) *= M.coeff(innerShell[i], innerShell[i])     / massscaling;
			constrValuesDEC(outerShell.size() + i, 0) *= M_dec.coeff(innerShell[i], innerShell[i]) / massscaling;
		}
	}

	std::cout << "SHAPES " <<std::endl;
	std::cout << constrValues.size() << std::endl;
	std::cout << B.rows() << ", " << B.cols()  << std::endl;
	std::cout << "FEM: " << std::endl;
	std::cout << A_fem.rows() << ", " << A_fem.cols()  << std::endl;
	std::cout << "DEC: " << std::endl;
	std::cout << A_dec.rows() << ", " << A_dec.cols()  << std::endl;

	std::cout << "n: " << n << std::endl;
	std::cout << "outerShell.size(): " << outerShell.size() << std::endl;
	std::cout << "innerShell.size(): " << innerShell.size() << std::endl;
	std::cout << "constrValues.size(): " << constrValues.size() << std::endl;

	std::cout << "...solve fem" << std::endl;
	solveConstrainedSymmetric(A_fem, B, constrIndices, constrValuesFEM, h_fem);
	std::cout << h_fem.rows() << ", " << h_fem.cols() << std::endl;

	std::cout << "...solve dec" << std::endl;
	solveConstrainedSymmetric(A_dec, B, constrIndices, constrValuesDEC, h_dec);
	std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;

	if (useMassMatrices) {
		for (int i = 0; i < h_dec.size(); ++i) {
			h_fem(i) /= M.coeff(i,i)     / massscaling;
			h_dec(i) /= M_dec.coeff(i,i) / massscaling;
		}	
	}


	if (addOpt) {
		std::cout << "...solve opt" << std::endl;
		solveConstrainedSymmetric(A_opt, B, constrIndices, constrValues, h_opt);
		std::cout << h_opt.size() << ", " << h_opt.cols() << std::endl;
	}

	if (addBoundaryOpt) {
		solveConstrainedSymmetric(A_bndropt, B, constrIndices, constrValues, h_bndropt);
	}

	bool dbg = true;
	if (dbg) {
		std::cout << "--------- Shell Index sizes: --------" << std::endl;
		std::cout << "Ignore Indices: " << ignoreIndices.size() << std::endl;
		std::cout << "innerShell    : " << innerShell.size() << std::endl;
		std::cout << "middleShell    : " << middleShell.size() << std::endl;
		std::cout << "outerShell    : " << outerShell.size() << std::endl;

		std::cout << "--------- LP TEST ---------" << std::endl;
		std::cout << " DEC: " << std::endl;
		if (!checkLP(tri, L_dec, ignoreIndices)) return false;
		std::cout << "...done" << std::endl;
		if (addOpt) {
			std::cout << " OPT: " << std::endl;
			if (!checkLP(tri, L_opt, ignoreIndices)) return false;
			std::cout << "...done" << std::endl;
		}
		if (useregular) {
			std::cout << " DECREG: " << std::endl;
			if (!checkLP(tri, L_r, ignoreIndices)) return false;
			std::cout << "...done" << std::endl;
		}
		std::cout << "--------- /LP TEST ---------" << std::endl;
	}


	if (useregular) {	
		
		std::cout << "------------------- REG L CHECK ------------------" << std::endl;
		std::cout << "(L_dec - L_r).squaredNorm(): ";
		std::cout << (L_dec - L_r).squaredNorm() << std::endl;
		std::cout << "(L_dec - L_fem).squaredNorm(): ";
		std::cout << (L_dec - L_fem).squaredNorm() << std::endl;
		std::cout << "------------------- /REG L CHECK -----------------" << std::endl;

		A_decreg = - L_r;
		std::cout << "DECreg: " << std::endl;
		std::cout << A_decreg.rows() << ", " << A_decreg.cols()  << std::endl;
		solveConstrainedSymmetric(A_decreg, B, constrIndices, constrValues, h_decreg);
		std::cout << h_dec.size() << ", " << h_dec.cols() << std::endl;
	}

	return true;
}

void testOptimizedLaplace(CGALTriangulation<Kernel>& tri, double stepsize, int maxits){
	std::cout << "Test opt Laplace" << std::endl;
	Eigen::SparseMatrix<double> L_dec, L_optimized, M;
	tri.DECLaplacian(L_dec, &M);
	L_optimized = L_dec;
	std::cout << "run optlaplace" << std::endl;
	std::vector<int> ignoreIndices;
	tri.DECLaplacianOptimized(L_optimized, stepsize, maxits, 0, ignoreIndices, "");
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

	bool heatdiffusion = true;
	bool dirichlet     = false;

	double regnoise;
	std::string run_postfix = "";
	// no gui output
	bool silent = true;

	std::string run_folder = argv[1];
	bool addOpt					  = false;
	bool addBoundaryOpt			  = false;
	bool includeMass              = true;

	//double stepsize			   = std::stod(argv[2]);
	//int    maxits              = atoi(argv[3]);
	static bool showHeat   = true;

	bool meshwrite      = false;
	bool regwrite       = false;
	bool saveLaplacians = true ;
	bool output_mel     = false;
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

	double boundaryCellMinVol = 0.;

	double heat_timestep= 0.005;
	if (argc >= 5) {
		heat_timestep = std::stod(argv[4]);
	}

	if (argc >= 6) {
		if (atoi(argv[5])) silent = false;
	}

	std::cout << "Silent: " << silent << std::endl;

	int maxits = 500;
	/*
	if (argc >= 5) {
		maxits = atoi(argv[4]);	
		run_postfix += std::string("MI") +  std::to_string(maxits) + std::string("_");
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

	*/

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

		// check the sphere size (important for loading the shell indices)
		double outerrad = 2.;
		double eps      = 1e-4;
		if (run_name.rfind("OriginSphere", 0) == 0){
			eps = 1e-2;
			outerrad = 1.;
		}

		int originind = -1;
		loadMeshWithShellIndices(tri, innerShell, middleShell, outerShell, filepath, originind, outerrad, eps);

		bool singleSphere = false;
		bool embeddedDoubleSphere = false;
		bool removeInner = true;
		if (run_name.rfind("SingleSphere", 0) == 0) {
			singleSphere = true;	
			removeInner  = false;

			middleShell = innerShell;
			innerShell.clear();
			innerShell.push_back(originind);

		} else if (run_name.rfind("EmbeddedDoubleSphere", 0) == 0) {
			embeddedDoubleSphere = true;
			removeInner = false;
		} else if (run_name.rfind("OriginSphere", 0) == 0) {
			singleSphere = true;
			removeInner = false;
			innerShell.clear();
			innerShell.push_back(originind);
		}
		std::cout << "..loaded Singlesphere : " << singleSphere << std::endl;
		std::cout << "Sphere Sizes: " << innerShell.size() << "," << middleShell.size() << "," << outerShell.size() << std::endl;

		double mels = tri.meanEdgeLengthSquared();
		if (regnoise >= 0) {
			// #################################
			// Replace by regular triangulation
			// #################################
			
			// set noise ifo mels
			double variance = mels * regnoise;  // (mels*mels*regnoise) * (mels*mels*regnoise);

			std::cout << "Number of Cells pre: " << tri.mesh.number_of_finite_cells() << std::endl;

			// replace the triangulation by the regular one

			int maxtries = 100;
			for (int tint=0; tint<maxtries; ++tint) { 
				// sometimes the noise process loses the origin (for the singlesphere)
				// in this case try maxtries times if another sample works

				// generate the regular triangulation
				regtri = tri.generateRandomRegular(variance);

				if(singleSphere) {
					// check that origin is contained
					bool containsorigin=false;
					for (auto vh: regtri.finite_vertex_handles()){
						if (vh->info() == innerShell[0]) {
							containsorigin=true;
						}
					} 
					if (containsorigin) break;
					std::cout << "... error with shell points, resample (try " << tint+1 << "/" << maxtries << ")" << std::endl;
				} else {
					break;
				}
			}

			reg    = &regtri;
			tri.replaceMeshByRegular(*reg, innerShell, middleShell, outerShell, boundaryCellMinVol, true, removeInner);
			if (innerShell.size() < 1 || middleShell.size() < 1 || outerShell.size() < 1) {
				std::cout << "spheresizes post: " << innerShell.size() << ", " << middleShell.size() << ", " << outerShell.size() << std::endl;
			}

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

			if (regwrite && reg) {
				std::string regweightsoutpath = run_folder + run_name  + run_postfix + "_regweights.csv";
				std::ofstream regweightsfile;
				regweightsfile.open(regweightsoutpath);
				regweightsfile << "regweight" << std::endl;
				int cntr = 0;
				for(auto it = reg->finite_vertices_begin(); it != reg->finite_vertices_end(); ++it)
				{
					if (it->info() != cntr) {
						regweightsfile << "ERROR: it->info() not in order! " << std::endl;
						break;
					}
					++cntr;
					regweightsfile << it->point().weight() << std::endl; 
				}
				regweightsfile.close();
			}

		} /*else if (regnoise == -2.) {
			// load reg tri with weights from file

			std::string regfilepath = run_folder + run_name + "_regweights.csv";
			regtri = tri.generateRegularFromWeightsfile(regfilepath);
			reg    = &regtri;
		
		}
		*/

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

		Eigen::VectorXd Dflags, CCCflags;
		tri.calcIsDelaunayFlagAllCells(Dflags);
		tri.calcContainsCircumcenterFlagAllCells(CCCflags);

		if (regnoise >= 0) {
			// write values to file again (reg might have changed them)
			std::ofstream feil;
			std::string metrics_out_path = run_folder + run_name + run_postfix + "metrics.csv";
			feil.open(metrics_out_path);
			feil << metric_names[minangle] << "," << metric_names[amips] << "," << metric_names[volume] << "," << "delaunayflag,containsccflag" << std::endl;
			for(int i; i < cell_metrics[volume].size(); i++) {
				feil << cell_metrics[minangle](i) << "," << cell_metrics[amips](i) << "," << cell_metrics[volume](i) <<
					"," << Dflags(i) << "," << CCCflags(i) << std::endl;
			}
			feil.close();
		}


		if (dirichlet) {
			std::cout << "DIRICHLET EXPERIMENT" << std::endl;

			/* ################## HEAT DIFFUSION ################ */
			if (innerShell.size() < 1 || middleShell.size() < 1 || outerShell.size() < 1) {
				std::cout << "Error with the shell points, abort." << std::endl;	
				std::cout << "innerShell.size: " << innerShell.size() << std::endl;	
				std::cout << "middleShell.size: " << middleShell.size() << std::endl;	
				std::cout << "outerShell.size: " << outerShell.size() << std::endl;	
				// return 1;
			} else {

				// testSurfaceVertices(tri);

				Eigen::MatrixXd h_fem, h_dec, h_opt, h_bndropt, h_decreg;
				Eigen::SparseMatrix<double> M, M_dec, M_decreg;
				std::string laplacesavepath = "";
				std::string regloadpath = "";
				if (saveLaplacians) {
					laplacesavepath = run_folder + run_name  + run_postfix;
				}
				if (regnoise == -2) {
					regloadpath  = run_folder + run_name;
				}

				if (!solveDirichletProblem(tri, reg, innerShell, middleShell, outerShell, h_fem, M, h_dec, M_dec, h_opt, h_bndropt, h_decreg, M_decreg, laplacesavepath, maxits, removeInner, addOpt, addBoundaryOpt, run_folder + run_name + run_postfix, false, regloadpath)){
					std::cout << "Error in Dirichlet Problem solving for " << run_name + run_postfix << std::endl;	
					return 1;
				}

				std::cout << "...write heat vals to file... " << std::endl;
				std::string res_out_path = run_folder + run_name + run_postfix + "heatvals.csv";
				std::ofstream feil;
				feil.open(res_out_path);
				feil << "middle shell indices" << std::endl;
				for(int i=0; i < middleShell.size() - 1; i++) feil << middleShell[i] << ", ";
				feil << middleShell[middleShell.size()-1] << std::endl;
				// -----------------------
				feil << "h_fem" << "," << "h_dec"; 
				if (addOpt)         feil << "," << "h_opt";
				if (addBoundaryOpt) feil << "," << "h_bndropt";
				if (reg || regnoise == -2)            feil << "," << "h_decreg";
				if (includeMass){
					feil << "," << "M" << "," << "M_dec";
					if (reg || regnoise == -2) feil << "," << "M_decreg";
				}
				feil << std::endl;
				for (int r = 0; r < h_fem.rows(); ++r) {
					feil << h_fem(r) << "," << h_dec(r);
					if (addOpt)         feil << "," << h_opt(r);
					if (addBoundaryOpt) feil << "," << h_bndropt(r);
					if (reg || regnoise == -2) feil << "," << h_decreg(r);
					if (includeMass) {
						feil << "," << M.coeff(r,r) << "," << M_dec.coeff(r,r);
						if (reg || regnoise == -2) feil << "," << M_decreg.coeff(r,r);
					}
					feil << std::endl;
				}
				feil.close();
				std::cout << "Finished the feil" << std::endl;
				// -----------------------
				// normalize the heat values:
				x = normalizeHeatValues(h_decreg);

			}

		} else if (heatdiffusion) {

			std::cout << "HEAT DIFFUSION EXPERIMENT" << std::endl;

			/* ################## HEAT DIFFUSION ################ */
			if (innerShell.size() < 1 || outerShell.size() < 1) {
				std::cout << "Heat Diffusion: Error with the shell points, abort." << std::endl;	
				std::cout << "innerShell.size: " << innerShell.size() << std::endl;	
				std::cout << "outerShell.size: " << outerShell.size() << std::endl;	
				// return 1;
			} else {

				Eigen::MatrixXd h_fem, h_dec, h_opt, h_bndropt, h_decreg;
				Eigen::SparseMatrix<double> M, M_dec, M_decreg;
				std::string laplacesavepath = "";
				std::string regloadpath = "";
				if (saveLaplacians) {
					laplacesavepath = run_folder + run_name  + run_postfix;
				}
				if (regnoise == -2) {
					regloadpath  = run_folder + run_name;
				}

				std::cout << "...solve heat diffusion" << std::endl;
				//double t = 0.001;// 0.005;
				solveHeatProblem(tri, reg, innerShell, outerShell, h_fem, h_dec, h_decreg, M, M_dec, M_decreg, regloadpath, heat_timestep);

				std::string heatgradOutPath_cells = run_folder + run_name + run_postfix +    "heatgradients_cells.csv";
				std::string heatgradOutPath_vertices = run_folder + run_name + run_postfix + "heatgradients_vertices.csv";

				bool includereg = (reg || regnoise == -2);
				calcAndWriteHeatGradsAndCentroidsAllCells(tri, h_fem, h_dec, h_decreg, includereg, heatgradOutPath_cells);
				calcAndWriteHeatGradsAllVertices(tri, h_fem, h_dec, h_decreg, includereg, heatgradOutPath_vertices);

				std::cout << "...calc distances to origin" << std::endl;
				Eigen::VectorXd D;
				tri.calcDistToPointAllVertices(D, CGAL::ORIGIN);
				std::cout << "D.size = " << D.size() << std::endl;
				std::cout << " h_fem.rows() = " << h_fem.rows() << std::endl;
				std::cout << " h_fem(0) " << h_fem(0) << std::endl;
				std::cout << " h_dec(0) " << h_dec(0) << std::endl;
				std::cout << " D(0) " << D(0) << std::endl;

				std::cout << "...write heat vals to file... " << std::endl;
				std::string res_out_path = run_folder + run_name + run_postfix + "heatvals.csv";
				std::ofstream feil;
				feil.open(res_out_path);
				// -----------------------
				feil << "h_fem" << "," << "h_dec"; 
				if (addOpt)         feil << "," << "h_opt";
				if (addBoundaryOpt) feil << "," << "h_bndropt";
				if (reg || regnoise == -2)            feil << "," << "h_decreg";
				if (includeMass){
					feil << "," << "M" << "," << "M_dec";
					if (reg || regnoise == -2) feil << "," << "M_decreg";
				}
				feil << ","<< "pointnorm";
				feil << std::endl;
				for (int r = 0; r < h_fem.rows(); ++r) {
					feil << h_fem(r) << "," << h_dec(r);
					if (addOpt)         feil << "," << h_opt(r);
					if (addBoundaryOpt) feil << "," << h_bndropt(r);
					if (reg || regnoise == -2) feil << "," << h_decreg(r);
					if (includeMass) {
						feil << "," << M.coeff(r,r) << "," << M_dec.coeff(r,r);
						if (reg || regnoise == -2) feil << "," << M_decreg.coeff(r,r);
					}
					feil << "," << D(r);
					feil << std::endl;
				}
				feil.close();
				std::cout << "Finished the feil" << std::endl;
				// -----------------------
				// normalize the heat values:
				x = normalizeHeatValues(h_fem);

			}



		
		} 

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

			viewer.launch();
		}
	run_idx++;
	}
}


