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


#include <Eigen/Dense>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>    

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
tetgenMeshSphere(CGALTriangulation<Kernel>& tri, int num_samples=30, int n_orbitpoints=3, std::string tetgenoptstring="")//, float q_maxre=-1, float q_inda=-1, float a_val=-1, int Oops=-1, int Olevel=-1)
{
    using namespace std;
    using namespace std::chrono;

	const double area = .1;
	double orbit_height=.5;

	std::cout << "Generate random samples" << std::endl;
	// randomPoints method has been moved to CGALMeshPolyhedron
	Eigen::MatrixXd randomsamples = randomPoints(num_samples);
	Eigen::MatrixXd orbitpoints   = randomPoints(n_orbitpoints);
    
    tetgenio in, out, out0;
    
    in.firstnumber = 0;

    in.numberofpoints = num_samples + n_orbitpoints + 1; 
    in.pointlist = new REAL[in.numberofpoints * 3];

	int offst = 0;

	// add origin
	in.pointlist[0]=0;
	in.pointlist[1]=0;
	in.pointlist[2]=0;
	offst+=3;

	// add normalized samples (on 3d shpere)
	for(int i=0; i< randomsamples.cols(); i++){
		for(int j=0; j<3; j++){
			in.pointlist[offst+j] = randomsamples(j,i) / randomsamples.col(i).norm();
		}
		offst+=3;
	}
	// add orbit points (on a mini sphere inside the other one)
	for(int i=0; i< orbitpoints.cols(); i++){
		for(int j=0; j<3; j++){
			in.pointlist[offst+j] = orbit_height * orbitpoints(j,i) / orbitpoints.col(i).norm();
		}
		offst+=3;
	}

	// GENERATE
	std::cout << "Step 1: generate starmesh" << std::endl;

    tetgenbehavior settings;
    string opts = string(""); //string("q1.414a") + to_string(area);
    settings.parse_commandline((char*)opts.c_str());
    settings.quiet = 1;
    tetrahedralize(&settings, &in, &out0);
	std::cout << "...done" << std::endl;


	// OPTIMIZE
	opts= "r" + tetgenoptstring;
    settings.parse_commandline((char*)opts.c_str());
    //settings.quiet = 1;
    tetrahedralize(&settings, &out0, &out);
    
    IndexedTetMesh mesh;
    
    mesh.vertices.resize(out.numberofpoints);
    copy_n(out.pointlist, 3 * out.numberofpoints, (double*)mesh.vertices.data() );
    
    mesh.tets.resize(out.numberoftetrahedra);
    copy_n(out.tetrahedronlist, 4 * out.numberoftetrahedra, (int*)mesh.tets.data() );
    
    mesh.convert(tri);
}

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

void replaceMeshByRegular(CGALTriangulation<Kernel> &tri, double variance, std::vector<int> &orbitinds, int &originind, double minVolume, bool boundary_only) {
		
	typedef CGAL::Regular_triangulation_vertex_base_3<Kernel> Vb0;
	typedef typename CGAL::Triangulation_vertex_base_with_info_3<int, Kernel, Vb0> VB;
	typedef CGAL::Regular_triangulation_cell_base_3<Kernel> Cb0;
	typedef typename CGAL::Triangulation_cell_base_with_info_3<int, Kernel, Cb0> CB;
	typedef CGAL::Triangulation_data_structure_3<VB, CB>  TriangulationDS;
	//typedef CGAL::Triangulation_3<Kernel, TriangulationDS_> Triangulation;
	typedef CGAL::Regular_triangulation_3<Kernel, TriangulationDS> Regular;
	typedef Kernel::Weighted_point_3 WPoint;

	std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, variance};
	
	std::vector< std::pair<WPoint,unsigned> > points;
    for (auto vh : tri.mesh.finite_vertex_handles()) {
		points.push_back( std::make_pair(WPoint(vh->point(),fabs(d(gen) )),vh->info()) );
    }

	Regular reg;
	reg.insert(points.begin(), points.end());
	std::cout << reg.is_valid() << std::endl;

	//handle infos: (finite vertex infos have been set when inserting them)
    reg.infinite_vertex()->info() = -1;
	int cnt = 0;
	for(auto it = reg.cells_begin(); it != reg.cells_end(); ++it)
	{
		if(reg.is_infinite(it)) it->info() = -1;
		else it->info() = cnt++;
	}

	// Translate to IndexedTetmesh
    IndexedTetMesh ret;
	int nv = reg.number_of_vertices();

	ret.vertices.resize(nv);
	std::unordered_map<int, int> idconversion;

	int inscounter = 0;
    for(auto it = reg.vertices_begin(); it != reg.vertices_end(); ++it)
        if(it->info() != -1){
			ret.vertices[inscounter][0] = it->point().x();
			ret.vertices[inscounter][1] = it->point().y();
			ret.vertices[inscounter][2] = it->point().z();
			idconversion[it->info()] = inscounter;
			inscounter++;
		}

	std::cout << "###############################################################" << std::endl;
	std::cout << "BouNDARY only: " << boundary_only << std::endl;
	std::cout << "###############################################################" << std::endl;

    for(auto it: reg.finite_cell_handles()){ // = reg.cells_begin(); it != reg.cells_end(); ++it)

		bool addCell = true;

		// check each cell if (!boundary_only):
		bool checkvolcell = false;
		if (!boundary_only){
			std::cout << "Check all Cells" << std::endl;		
			checkvolcell = true;
		}

		// if boundary only: check if cell is on boundary
		if (!checkvolcell) {
			for (int i = 0; i < 4 ; ++i) {
				if (reg.mirror_vertex(it, i)->info() == -1) {
					//std::cout << "Boundary cell " << std::endl;
					checkvolcell = true;
				}
			}
		}

		if (checkvolcell) {
			auto tet = reg.tetrahedron(it);
			double vol = tet.volume();
			//std::cout << "Vol:    " << vol               << std::endl;
			//std::cout << "Minvol: " << minVolume << std::endl;
			if (vol < minVolume){
				//std::cout << "Raus damit!" << std::endl;	
				addCell = false;
			}
		}

        if(it->info() == -1) addCell = false;

		if(addCell)
		{
            ret.tets.push_back(std::array<unsigned int, 4>{(unsigned int)idconversion[it->vertex(0)->info()],
															(unsigned int)idconversion[it->vertex(1)->info()],
															(unsigned int)idconversion[it->vertex(2)->info()],
															(unsigned int)idconversion[it->vertex(3)->info()] });
		}
	}
    
	CGALTriangulation<Kernel> newtri;
	ret.convert(newtri);

	// update
	tri = newtri;

	// convert origin and orbit indices
	if (originind >= 0) {
		int new_originind = -1;
		if (idconversion.find(originind) != idconversion.end()) {
			new_originind = idconversion[originind];
		} else {
			std::cout << "origin lost during noising" << std::endl;	
		}
		originind = new_originind;
	}

	std::vector<int> new_orbitinds;
	for (int i: orbitinds) {
		if (idconversion.find(i) != idconversion.end()) {
			new_orbitinds.push_back(idconversion[i]);	
		}
	}
	orbitinds = new_orbitinds;

}


// SET GLOBAL TO HANDLE IN CALLBACKS
enum Metric {minangle=0, amips, volume};
std::string FILENAME_base = "";
std::string FILENAME="";
Eigen::MatrixXd facecolors;
std::map<Metric, Eigen::MatrixXd> cellcolors;
std::vector<int>  faceids; 
Metric metric_shown;

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
	  for (int i=0; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
	  viewer.data().set_colors(facecolors);

	  if(metric_shown == minangle) {
		viewer.launch_shut(); 
	  } 
  }

  if (key == '0') {
	std::cout << "View angle: " << std::endl;
	std::cout << viewer.core().camera_view_angle << std::endl; 
	std::cout << "Camera Up:" << std::endl;
	std::cout << viewer.core().camera_up << std::endl; 
	std::cout << "Camera Eye:" << std::endl;
	std::cout << viewer.core().camera_eye << std::endl; 
	std::cout << "Camera base trans:" << std::endl;
	std::cout << viewer.core().camera_base_translation << std::endl; 
	std::cout << viewer.core().camera_translation << std::endl; 
	std::cout << viewer.core().camera_center << std::endl; 
	std::cout << viewer.core().trackball_angle.coeffs() << std::endl;
  }

  if (key == '9') {
	  viewer.core().camera_view_angle += 5;
  }
  if (key == '8') {
	  viewer.core().camera_translation *= 1.1;
  }
}

bool
adjustPointsOriginAndOrbit(CGALTriangulation<Kernel> &tri, std::vector<int> &changed_inds, std::vector<double> &changed_by, double orbit_dist =.5, double epsfactor=.1) {

	int nv = tri.mesh.number_of_vertices();
    std::vector<double> dists(nv, 0.);
	changed_inds.clear();
	changed_by.clear();

    typedef CGALTriangulation<Kernel>::Point Point;
    typedef CGALTriangulation<Kernel>::Cell_handle Cell_handle;
    Point origin = Point(0., 0., 0.);

    double dist;
    int minind = 0;

	double mean_edge_length = sqrt(tri.meanEdgeLengthSquared());
	double eps = epsfactor * mean_edge_length;

    // Go through all points and set the points that lie in the orbit region to the orbit distance
    // Also save the index of the point closest to the origin to set it exactly to the origin later
    for (auto a: tri.mesh.finite_vertex_handles()) {
        int v_ind = a->info();
        dist = sqrt(CGAL::squared_distance(origin, a->point()));
        dists[v_ind] = dist;
        if (dist <= dists[minind]) minind = v_ind;
        if (orbit_dist - eps <= dist  && dist < orbit_dist + eps ) {
            CGALTriangulation<Kernel>::Point normed_p = Point(a->point().x() / dist * orbit_dist, a->point().y() / dist * orbit_dist, a->point().z() / dist * orbit_dist);

			// check that the move does not flip a cell:
			std::vector<Cell_handle> oui;
			tri.mesh.incident_cells(a, std::back_inserter(oui));

			bool allsamesign = true;
			for (auto c: oui) {
				Point p[4];
				int vertex_ind_in_cell = -1;
				for (int i = 0; i < 4; ++i) {
					p[i] = c->vertex(i)->point();
					if (p[i] == a->point()) {
						vertex_ind_in_cell = i;
					}
				}
				Kernel::FT vpre  = CGAL::volume(p[0], p[1], p[2], p[3]);
				p[vertex_ind_in_cell] = normed_p;
				Kernel::FT vpost = CGAL::volume(p[0], p[1], p[2], p[3]);

				bool samesign = (vpre < 0 && vpost < 0) || ( vpre>0 && vpost > 0) || (vpre == 0 && vpost == 0);
				if (!samesign) allsamesign = false;
			}

			if (allsamesign) {
				double cb = sqrt(CGAL::squared_distance(normed_p, a->point())); 
				std::cout << "...move orbitpoint " << v_ind <<" by " << cb << std::endl;
				a->set_point(normed_p);
				dists[v_ind] = sqrt(CGAL::squared_distance(origin, a->point()));
				changed_inds.push_back(v_ind);
				changed_by.push_back(cb);
			} else {
				std::cout << ".....skipped moving an orbitpoint (it would flip a cell)" << std::endl;	
			}
        }
    }

    int n_orbit_points = changed_inds.size();
    std::cout << "Found " << n_orbit_points << " orbit points" <<  std::endl;

	/*
	if (n_orbit_points < 1) {
		std::cout << "Error: orbit eps too small, no points in range, abort" << std::endl;	
		return false;
	}
	*/

	double cb;
    // Set closest point to (0,0,0) to the origin
    for (auto a: tri.mesh.finite_vertex_handles()) {
        int v_ind = a->info();
        if (v_ind == minind){
			cb =  sqrt(CGAL::squared_distance(origin, a->point()));
            std::cout << "set point " << a->point() << "to the origin (moved it by " << cb << ")" << std::endl;
            a->set_point(origin);
            dists[v_ind] = sqrt(CGAL::squared_distance(origin, a->point()));
			std::cout << "Minind: " << minind << std::endl;
			//add origin to the orbit inds
        }
    }

    // Make sure none of the orbit points has been set to the origin
    for (auto o_ind : changed_inds){
        if (o_ind == minind) {
            std::cout << "ERROR: an orbit point is in the origin" << std::endl;
            return false;
        }
    }
	changed_inds.push_back(minind);
	changed_by.push_back(cb);

	return true;
}


void normalizeSphereBorder(CGALTriangulation<Kernel> &tri) {
	
	typedef typename Kernel::Point_3 Point;
	typedef typename CGALTriangulation<Kernel>::Vertex_handle Vertex_handle;

	std::vector<Vertex_handle> out;
	tri.mesh.finite_adjacent_vertices(tri.mesh.infinite_vertex(), std::back_inserter(out));
	for (auto p : out) {
		double length = sqrt(CGAL::squared_distance(Point(0., 0., 0.), p->point()));	
		p->point() = Point(p->point().x() / length, p->point().y() / length, p->point().z() / length);
	}

	/*
	for (auto p : tri.mesh.finite_vertex_handles()) {
		double length = sqrt(CGAL::squared_distance(Point(0., 0., 0.), p->point()));	
		if (fabs(1 - length) < 1e-3) {
			std::cout << "ar" << std::endl;
			//p->point() = Point(p->point().x() / length, p->point().y() / length, p->point().z() / length);
			p->point() = Point(3., p->point().y() / length, p->point().z() / length);
		}
	}
	*/

}


int main(int argc, char *argv[])
{

	if (argc < 8) {
		std::cout << "Usage: " << argv[0] << " cellSize reRatio facetSize use_lloyd use_perturb use_exude regnoise silent minvol filter_border_only (ma_min_thr ma_max_thr invert)" << std::endl;
		return EXIT_FAILURE;
	}

    CGALTriangulation<Kernel> tri;
    typedef CGALTriangulation<Kernel>::Point Point;
    
	// sphere generation options (with default values)
	double  edgeprob		=   0.; 

	// #########################################
	std::cout << "SPHERE CREATION" << std::endl;
	// #########################################

	bool silent=true;

	// tetlib CELLSIZE CERATIO LLOYD PERTURB EXUDE NFLIPS
	meshingOptions mOptions;
	double cellSize    = std::stod(argv[1]);
	mOptions.cell_size  = cellSize;
	mOptions.cell_radius_edge_ratio = std::stod(argv[2]);
	double facetSize   = std::stod(argv[3]);
	mOptions.facet_size = facetSize;
	if (atoi(argv[4])) mOptions.opt_lloyd   = true;
	if (atoi(argv[5])) mOptions.opt_perturb = true;
	if (atoi(argv[6])) mOptions.opt_exude   = true;
	double regnoise = std::stod(argv[7]);
	if (atoi(argv[8])) silent = false;

	double minVolume     = 0.;
	bool   boundary_only = true;
	if (argc >= 10) {
		minVolume = std::stod(argv[9]);
	}
	if (argc >= 11) {
		if (!atoi(argv[10])){
			boundary_only = false;	
			std::cout << "...aply minvol filter to all cells" << std::endl;
		}
	}

	// minangle colormap mapping options
	double minangle_min_threshold = 0.;
	double minangle_max_threshold = 70.5;
	bool invert_minangle_cmap = true;

	if (argc == 14) {
		minangle_min_threshold = std::stod(argv[11]);	
		minangle_max_threshold = std::stod(argv[12]);	
		if (!atoi(argv[13])) invert_minangle_cmap = false;
	}


	bool addOrbitpointsMEL = true;

	//int min_orbitpoints  	= std::atoi(argv[7]);
	//int n_flips				= std::atoi(argv[8]);
	bool addorigin = false;
	int min_orbitpoints = -1;
	int n_flips = 0;

	// TODO: torus argument not working in meshpolyhedron
	bool torus=false;
	//if (atoi(argv[7])) torus = true;

	FILENAME_base = "Sphere_";
	if (torus) {
		FILENAME_base = "Torus_";
	}

	for (int i=0; i < 7; ++i) FILENAME_base += argv[i+1] + std::string("_");

	meshSphere<CGAL::Exact_predicates_inexact_constructions_kernel>(tri,mOptions,torus);

	if (!torus) {
		normalizeSphereBorder(tri);	
	}

	int origin_ind = -1;
	double origin_changedby=-1;
	std::vector<int> orbit_indices;
	std::vector<double> orbit_changedby;

	if (addOrbitpointsMEL) {
		// #########################################
		std::cout << "Add origin and orbitpoints ifo MEL" << std::endl;
		// #########################################
		adjustPointsOriginAndOrbit(tri, orbit_indices, orbit_changedby);
		std::cout << "Added " << orbit_indices.size() - 1 << " orbitpoints (MELMELMEL)" << std::endl;

		origin_ind = orbit_indices[orbit_indices.size()-1];
		origin_changedby = orbit_changedby[orbit_indices.size()-1];
		orbit_indices.pop_back();
		orbit_changedby.pop_back();

		// Write orbitpoints and amount they were moved to a file
		std::ofstream feil;
		feil.open("out/" + FILENAME_base + "orbitpoints.csv");
		feil << "orbit_idx" << ", " << "moved_by" << std::endl;
		// add origin
		feil << origin_ind << ", " << origin_changedby << std::endl;
		// add orbitpoints
		for (int i=0; i < orbit_indices.size(); ++i) {
			feil      << orbit_indices[i] << ", " << orbit_changedby[i] << std::endl; 		
		}
		feil.close();

	}

	if (regnoise >= 0) {
		// #########################################
		std::cout << "Noise with Reg Triangulation" << std::endl;
		// #########################################

		// set variance st that the std deviation of the reg weights is in the magnitude of mels * rgnoise
		double mels = tri.meanEdgeLengthSquared();
		double variance = mels * regnoise;  // (mels*mels*regnoise) * (mels*mels*regnoise);
		//std::cout << "MELS: " << mels << std::endl;
		//std::cout << "VAR : " << variance << std::endl;
		//double variance = regnoise;

		double origin_ind_pre = origin_ind;
		std::vector<int> orbit_indices_pre = orbit_indices;

		// old:
		//replaceMeshByRegular(tri, variance, orbit_indices, origin_ind, minVolume, boundary_only);
		// new:
		tri.replaceMeshByRegular(variance, orbit_indices, origin_ind, minVolume, boundary_only);

		std::cout << "Origin index pre : " << origin_ind_pre << std::endl;
		std::cout << "Origin index post: " << origin_ind     << std::endl;

		std::cout << "Orbit indices.size pre : " << orbit_indices_pre.size() << std::endl;
		std::cout << "Orbit indices.size post: " << orbit_indices.size()     << std::endl;

		// Write orbitpoints and amount they were moved to a file
		std::ofstream feil;
		feil.open("out/" + FILENAME_base + "orbitpoints.csv");
		feil << "orbit_idx" << ", " << "moved_by" << std::endl;
		// add origin
		feil << origin_ind << ", " << origin_changedby << std::endl;
		// add orbitpoints
		for (int i=0; i < orbit_indices.size(); ++i) {
			feil      << orbit_indices[i] << ", " << orbit_changedby[i] << std::endl; 		
		}
		feil.close();
	}

	/*
	if (false) {
	//if (addorigin and min_orbitpoints < 0) {
		// #########################################
		std::cout << "ADD ORIGIN AND ORIGIN ONLY" << std::endl;
		// #########################################

		for (auto a: tri.mesh.finite_vertex_handles()) {
			int v_ind = a->info();
			dist = sqrt(CGAL::squared_distance(origin, a->point()));
			dists[v_ind] = dist;
			if (dist <= dists[minind]) minind = v_ind;
			if (orbit_dist - eps <= dist  && dist < orbit_dist + eps ) {
				CGALTriangulation<Kernel>::Point normed_p = Point(a->point().x() / dist * orbit_dist, a->point().y() / dist * orbit_dist, a->point().z() / dist * orbit_dist);

				// check that the move does not flip a cell:
				std::vector<Cell_handle> oui;
				tri.mesh.incident_cells(a, std::back_inserter(oui));

				bool allsamesign = true;
				for (auto c: oui) {
					Point p[4];
					int vertex_ind_in_cell = -1;
					for (int i = 0; i < 4; ++i) {
						p[i] = c->vertex(i)->point();
						if (p[i] == a->point()) {
							vertex_ind_in_cell = i;
						}
					}
					Kernel::FT vpre  = CGAL::volume(p[0], p[1], p[2], p[3]);
					p[vertex_ind_in_cell] = normed_p;
					Kernel::FT vpost = CGAL::volume(p[0], p[1], p[2], p[3]);

					bool samesign = (vpre < 0 && vpost < 0) || ( vpre>0 && vpost > 0) || (vpre == 0 && vpost == 0);
					if (!samesign) allsamesign = false;
				}

				if (allsamesign) {
					double cb = sqrt(CGAL::squared_distance(normed_p, a->point())); 
					std::cout << "...move orbitpoint " << v_ind <<" by " << cb << std::endl;
					a->set_point(normed_p);
					dists[v_ind] = sqrt(CGAL::squared_distance(origin, a->point()));
					changed_inds.push_back(v_ind);
					changed_by.push_back(cb);
				} else {
					std::cout << ".....skipped moving an orbitpoint (it would flip a cell)" << std::endl;	
				}
			}
		}

	}
	*/

	if (false) {
	//if (n_flips > 0) {
		// #########################################
		std::cout << "FLIPPING EDGES" << std::endl;
		// #########################################
		tri.performRandomFlips(n_flips, 2*n_flips, edgeprob);
	}

	if (false) {
	//if (min_orbitpoints >= 0) {
		// #########################################
		std::cout << "GENERATE ORBIT POINTS" << std::endl;
		// #########################################
		std::cout << "look for orbit points in small radius and shift them" << std::endl;	
		std::cout << "Start with an epsilon of cellSize^(1_3) / 500" << std::endl;
		double epsilon = cbrt(cellSize) / 500;

		origin_ind = -1;
		origin_changedby=-1;
		std::vector<int> changed_indices;
		std::vector<double> changed_by;

		int n_orbitpoints=0;
		bool any_new_ones=true;
		while (n_orbitpoints < min_orbitpoints && epsilon <= 0.5 && any_new_ones) {
			adjustPointsOriginAndOrbit(tri, changed_indices, changed_by, .5, epsilon);
			n_orbitpoints = changed_indices.size()-1;
			for (double cb : changed_by) std::cout << cb << " ";
			std::cout << std::endl;

			if (-1  == origin_ind) {
				// origin is only set once
				//std::cout << "origin setter " << changed_by[changed_indices.size()-1] << std::endl;
				origin_ind = changed_indices[changed_indices.size()-1];
				origin_changedby = changed_by[changed_indices.size()-1];
			}
			
			// update orbit index and changedby for all NEW orbitpoints
			any_new_ones = false;
			for (int j =0; j < n_orbitpoints; ++j){
				int cind = changed_indices[j];
				bool new_one=true;
				for (int i: orbit_indices) {
					if (i == cind) {
						new_one = false;	
						break;
					}
				}
				if (new_one) {
					any_new_ones=true;
					//std::cout << "New point: (" << cind << "|" << changed_by[j] << ")" << std::endl;
					orbit_indices.push_back(cind);
					orbit_changedby.push_back(changed_by[j]);
				}
			}

			if (!any_new_ones) std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!! No New points can be found with larger radius, stop" << std::endl;

			if(n_orbitpoints < min_orbitpoints) epsilon += cbrt(cellSize) / 2000;
		}


		// Write orbitpoints and amount they were moved to a file
		std::ofstream feil;
		feil.open("out/" + FILENAME_base + "orbitpoints.csv");
		feil << "orbit_idx" << ", " << "moved_by" << std::endl;
		// add origin
		feil << origin_ind << ", " << origin_changedby << std::endl;
		// add orbitpoints
		for (int i=0; i < orbit_indices.size(); ++i) {
			feil      << orbit_indices[i] << ", " << orbit_changedby[i] << std::endl; 		
		}
		feil.close();

		if (false) {
		// TEST WRITEOUT: --------------------------------
		// - save orbitpoints
		std::vector<Point> orbitpnts;
		for (int i : changed_indices){
			for (auto vh : tri.mesh.finite_vertex_handles()) {
				if (vh->info() == i) {
					orbitpnts.push_back(vh->point());
				}
			}
		} 

		// store mesh
		std::string outfilename = "out/" + FILENAME_base;
		tri.write(outfilename);
		CGALTriangulation<Kernel> tri_l;
		tri_l.read(outfilename);
		
		// - compare  orbitpoints
		int j = 0;
		for (auto i :changed_indices){
			for (auto vh : tri_l.mesh.finite_vertex_handles()) {
				if (vh->info() == i) {
					if ( (vh->point() - orbitpnts[j]).squared_length() < 1e-8) {
						std::cout << " Points the same" << std::endl;
					} else {
						std::cout << " Points NOT the same: " << vh->point() << orbitpnts[j] << std::endl;
					}
					j++;
				}
			}
		} 
		
		// -----------------------------------------------
		}
	}

	// #########################################
	std::cout << "WRITE TO FILE"  << std::endl;
	// #########################################
	std::string outfile = "out/" + FILENAME_base + ".meshfile";
	tri.write(outfile);

	std::ofstream mNfile;
	std::string meshNamesFile= "out/meshes.txt";
	mNfile.open(meshNamesFile, std::ios_base::app); // append instead of overwrite
	mNfile << FILENAME_base << std::endl; 
	mNfile.close();


	// #########################################
	std::cout << "METRICS"  << std::endl;
	// #########################################
	std::cout << "Calc metrics" << std::endl;
	metric_shown = minangle;
	std::map<Metric, Eigen::VectorXd> cell_metrics;
	std::map<Metric, std::string> metric_names;
	metric_names[minangle] = "minangle";
	metric_names[amips] = "amips";
	metric_names[volume] = "volume";
	FILENAME = FILENAME_base + metric_names[metric_shown];

	Eigen::VectorXd Vol, Minang, Amips; 
	tri.calcVolumeAllCells(Vol);
	tri.calcMinAngleAllCells(Minang);
	tri.calcAMIPSAllCells(Amips);
	cell_metrics[volume]=Vol;
	cell_metrics[minangle]=Minang;
	cell_metrics[amips]=Amips;

	// write values to file
	std::ofstream feil;
	feil.open("out/" + FILENAME_base + "metrics.csv");
	feil << metric_names[minangle] << "," << metric_names[amips] << "," << metric_names[volume] << std::endl;
	for(int i; i < cell_metrics[volume].size(); i++) {
		feil << cell_metrics[minangle](i) << "," << cell_metrics[amips](i) << "," << cell_metrics[volume](i) << std::endl;
	}
	feil.close();

	if (!silent) {

		bool normalize=false;
		double amips_max = 50;


		if (!normalize) {
			for (int i=0; i < cell_metrics[minangle].size(); i++){

				double value = cell_metrics[minangle][i];
				// clip
				value = (value > minangle_min_threshold)? value : minangle_min_threshold;
				value = (value < minangle_max_threshold)? value : minangle_max_threshold;
				// normalize
				value = (value-minangle_min_threshold) / (minangle_max_threshold- minangle_min_threshold)  ;
				// (invert?) and save	
				cell_metrics[minangle][i] = (invert_minangle_cmap)? 1 - value : value;
			}

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
		for (int i =0; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
		
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
			//   menu.draw_viewer_menu();
			
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
					for (int i=0; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
					viewer.data().set_colors(facecolors);

					FILENAME = FILENAME_base + metric_names[metric_shown];
				}
			}
			
			// Add new group
			if (ImGui::CollapsingHeader("Cut View", ImGuiTreeNodeFlags_DefaultOpen))
			{
				double mi = -1.8;
				double ma = 1.8;
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
						for (int i=0; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
						viewer.data().set_colors(facecolors);
						//viewer.data(cutMeshId).uniform_colors(ambient, diffuse, specular);
						//viewer.data(cutMeshId).show_texture = 1;
						/*
						if (showValues) {	
							for(int i = 0; i < F.rows(); ++i) {
								const Eigen::Vector3d FaceCenter( (V(F(i,0), 0)+ V(F(i,1), 0)+ V(F(i,2), 0))/ 3.,  (V(F(i,0), 1)+ V(F(i,1), 1)+ V(F(i,2), 1))/ 3.,  (V(F(i,0), 2)+ V(F(i,1), 2)+ V(F(i,2), 2))/ 3.);
								viewer.data().add_label(FaceCenter,std::to_string(cell_metrics[metric][faceids[i]]));
							}
						} 
						*/
						
						Eigen::VectorXd x2(ids.size());
						
						for(int i = 0; i < ids.size(); ++i)
							x2(i) = x(ids[i]);
						
						//  Eigen::MatrixXd color;
						//  igl::colormap(igl::COLOR_MAP_TYPE_INFERNO, x2, false, color);
						//viewer.data().set_colors(color);
						
						/*
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

			
		/* Add labels for debugging
		for(int i = 0; i < F.rows(); ++i) {
			const Eigen::Vector3d FaceCenter( (V(F(i,0), 0)+ V(F(i,1), 0)+ V(F(i,2), 0))/ 3.,  (V(F(i,0), 1)+ V(F(i,1), 1)+ V(F(i,2), 1))/ 3.,  (V(F(i,0), 2)+ V(F(i,1), 2)+ V(F(i,2), 2))/ 3.);
			viewer.data().add_label(FaceCenter,std::to_string(faceids[i]));
		}
		*/

		//viewer.data().uniform_colors(ambient, diffuse, specular);
		//viewer.data().show_texture = 1;
		/*
		setTexture("../data/tex.png", viewer.data());
		Eigen::MatrixXd UV(V.rows(), 2);
		viewer.data().set_uv(UV.setZero());
		*/


		/*
		viewer.core().trackball_angle.x() = 0.121685;
		viewer.core().trackball_angle.y() = 0.335208;
		viewer.core().trackball_angle.z() = 0.0437081;
		viewer.core().trackball_angle.w() = 0.93323;
		*/

		viewer.core().trackball_angle.x() = 0.0727511;
		viewer.core().trackball_angle.y() = 0.861667;
		viewer.core().trackball_angle.z() = 0.129162;
		viewer.core().trackball_angle.w() = 0.485339;
 
		//screenshot(viewer, "dadada");

		viewer.launch();
	}

}
