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

std::string FILENAME="";

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    // Allocate temporary buffers
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

    // Draw the scene in the buffers
    viewer.core().draw_buffer(
      viewer.data(),false,R,G,B,A);

    // Save it to a PNG
    igl::png::writePNG(R,G,B,A, "out/" + FILENAME +  "out.png");
  }
}

bool
adjustPointsOriginAndOrbit(CGALTriangulation<Kernel> &tri, std::vector<int> &changed_inds, double orbit_dist =.5, double eps=.001) {

	int nv = tri.mesh.number_of_vertices();
    std::vector<double> dists(nv, 0.);
	changed_inds.clear();

    typedef CGALTriangulation<Kernel>::Point Point;
    typedef CGALTriangulation<Kernel>::Cell_handle Cell_handle;
    Point origin = Point(0., 0., 0.);

    double dist;
    int minind = 0;

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
				std::cout << "...move orbitpoint by " << sqrt(CGAL::squared_distance(normed_p, a->point())) << std::endl;
				a->set_point(normed_p);
				dists[v_ind] = sqrt(CGAL::squared_distance(origin, a->point()));
				changed_inds.push_back(v_ind);
			} else {
				std::cout << ".....skipped moving an orbitpoint (it would flip a cell)" << std::endl;	
			}
        }
    }
    int n_orbit_points = changed_inds.size();
    std::cout << "Found " << n_orbit_points << " orbit points" <<  std::endl;
	if (n_orbit_points < 1) {
		std::cout << "Error: orbit eps too small, no points in range, abort" << std::endl;	
		return false;
	}

    // Set closest point to (0,0,0) to the origin
    for (auto a: tri.mesh.finite_vertex_handles()) {
        int v_ind = a->info();
        if (v_ind == minind){
            std::cout << "set point " << a->point() << "to the origin (moved it by " << sqrt(CGAL::squared_distance(origin, a->point())) << ")" << std::endl;
            a->set_point(origin);
            dists[v_ind] = sqrt(CGAL::squared_distance(origin, a->point()));
        }
    }

    // Make sure none of the orbit points has been set to the origin
    for (auto o_ind : changed_inds){
        if (o_ind == minind) {
            std::cout << "ERROR: an orbit point is in the origin" << std::endl;
            return false;
        }
    }

	//add origin to the orbit inds
	changed_inds.push_back(minind);

	return true;
}


int main(int argc, char *argv[])
{
    CGALTriangulation<Kernel> tri;
    typedef CGALTriangulation<Kernel>::Point Point;
    
	// sphere generation options (with default values)
	int n_orbitpoints		=   -1; 
	int n_flips				=   0 ;
	double  edgeprob		=   0.; 
	double maxPointMove     =   0.001;
	
	std::cout << "START" << std::endl;

	std::string FILENAME_base = "";
	for (int i=0; i < argc-1; ++i) FILENAME_base += argv[i+1] + std::string("_");

	// tetlib NORBITPOINTS CELLSIZE RTOERATIO NFLIPS MAXPOINTMOVE LLOYD PERTURB EXUDE
	meshingOptions mOptions;
	if(argc >= 2){
		mOptions.n_orbitpoints = -1;  //atoi(argv[1]);
		if (argc>=3) {
			mOptions.cellSize = std::stod(argv[2]);
			if (argc>=4) {
				mOptions.cell_radius_edge_ratio = std::stod(argv[3]);
				if(argc>=5) {
					n_flips = std::atoi(argv[4]);
					if(argc>=6) {
							maxPointMove= std::stod(argv[5]);
						if(argc>=7) {
							if (atoi(argv[6])) mOptions.opt_lloyd=true;
							if(argc>=8) {
								if (atoi(argv[7])) mOptions.opt_perturb=true;
								if(argc>=9) {
									if (atoi(argv[8])) mOptions.opt_exude=true;
								}
							}
						}
					}
				}
			}
		}
	}

	meshSphere<CGAL::Exact_predicates_inexact_constructions_kernel>(tri,mOptions);


	std::cout << "Finished sphere creation" << std::endl;

	if (n_orbitpoints < 0) {
		std::cout << "n_orbit points < 0, look for orbit points in small radius and shift them" << std::endl;	

		std::vector<int> changed_indices;
		// write config - orbit points to file
		std::ofstream cfile;
		cfile.open("out/" "configs.txt", std::ios_base::app);
		for (int i=0; i < argc; ++i) cfile << argv[i] << ","; 
		cfile << std::endl;

		if (adjustPointsOriginAndOrbit(tri, changed_indices, .5, maxPointMove)) {
			std::cout << "changed " << changed_indices.size() << "points" << std::endl;

			cfile << changed_indices.size()-1 << std::endl;
			cfile.close();

			FILENAME_base += std::to_string( changed_indices.size() - 1) + "_";

		} else {
			std::cout << "No orbit points found, exit." << std::endl;
			cfile << "0" << std::endl;
			cfile.close();

			return EXIT_FAILURE;
		}

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


	/*
	std::cout << "Finished Sphere Creation (" << mode << ")" << std::endl;
	for (auto v: tri.mesh.finite_vertex_handles()) {
		std::cout << v->point() << std::endl;
	}
	*/

	//std::cout << "... perform flips " << std::endl;
	tri.performRandomFlips(n_flips, 2*n_flips, edgeprob);

	std::cout << "Calc metrics" << std::endl;
	enum Metric {minangle=0, amips, volume};
	Metric metric_shown = minangle;
	std::map<Metric, Eigen::VectorXd> cell_metrics;
	std::map<Metric, Eigen::MatrixXd> cellcolors;
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


	Eigen::MatrixXd facecolors;

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
	std::vector<int>  faceids = cmreturn[1];
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
				for (int i; i < faceids.size(); i++) facecolors.row(i) = cellcolors[metric_shown].row(faceids[i]);
				viewer.data().set_colors(facecolors);

				FILENAME = FILENAME_base + metric_names[metric_shown];
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

    viewer.launch();
}
