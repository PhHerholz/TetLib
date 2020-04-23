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
    igl::png::writePNG(R,G,B,A, filename +  "out.png");

}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{

	if (key == '0') {
		std::cout << "Trackball coords: ";
		std::cout << viewer.core().trackball_angle.x() << " "; 
		std::cout << viewer.core().trackball_angle.y() << " "; 
		std::cout << viewer.core().trackball_angle.z() << " "; 
		std::cout << viewer.core().trackball_angle.w() << std::endl; 

		std::cout << viewer.core().trackball_angle.coeffs() << std::endl;
	}


	if (key == '2')
	{
		std::cout << FILENAME_base << std::endl;
		screenshot(viewer, FILENAME);

		/*
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

double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}


void replaceMeshByRegular(CGALTriangulation<Kernel> &tri, std::vector<int> &orbitinds, int originind, double variance) {
		
	typedef CGAL::Regular_triangulation_vertex_base_3<Kernel> Vb0;
	typedef typename CGAL::Triangulation_vertex_base_with_info_3<int, Kernel, Vb0> VB;
	//typedef typename CGAL::Triangulation_vertex_base_with_info_3<int, Kernel> VB_;
	typedef CGAL::Regular_triangulation_cell_base_3<Kernel> Cb0;
	typedef typename CGAL::Triangulation_cell_base_with_info_3<int, Kernel, Cb0> CB;
	//typedef typename CGAL::Triangulation_cell_base_with_info_3<int, Kernel> CB_;
	typedef CGAL::Triangulation_data_structure_3<VB, CB>  TriangulationDS;
	//typedef CGAL::Triangulation_data_structure_3<VB, CB>  TriangulationDS_;
	//typedef CGAL::Triangulation_3<Kernel, TriangulationDS_> Triangulation;
	typedef CGAL::Regular_triangulation_3<Kernel, TriangulationDS> Regular;
	typedef Kernel::Weighted_point_3 WPoint;

	std::vector< std::pair<WPoint,unsigned> > points;
	std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, variance};
	
    for (auto vh : tri.mesh.finite_vertex_handles()) {
		points.push_back( std::make_pair(WPoint(vh->point(),fabs(d(gen) )),vh->info()) );
    }

    
	Regular reg;
	reg.insert(points.begin(), points.end());
	std::cout << reg.is_valid() << std::endl;

    reg.infinite_vertex()->info() = -1;
    
	/*
	std::vector<Point> pts;
	std::vector<int>   vinds;
	for (auto vh : tri.mesh.finite_vertex_handles()) {
		pts.push_back(vh->point());
		vinds.push_back(vh->info());
	}

	std::cout << "Vinds.size = " << vinds.size() << std::endl;

	std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0, variance};

	std::cout <<" 1" << std::endl;

	Regular reg;
	int cnt = 0;
	reg.infinite_vertex()->info() = -1;


	for(int i=0; i<pts.size(); ++i)
	{
		std::cout << "I " << i << std::endl;
		float rndm = fabs(d(gen));
		std::cout << "rndm:  " << rndm << std::endl;
		std::cout << "vs: " << vinds.size()<< std::endl;
		std::cout << "vpts: " << pts.size()<< std::endl;
		std::cout << "i: " <<  i  << std::endl;

		//std::cout << vinds[i] << std::endl;
		std::cout << pts[i] ;

		std::cout << "pst" << std::endl;
		//Point p = pts[i];
		std::cout << "pr-";
		WPoint wp = WPoint(pts[i],rndm);
		std::cout << "r";
		auto h = reg.insert(wp);
		std::cout << "-i";
		h->info() = vinds[i];
		std::cout << "-w" << std::endl;;
	}
	std::cout <<" 2" << std::endl;
	*/

    int cnt = 0;
    for(auto it = reg.cells_begin(); it != reg.cells_end(); ++it)
    {
        if(reg.is_infinite(it)) it->info() = -1;
        else it->info() = cnt++;
    }


	//std::cout << "Inserted " << std::endl;
	//std::cout << " Isvalid: " << reg.is_valid() << std::endl;

	cnt = 0;
	for(auto it = reg.cells_begin(); it != reg.cells_end(); ++it)
	{
		if(reg.is_infinite(it)) it->info() = -1;
		else it->info() = cnt++;
	}
	// tri.mesh = reg;
	//std::cout << "Converted Mesh to a basic Regular Delaunay" << std::endl;

	/*
    for(auto it = reg.vertices_begin(); it != reg.vertices_end(); ++it) {
		std::cout << it->info() << " ";//  << std::endl;	
	}
	std::cout << std::endl;	
	*/

	//CGAL::draw(reg);

    IndexedTetMesh ret;
	int nv = reg.number_of_vertices();
	//std::cout << "NV: " << nv << std::endl;

	ret.vertices.resize(nv);

	std::unordered_map<int, int> idconversion;

	bool contains_origin = false;
    
	int inscounter = 0;
    for(auto it = reg.vertices_begin(); it != reg.vertices_end(); ++it)
        if(it->info() != -1){
			ret.vertices[inscounter][0] = it->point().x();
			ret.vertices[inscounter][1] = it->point().y();
			ret.vertices[inscounter][2] = it->point().z();
			idconversion[it->info()] = inscounter;
			inscounter++;
            //ret.vertices.push_back(std::array<double, 3>{it->point().x(), it->point().y(), it->point().z()});

			if (it->info() == originind) {
				std::cout << "ORIGIN IN !!!!!!!!!!!!!!!!!!!!!!!! " << std::endl;
				contains_origin=true;
			}
		}

	//std::cout << "Reg Vertices: " << ret.vertices.size() << std::endl;
	//std::cout << "Finite Cells: " << reg.number_of_finite_cells() << std::endl;

    for(auto it: reg.finite_cell_handles()) // = reg.cells_begin(); it != reg.cells_end(); ++it)
        if(it->info() != -1)
		{
			//std::cout << it->info();
            ret.tets.push_back(std::array<unsigned int, 4>{(unsigned int)idconversion[it->vertex(0)->info()],
															(unsigned int)idconversion[it->vertex(1)->info()],
															(unsigned int)idconversion[it->vertex(2)->info()],
															(unsigned int)idconversion[it->vertex(3)->info()] });
			//std::cout << "-+" << std::endl;
		}

	//std::cout << "built map" << std::endl;

	std::vector<int> new_orbitinds;
	for (int i: orbitinds) {
		if (idconversion.find(i) != idconversion.end()) {
			new_orbitinds.push_back(idconversion[i]);	

			std::cout << "wrr" << std::endl;
			std::cout << (idconversion.find(i) != idconversion.end() ) << std::endl;
		}
	}

	if (contains_origin) {
		int new_originind = idconversion[originind];
		std::cout << "New originind: " << new_originind << std::endl;
	} else {
		std::cout << "ARGH " << idconversion[originind] << std::endl;	
		std::cout << (idconversion.find(originind) == idconversion.end()) << std::endl;
		std::cout << (idconversion[-4]) << std::endl;
	}

	//std::cout << "Converted Regular Delauney to IndexedTetmesh " << std::endl;

    //ret.write("../dbg.tet");
    
	CGALTriangulation<Kernel> newtri;
	ret.convert(newtri);

	// update
	tri = newtri;
	orbitinds = new_orbitinds;
}

int main(int argc, char *argv[])
{

	// ################
	// # ARG HANDLING # 
	// ################
	if (argc < 2) {
		std::cout << "usage: argv[0] file (minangle_max_thr minangle_min_thr minangle_invert_cmap)" << std::endl;  //run_folder variance r_postfix (1 for gui output)" << std::endl;
		return 0;
	}

	std::string filepath = argv[1];
	FILENAME_base = filepath;

	// minangle colormap mapping options
	double minangle_min_threshold = 0.;
	double minangle_max_threshold = 70.5;
	bool invert_minangle_cmap = true;
	double mindecentry_maxabs = 5;

	if (argc == 5) {
		minangle_min_threshold = std::stod(argv[2]);	
		minangle_max_threshold = std::stod(argv[3]);	
		if (!atoi(argv[4])) invert_minangle_cmap = false;
	}

	if (argc >= 6) {
		mindecentry_maxabs = std::stod(argv[5]);	
	}

	// #############
	// # LOAD MESH #
	// #############
	CGALTriangulation<Kernel> tri;
	std::vector<int> orbitinds;

	if(loadMeshWithOrbitpoints(tri, orbitinds, filepath)) {
		std::cout << "...loaded mesh with " << orbitinds.size() << " orbitinds" << std::endl;
	} else {
		std::cout << "Something went wrong loading the mesh" << std::endl;	
		return 0;
	}

	// #########################################
	std::cout << "METRICS"  << std::endl;
	// #########################################

	enum Metric {minangle=0, amips, volume, minentrydecl};
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
	metric_names[minentrydecl] = "minentrydecl";
	FILENAME = FILENAME_base + metric_names[metric_shown];

	Eigen::VectorXd Vol, Minang, Amips, Minentrydecl; 

	tri.calcMinDECEdgeContributionAllCells(Minentrydecl);
	tri.calcVolumeAllCells(Vol);
	tri.calcMinAngleAllCells(Minang);
	tri.calcAMIPSAllCells(Amips);

	cell_metrics[volume]=Vol;
	cell_metrics[minangle]=Minang;
	cell_metrics[amips]=Amips;
	cell_metrics[minentrydecl]=Minentrydecl;

	std::cout << "clcd" << std::endl;
	std::cout << "Minenrydecls: " << std::endl;
	double minval = std::numeric_limits<double>::max();
	for (int i = 0; i < cell_metrics[minentrydecl].size(); i++) {
		std::cout << cell_metrics[minentrydecl][i] << std::endl;	
		if (cell_metrics[minentrydecl][i] < minval) minval = cell_metrics[minentrydecl][i];
	}
	std::cout << "MINVAL: " << minval << std::endl;

	
	/*
	// verify that vertices and tets are in correct order
	std::cout << "Vertices: ";
    for(auto it = tri.mesh.vertices_begin(); it != tri.mesh.vertices_end(); ++it){
        if(it->info() != -1) {
			std::cout << it->info() << " ";	
		}
	}
	std::cout << std::endl;
	std::cout << "Tets: ";
    for(auto it = tri.mesh.cells_begin(); it != tri.mesh.cells_end(); ++it) {
        if(it->info() != -1) {
			std::cout << it->info() << " ";	
		}
	}
	std::cout << std::endl;
	*/
	
	std::cout << "Write Metrics to File " << std::endl;

	std::ofstream feil;
	feil.open(FILENAME_base + "metrics.csv");
	feil << metric_names[minangle] << "," << metric_names[amips] << "," << metric_names[volume] << "," << metric_names[minentrydecl] << std::endl;
	for(int i; i < cell_metrics[volume].size(); i++) {
		feil << cell_metrics[minangle](i) << "," << cell_metrics[amips](i) << "," << cell_metrics[volume](i) << "," << cell_metrics[minentrydecl](i) << std::endl;
	}
	feil.close();

	tri.write(FILENAME_base + "meshfile.tet");


	bool normalize=false;
	double amips_max = 100;
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

		// normalize the min dec entry to show only negative entries
		for (int i=0; i < cell_metrics[minentrydecl].size(); i++){
			double value = cell_metrics[minentrydecl][i];
			cell_metrics[minentrydecl][i] = (value < 0)? -value / mindecentry_maxabs : 0;
		}
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
	Eigen::MatrixXd cellcolors_minentrydecl; 
	igl::colormap(igl::COLOR_MAP_TYPE_PLASMA, cell_metrics[minentrydecl], normalize, cellcolors_minentrydecl);
	cellcolors[minentrydecl] = cellcolors_minentrydecl;

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
			ImGui::Combo("Metric", (int *)(&metric_shown), "MinAngle\0AMIPS\0Volume\0MINDECL\0");

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

	viewer.core().trackball_angle.x() = 0.0774649 ;
	viewer.core().trackball_angle.y() = 0.493061  ;
	viewer.core().trackball_angle.z() = 0.0441348 ;
	viewer.core().trackball_angle.w() = 0.865415  ;

	viewer.launch();
}
