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

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
#include <Eigen/Dense>

int number_of_included_elements(std::vector<int> A, std::vector<int> B) {
	int fnd = 0;
	for (int el: B) {
		if(std::find(A.begin(), A.end(), el) != A.end()) {
			fnd++;
		}
	}
	return fnd;
}

void testSurfaceVertices(CGALTriangulation<Kernel> tri, std::vector<int> innerShell, std::vector<int> middleShell, std::vector<int> outerShell){
	std::vector<int> surfA, surfB;
	int nfv = 0;
	for(auto h : tri.mesh.finite_vertex_handles()) {
		nfv++;
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
	int fnd =  number_of_included_elements(surfA, surfB);
	std::cout << "A contains " << fnd << "/" << surfB.size() << " elements from B" << std::endl;

	std::cout << "NFV = " << nfv << std::endl;
	std::cout << "NV  = " << tri.mesh.tds().number_of_vertices() << std::endl;

	std::cout << "Inner.size = " << innerShell.size() << std::endl;
	std::cout << "Middle.size = " << middleShell.size() << std::endl;
	std::cout << "Outer.size = " << outerShell.size() << std::endl;

	std::cout << "Inner  in A " << number_of_included_elements(innerShell, surfA) << std::endl;
	std::cout << "Middle in A " << number_of_included_elements(middleShell, surfA) << std::endl;
	std::cout << "Outer  in A " << number_of_included_elements(outerShell, surfA) << std::endl;

	std::cout << "Inner  in B " << number_of_included_elements(innerShell, surfB) << std::endl;
	std::cout << "Middle in B " << number_of_included_elements(middleShell, surfB) << std::endl;
	std::cout << "Outer  in B " << number_of_included_elements(outerShell, surfB) << std::endl;
}

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


int boundary_example(int argc, char *argv[])
{
	// -------- ARG HANDLING ---------
	if (argc < 2) {
		std::cout << "usage: argv[0] YOURMESH.meshfile" << std::endl;
		return 0;
	}
	std::string filepath = argv[1];	
	CGALTriangulation<Kernel> tri;

	std::vector<int> innerShell;
	std::vector<int> middleShell;
	std::vector<int> outerShell;

	// tri.read(filepath);
	loadMeshWithShellIndices(tri, innerShell, middleShell, outerShell, filepath);
	testSurfaceVertices(tri, innerShell, middleShell, outerShell);
	return 0;
}

// --------------------------------------------------------------------------------------------------------------------------

#include "IndexedTetMesh.hpp"
typedef CGAL::Exact_predicates_inexact_constructions_kernel TKernel;
using namespace CGAL::parameters;
typedef FT_to_point_function_wrapper<typename TKernel::FT, typename TKernel::Point_3> Function;
typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function> Function_wrapper;
typedef typename Function_wrapper::Function_vector Function_vector;
typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef typename Mesh_criteria::Facet_criteria    Facet_criteria;
typedef typename Mesh_criteria::Cell_criteria     Cell_criteria;
typedef typename TKernel::Point_3 Point;

/*
namespace internal
{
    template<class Kernel, class TMesh>
    IndexedTetMesh
    extractIndexed(TMesh& mesh)
    {
        using namespace std;
        auto& tr = mesh.triangulation();
        map<typename TMesh::Vertex_handle, unsigned> idMap;
        IndexedTetMesh tm;
        unsigned cnt = 0;
        
        for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
        {
            idMap[it] = cnt++;
            tm.vertices.push_back({it->point().x(), it->point().y(), it->point().z() });
        }

		std::cout << "tr.number_of_vertices(): " << tr.number_of_vertices() << std::endl;
		std::cout << "tr.number_of_cells(): "    << tr.number_of_cells()    << std::endl;
        
		int ciccnt = 0;
        for(auto it =  mesh.cells_in_complex_begin(); it !=  mesh.cells_in_complex_end(); ++it)
        {

            if(idMap.find(it->vertex(0)) != idMap.end()
               && idMap.find(it->vertex(1)) != idMap.end()
               && idMap.find(it->vertex(2)) != idMap.end()
               && idMap.find(it->vertex(3)) != idMap.end())
            {
                tm.tets.push_back({
                    idMap[it->vertex(0)],
                    idMap[it->vertex(1)],
                    idMap[it->vertex(2)],
                    idMap[it->vertex(3)]
                });
				ciccnt++;
            }
        }
		std::cout << "transfered " << ciccnt << " cells from the complex to IndexedTetMesh" << std::endl;
        
        return tm;
    }
}
*/


C3t3 
meshEmbeddedDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions)
{

	// sphere_function (const typename Point& p, double rad)

	// Define functions
	Function f1( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 3.,  0., 0., 0. );} );
	Function f2( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1., 0., 0., 0. );} );
	Function f3( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.5,  0., 0., 0. );} );
	Function f4( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 2.,  0., 0., 0. );} );

	Function_vector v;
	v.push_back(f1);
	v.push_back(f2);
	v.push_back(f3);
	v.push_back(f4);

	//std::vector<std::string> vps;
	//vps.push_back("++-");
	//vps.push_back("+--");

	// DOMAIN
	Mesh_domain domain(function = Function_wrapper(v), //, vps),
					 // bounding_object = CGAL::Bbox_3(-3, -3, -3, 3, 3, 3), 
					 bounding_object = typename TKernel::Sphere_3(CGAL::ORIGIN, mOptions.boundingRad*mOptions.boundingRad),
					 relative_error_bound = 1e-6);

	// CRITERIA
	
	// Sizing field
	struct non_symmetric_sizing_field
	{
		non_symmetric_sizing_field(double cSize=0.2, double mSize=0.02) :	
										cellSize(cSize),
										minSize(mSize) {}
		double cellSize;
		double minSize;

		typedef typename TKernel::FT FT;
		typedef typename Mesh_domain::Index Index;

		FT operator()(const Point& p, const int, const Index&) const
		{
			// distance from  the plane orthogonal to y and going through (0 -1 0) - 
			return (p.y() + 2) * ((cellSize - minSize) / 4) + minSize;
		}
	};

	Facet_criteria facet_criteria(30, mOptions.facet_size, mOptions.approx_val); // angle, size, approximation
	std::cout << "use sizing field: " << mOptions.use_sizing_field << std::endl;
	//std::cout << "used: " (mOptions.use_sizing_field)? "ja":"nein" << std::endl;

	//non_symmetric_sizing_field size(mOptions.cell_size, mOptions.minSize);
	Cell_criteria cell_criteria(mOptions.cell_radius_edge_ratio, mOptions.cell_size); //  (mOptions.use_sizing_field)? size : mOptions.cell_size); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=10, time_limit=10);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);
}
// --------------------------------------------------------------------------------------------------------------------------


/*
int main (int argc, char *argv[]) 
{
    CGALTriangulation<Kernel> tri;

	bool singleSphere = false; // works for facet_size = 1., approx_val = 0.0067;
	bool embeddedDoubleSphere = false;


	meshingOptions mOptions;
	mOptions.cell_size              = std::stod(argv[1]);
	if (argc >=4) {
		if (atoi(argv[3])) embeddedDoubleSphere = true;	
		mOptions.facet_size             = std::stod(argv[4]); //1.; // 0.1 works with 0.02 approx val
		mOptions.approx_val             = std::stod(argv[5]); //0.0067;
	}
	mOptions.opt_lloyd        = true;
	mOptions.opt_perturb      = false;
	mOptions.opt_exude        = true;
	mOptions.use_sizing_field = false;

	std::cout << "######################" << std::endl;
	std::cout << "    MESH GENERATION   " << std::endl;
	std::cout << "######################" << std::endl;

	IndexedTetMesh indexed;
	C3t3 c3t3 = meshEmbeddedDoubleSphere(indexed, mOptions);
	indexed.convert(tri);
	
	Eigen::SparseMatrix<double> La, Lb, Ma, Mb;

	//tri.DECLaplacianRegular(tri, La,&Ma);
	tri.DECLaplacianRegular(c3t3,Lb,&Mb);

	std::cout << (La - Lb).norm() << std::endl;
}
*/

int main (int argc, char *argv[]) {
	std::cout << "commented out to compile" << std::endl;
}
