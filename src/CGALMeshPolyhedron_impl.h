
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Implicit_to_labeling_function_wrapper.h>

#include "IndexedTetMesh.hpp"
#include "CGALTriangulation.hpp"
#include "Timer.hpp"


#include <unordered_map>

#include <cxxabi.h>
 
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

template<class TPolyhedron>
void
meshPolyhedron(const TPolyhedron& p, IndexedTetMesh& indexed, const double cellSize)
{
    using namespace CGAL;
    using namespace CGAL::parameters;
    
    typedef typename TPolyhedron::Kernel Kernel;
    typedef typename Kernel::Point_3 Point3;
    typedef typename CGAL::Polyhedral_mesh_domain_3<typename TPolyhedron::Polyhedron, Kernel> Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default, CGAL::Sequential_tag>::type Tr;
    typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> TMesh;
    
    // Create domain
    Mesh_domain domain(p.poly);
    typedef  CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    
    // Mesh criteria (no cell_size set)
    Mesh_criteria criteria(facet_size=0.05, facet_distance=0.008, cell_size=cellSize);
    
    // Mesh generation
    auto c3t3 = CGAL::make_mesh_3<TMesh>(domain, criteria, no_perturb(), no_exude());

    indexed = ::internal::extractIndexed<Kernel>(c3t3);
}

template<class TPolyhedron, class TKernel>
void
meshPolyhedron(const TPolyhedron& p, CGALTriangulation<TKernel>& tri, const double cellSize)
{
    IndexedTetMesh indexed;
    meshPolyhedron(p, indexed, cellSize);
    indexed.convert(tri);
}

// ################################### SPHERE #############################################
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

template<class TKernel>
typename TKernel::FT
sphere_function (const typename TKernel::Point_3& p, double rad)
{ return CGAL::squared_distance(p, typename TKernel::Point_3(CGAL::ORIGIN))-rad*rad; }


double tf (double x, double y, double z) {
	double r = .5;
	double R = 1;
	double x2=x*x, y2=y*y, z2=z*z;
	double r2=r*r;
	//double x4=x2*x2, y4=y2*y2, z4=z2*z2;
	// return x4  + y4  + z4  + 2 *x2*  y2  + 2* x2*z2  + 2*y2*  z2  - 5 *x2  + 4* y2  - 5*z2+4;
	
	return pow(sqrt(x2 + y2) - R, 2) + z2 - r2;
}

double implicit_sphere_function(double x, double y, double z, double r, double mx, double my, double mz) {
	return (x-mx)*(x-mx) + (y-my)*(y-my) + (z-mz)*(z-mz) - r*r;
}

/*
double sf_outer (double x, double y, double z) {
	double r = 3;
	return x*x + y*y + z*z - r*r;
}
double sf_inner (double x, double y, double z) {
	double r = 0.01;
	return x*x + y*y + z*z - r*r;
}

template<double rad, double mx, double my, double mz> 
double
impl_circle_fun(double x double y double z){
	return (x-mx)*(x-mx) + (y-my)*(y-my) + (z-mz)*(z-mz) - rad*rad;
}
*/

template<class TKernel>
typename TKernel::FT
torus_fun (const typename TKernel::Point_3& p)
{ return tf(p.x(), p.y(), p.z());}

struct meshingOptions {
    meshingOptions() : cell_size(0.1),
					   cell_radius_edge_ratio(2.),
					   facet_size(1.),
					   approx_val(0.02),
					   boundingRad(5.),
					   minSize(0.02),
					   opt_lloyd(false),
					   opt_perturb(false), 
					   opt_exude(false),
					   use_sizing_field(false) {}
	double cell_size;
	double cell_radius_edge_ratio;
	double facet_size;
	double approx_val;
	double boundingRad;
	double minSize;
	bool opt_lloyd;
	bool opt_perturb;
	bool opt_exude;
	bool use_sizing_field;

};

template<class TKernel>
void 
meshSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath)
{
    using namespace CGAL;
    using namespace CGAL::parameters;

	//typedef typename TKernel Kernel;
    typedef typename TKernel::FT FT;
    typedef typename TKernel::Point_3 Point;
	typedef FT (Function)(const Point&);

	//typedef CGAL::Mesh_domain_with_polyline_features_3< CGAL::Labeled_mesh_domain_3<TKernel> > Mesh_domain;
	typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default, CGAL::Sequential_tag>::type Tr;
	typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> TMesh;
    typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

	//Mesh_domain domain;
	//if (!use_torus) {
		// SPHERE
	Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(
										[](Point p)->double{return sphere_function<TKernel>(p, 1.);},
										typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.)
										);

	Mesh_criteria criteria(facet_angle=30, facet_size=mOptions.facet_size, facet_distance=0.025,
                         cell_radius_edge_ratio=mOptions.cell_radius_edge_ratio, cell_size=mOptions.cell_size);

	// Mesh generation 
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
										no_perturb(), no_exude());
  
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=0, time_limit=0); //sliver_bound=10, time_limit=10);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);

	if (!regweightsoutpath.empty()) {
		std::ofstream regweightsfile;
		regweightsfile.open(regweightsoutpath);
		regweightsfile << "regweight" << std::endl;
        auto& tr = c3t3.triangulation();
        for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
        {
			regweightsfile << it->point().weight() << std::endl; 
        }
		regweightsfile.close();
	}
}


template <int Sq_radius>
double sphere_function (double x, double y, double z) // (c=(0,0,0), r=Sq_radius)
{
  double x2=x*x, y2=y*y, z2=z*z;
  return (x2+y2+z2)/Sq_radius - 1;
}

template <typename FT, typename P>
class FT_to_point_function_wrapper : public CGAL::cpp98::unary_function<P, FT>
{
  typedef FT (*Implicit_function)(FT, FT, FT);
  Implicit_function function;
public:
  typedef P Point;
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
};

template<class TKernel>
void 
meshDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath, std::string decreglaplacianoutpath)
{
	using namespace CGAL::parameters;
	// Domain
	typedef FT_to_point_function_wrapper<typename TKernel::FT, typename TKernel::Point_3> Function;
	typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
															Function_wrapper;
	typedef typename Function_wrapper::Function_vector Function_vector;
	typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
	// Triangulation
	typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
	// Mesh Criteria
	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
	typedef typename Mesh_criteria::Facet_criteria    Facet_criteria;
	typedef typename Mesh_criteria::Cell_criteria     Cell_criteria;

	typedef typename TKernel::Point_3 Point;

	// sphere_function (const typename Point& p, double rad)

	// Define functions
	Function f1( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.,  0., 0., 0. );} );
	Function f2( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.5, 0., 0., 0. );} );
	Function f3( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 2.,  0., 0., 0. );} );

	Function_vector v;
	v.push_back(f1);
	v.push_back(f2);
	v.push_back(f3);

	std::vector<std::string> vps;
	vps.push_back("++-");
	vps.push_back("+--");

	// DOMAIN
	Mesh_domain domain(function = Function_wrapper(v, vps),
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
	
	//non_symmetric_sizing_field size(mOptions.cell_size, mOptions.minSize);
	//std::cout << "USE SIZING FILED: " << ((mOptions.use_sizing_field) ? "Ja" : "Nein")  << std::endl;
	//auto sz =  (mOptions.use_sizing_field)? size : mOptions.cell_size;
	Cell_criteria cell_criteria(mOptions.cell_radius_edge_ratio, mOptions.cell_size); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=0, time_limit=0); // sliver_bound=10, time_limit=10

    indexed = ::internal::extractIndexed<TKernel>(c3t3);

	if (!regweightsoutpath.empty()) {
		std::ofstream regweightsfile;
		regweightsfile.open(regweightsoutpath);
		regweightsfile << "regweight" << std::endl;
        auto& tr = c3t3.triangulation();
        for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
        {
			regweightsfile << it->point().weight() << std::endl; 
        }
		regweightsfile.close();
	}

	// test the laplace calc
	if (!decreglaplacianoutpath.empty()) {
		std::cout << "calc dec laplacian regular from c3t3" << std::endl;	
		Eigen::SparseMatrix<double> L, M; 

		/*
		int status;
		char* realname;
		realname = abi::__cxa_demangle(c3t3.name(), 0, 0, &status);
		std::cout << c3t3.name() << "\t=> " << realname << "\t: " << status << '\n';
		free(realname);
		*/
		
		calcDECLaplacianRegularFromC3t3<TKernel>(c3t3,L,&M);
		Eigen::saveMarket(L, decreglaplacianoutpath + "decregLoptimized.mtx");
		Eigen::saveMarket(M, decreglaplacianoutpath + "decregMoptimized.mtx");
	}
}

template<class TKernel>
void 
meshEmbeddedDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath)
{
	using namespace CGAL::parameters;
	// Domain
	typedef FT_to_point_function_wrapper<typename TKernel::FT, typename TKernel::Point_3> Function;
	typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
															Function_wrapper;
	typedef typename Function_wrapper::Function_vector Function_vector;
	typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
	// Triangulation
	typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
	// Mesh Criteria
	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
	typedef typename Mesh_criteria::Facet_criteria    Facet_criteria;
	typedef typename Mesh_criteria::Cell_criteria     Cell_criteria;

	typedef typename TKernel::Point_3 Point;

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
	//std::cout << "use sizing field: " << mOptions.use_sizing_field << std::endl;
	//std::cout << "used: " (mOptions.use_sizing_field)? "ja":"nein" << std::endl;

	//non_symmetric_sizing_field size(mOptions.cell_size, mOptions.minSize);
	Cell_criteria cell_criteria(mOptions.cell_radius_edge_ratio, mOptions.cell_size); //  (mOptions.use_sizing_field)? size : mOptions.cell_size); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=0, time_limit=0); //sliver_bound=10, time_limit=10);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);

	if (!regweightsoutpath.empty()) {
		std::ofstream regweightsfile;
		regweightsfile.open(regweightsoutpath);
		regweightsfile << "regweight" << std::endl;
        auto& tr = c3t3.triangulation();
        for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
        {
			regweightsfile << it->point().weight() << std::endl; 
        }
		regweightsfile.close();
	}
}


template<class TKernel>
void 
meshSingleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath)
{
	using namespace CGAL::parameters;
	// Domain
	typedef FT_to_point_function_wrapper<typename TKernel::FT, typename TKernel::Point_3> Function;
	typedef CGAL::Implicit_multi_domain_to_labeling_function_wrapper<Function>
															Function_wrapper;
	typedef typename Function_wrapper::Function_vector Function_vector;
	typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
	// Triangulation
	typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
	// Mesh Criteria
	typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
	typedef typename Mesh_criteria::Facet_criteria    Facet_criteria;
	typedef typename Mesh_criteria::Cell_criteria     Cell_criteria;

	typedef typename TKernel::Point_3 Point;

	// sphere_function (const typename Point& p, double rad)

	// Define functions
	Function f1( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.,  0., 0., 0. );} );
	Function f2( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 2., 0., 0., 0. );} );
	//Function f3( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.5,  0., 0., 0. );} );

	std::cout << "fun test:" << std::endl;
	Point p1(0., 0., 0.1), p2(0., 0., 0.6), p3(0., 0., 1.1), p4(0., 0., 1.6);
	
	/*
	std::cout << " ORIGIN : " << std::endl;
	std::cout << CGAL::ORIGIN << std::endl;
	*/

	Function_vector v;
	v.push_back(f1);
	v.push_back(f2);
	//v.push_back(f3);
	std::vector<std::string> vps;
	//vps.push_back("+-");
	//vps.push_back("--");
	//vps.push_back("---");

	auto v_wrapped = Function_wrapper(v);//, vps);
	std::cout << v_wrapped(p1) << ", " << v_wrapped(p2) << ", " << v_wrapped(p3) << "," << v_wrapped(p4) << std::endl;


	// DOMAIN
	Mesh_domain domain(function = v_wrapped,  //Function_wrapper(v), //, vps),
					 // bounding_object = CGAL::Bbox_3(-3, -3, -3, 3, 3, 3), 
					 bounding_object = typename TKernel::Sphere_3(Point(0., 0., 0.), mOptions.boundingRad*mOptions.boundingRad),
					 relative_error_bound = 1e-6);
	//marc:
	//Mesh_domain domain(Function_wrapper(v,vps), typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);

	// CRITERIA
	//marc:
	//Facet_criteria facet_criteria(30, 0.02, 0.005); // angle, size, approximation
	Facet_criteria facet_criteria(30, mOptions.facet_size, mOptions.approx_val); // angle, size, approximation
	//marc:
	//Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
	Cell_criteria cell_criteria(mOptions.cell_radius_edge_ratio, mOptions.cell_size); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3,sliver_bound=0, time_limit=0); //sliver_bound=10, time_limit=10);
	
	double maxdist = 0;
	for (auto vh : c3t3.triangulation().finite_vertex_handles()) {
		double dist = 	sqrt(CGAL::squared_distance(vh->point(), typename TKernel::Point_3(CGAL::ORIGIN)));
		if (dist > maxdist) {
			maxdist = dist;
		}
	}
	std::cout << "Maxdist: " << maxdist << std::endl;

    indexed = ::internal::extractIndexed<TKernel>(c3t3);

	if (!regweightsoutpath.empty()) {
		std::ofstream regweightsfile;
		regweightsfile.open(regweightsoutpath);
		regweightsfile << "regweight" << std::endl;
        auto& tr = c3t3.triangulation();
        for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
        {
			regweightsfile << it->point().weight() << std::endl; 
        }
		regweightsfile.close();
	}

}



template<class TKernel2, class TKernel>
void
meshSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath)
{
    IndexedTetMesh indexed;
    meshSphere<TKernel2>(indexed, mOptions, regweightsoutpath);
    indexed.convert(tri);
}


template<class TKernel2, class TKernel>
void
meshDoubleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath, std::string decreglaplacianoutpath)
{
    IndexedTetMesh indexed;
    meshDoubleSphere<TKernel2>(indexed, mOptions, regweightsoutpath, decreglaplacianoutpath);
    indexed.convert(tri);
}

template<class TKernel2, class TKernel>
void
meshEmbeddedDoubleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath)
{
    IndexedTetMesh indexed;
    meshEmbeddedDoubleSphere<TKernel2>(indexed, mOptions, regweightsoutpath);
    indexed.convert(tri);
}

template<class TKernel2, class TKernel>
void
meshSingleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath)
{
    IndexedTetMesh indexed;
    meshSingleSphere<TKernel2>(indexed, mOptions, regweightsoutpath);
    indexed.convert(tri);
}


template<class TKernel>
void
calcDECLaplacianRegularFromC3t3(
		//CGAL::Mesh_complex_3_in_triangulation_3<CGAL::Mesh_triangulation_3<CGAL::Labeled_mesh_domain_3<TKernel>>> c3t3,
		C3T3T<TKernel> c3t3,
		Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M)
{

	typedef CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
	typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
	typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

	typedef typename TKernel::Point_3 Point;
	typedef typename TKernel::Plane_3 Plane;
	typedef typename TKernel::Line_3  Line;
	typedef typename TKernel::Segment_3  Segment;

    std::vector<Eigen::Triplet<double>> triplets;
    const int nv = c3t3.triangulation().number_of_vertices();

    Eigen::MatrixXd V(nv, 3);

	// id Map to  store to get the vertex indices withouth a ->info() construct
    auto& tr = c3t3.triangulation();
	std::map<typename C3t3::Vertex_handle, unsigned> idMap;
	unsigned cnt = 0;
	for(auto it = tr.finite_vertices_begin(); it != tr.finite_vertices_end(); ++it)
	{
		V(cnt, 0) = it->point().x();
		V(cnt, 1) = it->point().y();
		V(cnt, 2) = it->point().z();
		idMap[it] = cnt++;
	}
    
    // turn off some costly sanity tests
    bool dbg = true;
	bool dbg2 = false;
    bool dbg3 = false;
    
    std::vector<typename TKernel::Vector_3> vecs(nv);
    
    if(M)
    {
        M->resize(nv, nv);
        M->resizeNonZeros(nv);
        
        for(int i = 0; i < nv; ++i)
        {
            M->outerIndexPtr()[i] = i;
            M->innerIndexPtr()[i] = i;
            M->valuePtr()[i] = .0;
        }
        
        M->outerIndexPtr()[nv] = nv;
    }
    
    //for(auto h : tr.finite_cell_handles())
    for(auto h =  c3t3.cells_in_complex_begin(); h !=  c3t3.cells_in_complex_end(); ++h)
    {
        auto tet = tr.tetrahedron(h);
        double vol = tet.volume();
        
        for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            if( i != j )
            {
                const int k = Tr::next_around_edge(i, j);
                const int l = Tr::next_around_edge(j, i);

				auto tet_dual  = tr.dual(h);
				auto face_dual_res = tr.dual(h, l); // facet (tet[i], tet[j], tet[k]);
				Segment face_dual;
				assign(face_dual, face_dual_res);

				// auto ccf = CGAL::circumcenter(tet[i], tet[j], tet[k]);

				Line  fd_line   = face_dual.supporting_line();
				Plane tri_plane(tet[i], tet[j], tet[k]);
				auto ccf_res = intersection(tri_plane, fd_line);
				auto ccf = boost::get<Point>(&*ccf_res);

				// auto cce = CGAL::circumcenter(tet[i], tet[j]);
				auto edge = tet[j] - tet[i];
				Line edge_line(tet[i], edge);
				Plane dual_plane(tet_dual, edge);
				//auto cce = intersection(dual_plane, edge_line);
				
				// access point directly since intersection cannot be a line in this case
				auto cce_res = intersection(dual_plane, edge_line);
				auto cce = boost::get<Point>(&*cce_res);

				auto edge_normalized = edge / sqrt(edge.squared_length());
				Point cce_const = CGAL::ORIGIN + (edge_normalized * (tet_dual - tet[i])) * edge_normalized;

				//std::cout << "Calced" << std::endl;

				// ----------------------------------------------
				// construct ccf hopefully correctly?
				double tri_area = 0.5 * sqrt(CGAL::cross_product(tet[j] - tet[i], tet[k] - tet[i]).squared_length());

				Point ccf_const = tet[i]; 
				// contrib of point k
				auto e_k = (tet[j] - tet[i]) / sqrt((tet[j]-tet[i]).squared_length());
				auto normal_k   = tet[k] - (tet[i] + e_k * ((tet[k]- tet[i]) * e_k));
				normal_k = normal_k / sqrt(normal_k.squared_length()) * sqrt((tet[j] - tet[i]).squared_length());
				ccf_const += (((tet[k] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(k)->point().weight()) * normal_k) / (4. * tri_area) ; 
				//contrib of point j
				auto e_j = (tet[k] - tet[i]) / sqrt((tet[k]-tet[i]).squared_length());
				auto normal_j   = tet[j] - (tet[i] + e_j * ((tet[j]- tet[i]) * e_j));
				normal_j = normal_j / sqrt(normal_j.squared_length()) * sqrt((tet[k] - tet[i]).squared_length());
				ccf_const += (((tet[j] - tet[i]).squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * normal_j) / (4. * tri_area); 

				//auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - cce_const, tet_dual - cce_const); 
				auto nrml_const = 0.5 * CGAL::cross_product(ccf_const - *cce, tet_dual - *cce);
				double val = (nrml_const * edge) / edge.squared_length();
				// ----------------------------------------------
                const int r = idMap[h->vertex(i)]; // h->vertex(i)->info();
                const int s = idMap[h->vertex(j)];// h->vertex(j)->info();
            
                if(M)
                {
					// TODO: eddge.squared_length the correct thing here?
                    M->valuePtr()[r] += val * (edge.squared_length()) / 6.;
                    M->valuePtr()[s] += val * (edge.squared_length()) / 6.;
                }

				if (dbg2) {
					// compare with non-weighte circumcenter calls (only makes sense if weights are all 0)
                    auto cc_nr  = CGAL::circumcenter(tet);
                    auto ccf_nr = CGAL::circumcenter(tet[i], tet[j], tet[k]);
                    auto cce_nr = CGAL::circumcenter(tet[i], tet[j]);
            
                    auto nrml_nr = 0.5 * CGAL::cross_product(ccf_nr - cce_nr, cc_nr - cce_nr);
                    double val2 = (nrml_nr * edge) / edge.squared_length();

					Point cce_const = tet[i] + 0.5 * (edge.squared_length() + h->vertex(i)->point().weight() - h->vertex(j)->point().weight()) * (edge / sqrt(edge.squared_length()));

					if (std::abs(val - val2) > 1e-13) {

						std::cout << std::endl << "Val  : " << val << std::endl;
						std::cout << "Val2 : " << val2 << std::endl;
						std::cout << "CC  - CC_nr  = " << tet_dual - cc_nr << std::endl;
						std::cout << "CCF - CCF_nr = " << *ccf - ccf_nr << std::endl;
						std::cout << "CCE - CCE_nr = " << *cce - cce_nr << std::endl;

						std::cout << std::endl;
						std::cout << "CCF      : " << *ccf      << std::endl;
						std::cout << "CCF_nr   : " << ccf_nr    << std::endl;
						std::cout << "ccf_const: " << ccf_const << std::endl;

						std::cout << std::endl;
						std::cout << "CCE       : " << *cce      << std::endl;
						std::cout << "CCE_nr    : " << cce_nr    << std::endl;
						std::cout << "CCE_const : " << cce_const << std::endl;

					}

                }

				if (dbg3) {
					// compare to phillipps dec impl (only makes sense if weights are all 0)
					
					auto a2 = tet[i] - tet[l];
					auto b = tet[j] - tet[l];
					auto c = tet[k] - tet[l];
					
					auto b2 = tet[i] - tet[k];
					auto c2 = tet[j] - tet[i];
					auto a =  tet[k] - tet[j];
					
					 
					const double n = (a * b2) * (b2 * c2) * (c * a2) + (b2 * c2) * (c2 * a) * (a2 * b)
					+ (c2 * a) * (a * b2) * (b * c) + (a * b2) * (b2 * c2) * (c2 * a);
				  
					const double d = 192 * vol * 0.25 * CGAL::cross_product(b2, a).squared_length();
					const double fac = std::abs(d) < 1e-24 ? .0 : 1. / d;
					const double val3 = -n * (a * b2) * fac;

					
					if (std::abs(val - val3) > 1e-13) {
						std::cout << "ERROR: " << val << "!= " << val3 << std::endl;
					}
				}

				if (r > nv) {
					std::cout << "r " << r << " OUT!!!!!!" << std::endl;	
				}
				if (s > nv) {
					std::cout << "s " << s << " OUT!!!!!!" << std::endl;	
				}

                triplets.emplace_back(r, r, -val);
                triplets.emplace_back(s, s, -val);
                    
                triplets.emplace_back(r, s, val);
                triplets.emplace_back(s, r, val);

				/*
				if(std::abs(val - val2) > 1e-10) std::cout << "error: " << val << " " << val2 << std::endl;
			  
				vecs[r] += nrml;
				vecs[s] -= nrml;
				*/

            }
    }
   
    L.resize(nv, nv);
    L.setFromTriplets(triplets.begin(), triplets.end());

}
