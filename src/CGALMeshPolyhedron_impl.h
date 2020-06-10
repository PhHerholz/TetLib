
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
            }
        }
        
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
					   facet_size(0.1),
					   approx_val(0.005),
					   boundingRad(5.),
					   opt_lloyd(false),
					   opt_perturb(false), 
					   opt_exude(false)  {}
	double cell_size;
	double cell_radius_edge_ratio;
	double facet_size;
	double approx_val;
	double boundingRad;
	bool opt_lloyd;
	bool opt_perturb;
	bool opt_exude;
};

template<class TKernel>
void 
meshSphere(IndexedTetMesh& indexed, meshingOptions mOptions, bool use_torus)
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

	std::cout << "USE TORUS: " << use_torus << std::endl;

	//Mesh_domain domain;
	//if (!use_torus) {
		// SPHERE
	Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(
										[](Point p)->double{return sphere_function<TKernel>(p, 1.);},
										typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.)
										);
	/*
	if (use_torus) {
		domain = Mesh_domain::create_implicit_mesh_domain(
											torus_fun<TKernel>,
											typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.));
	}
	*/

	/*
	std::cout << "Using a sizing field" << std::endl;
    // Sizing field
    struct Spherical_sizing_field
    {
		//typedef ::FT FT;
		typedef Point Point_3;
        typedef typename Mesh_domain::Index Index;

        FT operator()(const Point_3& p, const int, const Index&) const
        {
            FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
            return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 5. + 0.025;
        }
    };

    Spherical_sizing_field size;
	Mesh_criteria criteria(facet_angle=30, facet_size=mOptions.facet_size, facet_distance=0.025,
                         cell_radius_edge_ratio=mOptions.cell_radius_edge_ratio, cell_size=size);
						 */

	Mesh_criteria criteria(facet_angle=30, facet_size=mOptions.facet_size, facet_distance=0.025,
                         cell_radius_edge_ratio=mOptions.cell_radius_edge_ratio, cell_size=mOptions.cell_size);

	// Mesh generation 
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
										no_perturb(), no_exude());
  
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=10, time_limit=10);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);
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
meshDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions)
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

	//Function f1( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 1.0, 0., 0., 0. );} );
	//Function f2( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 0.75, 0., 0., 0. );} );
	//Function f3( [](double x, double y, double z) -> double{return implicit_sphere_function(x, y, z, 0.5, 0., 0., 0. );} );

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
	//if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	//if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	//if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=10, time_limit=10);


	/*
	// marcs code: (works as expected)
	Function f1(&sphere_function<1>);
	Function f2(&sphere_function<2>);
	Function f3(&sphere_function<3>);
	Function_vector v;
	v.push_back(f1);
	v.push_back(f2);
	v.push_back(f3);
	std::vector<std::string> vps;
	vps.push_back("++-");
	vps.push_back("+--");
	// Domain (Warning: Sphere_3 constructor uses square radius !)
	Mesh_domain domain(Function_wrapper(v,vps), typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);
	// Set mesh criteria
	Facet_criteria facet_criteria(30, 0.02, 0.005); // angle, size, approximation
	Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);
	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());
	// Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
	CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
	// Exudation
	CGAL::exude_mesh_3(c3t3,12);
	*/

    indexed = ::internal::extractIndexed<TKernel>(c3t3);
}


template<class TKernel2, class TKernel>
void
meshSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, bool use_torus)
{
    IndexedTetMesh indexed;
    meshSphere<TKernel2>(indexed, mOptions, use_torus);
    indexed.convert(tri);
}


template<class TKernel2, class TKernel>
void
meshDoubleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions)
{
    IndexedTetMesh indexed;
    meshDoubleSphere<TKernel2>(indexed, mOptions);
    indexed.convert(tri);
}

