
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
  double x2=x*x, y2=y*y, z2=z*z;
  double x4=x2*x2, y4=y2*y2, z4=z2*z2;

  return x4  + y4  + z4  + 2 *x2*  y2  + 2*
    x2*z2  + 2*y2*  z2  - 5 *x2  + 4* y2  - 5*z2+4;
}

template<class TKernel>
typename TKernel::FT
torus_fun (const typename TKernel::Point_3& p)
{ return tf(p.x(), p.y(), p.z());}

struct meshingOptions {
    meshingOptions() : cellSize(0.1),
					   cell_radius_edge_ratio(2.),
					   n_orbitpoints(-1), 
					   opt_lloyd(false),
					   opt_perturb(false), 
					   opt_exude(false)  {}
	double cellSize;
	double cell_radius_edge_ratio;
	int n_orbitpoints;
	bool opt_lloyd;
	bool opt_perturb;
	bool opt_exude;
};

template<class TKernel>
void 
meshSphere(IndexedTetMesh& indexed, meshingOptions mOptions)
{
    using namespace CGAL;
    using namespace CGAL::parameters;

	//typedef typename TKernel Kernel;
    typedef typename TKernel::FT FT;
    typedef typename TKernel::Point_3 Point;
	//typedef typename CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
	typedef CGAL::Mesh_domain_with_polyline_features_3< CGAL::Labeled_mesh_domain_3<TKernel> > Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default, CGAL::Sequential_tag>::type Tr;
	typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> TMesh;
    typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

	std::vector<Point> addPoints;
	if (mOptions.n_orbitpoints >= 0) {
		// add the origin 
		Point origin = Point(0., 0., 0.);
		addPoints.push_back(origin);

		// add orbit points (on a mini sphere inside the other one)
		Eigen::MatrixXd orbitpoints = randomPoints(mOptions.n_orbitpoints);
		double orbit_height = 0.5;
		for(int i=0; i< orbitpoints.cols(); i++){
			double vecnorm = orbitpoints.col(i).norm();
			addPoints.push_back(Point(orbitpoints(0,i) / vecnorm * orbit_height, 
									  orbitpoints(1,i) / vecnorm * orbit_height, 
									  orbitpoints(2,i) / vecnorm * orbit_height));
		}
	}

	//Mesh domain (sphere)
	Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain( 
											[](Point p)->double{return sphere_function<TKernel>(p, 1.);}, //torus_fun<TKernel>,
											typename TKernel::Sphere_3(CGAL::ORIGIN, 5.*5.));

	/*
	// adaptive sizing field:
	// Sizing field
	struct Spherical_sizing_field
	{
	  FT operator()(const Point& p, const int, const typename Mesh_domain::Index&) const
	  {
		FT sq_d_to_origin = CGAL::squared_distance(p, Point(-2, -2, -2));
		return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 5. + 0.025; 
	  }
	};
	*/
	typedef FT (Function)(const Point&);

    // Sizing field
    struct Spherical_sizing_field
    {
        typedef Point Point_3;
        typedef typename Mesh_domain::Index Index;

        FT operator()(const Point_3& p, const int, const Index&) const
        {
            FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
            return CGAL::abs( CGAL::sqrt(sq_d_to_origin)-0.5 ) / 5. + 0.025;
        }
    };

	//Spherical_sizing_field cellSize_field;

	// negative number of orbit points means not to add the origin eather
	//add origin and orbit points as 0-dim features (called corners)
	if (mOptions.n_orbitpoints >= 0) domain.add_corners(addPoints.begin(), addPoints.end());

	/*
	// Mesh criteria
	Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025, cell_radius_edge_ratio=mOptions.cell_radius_edge_ratio, cell_size=cellSize_field); // mOptions.cellSize);
	*/
    Spherical_sizing_field size;
	Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025,
                         cell_radius_edge_ratio=2, cell_size=mOptions.cellSize);
  

	// Mesh generation 
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, 
										no_perturb(), no_exude());
  
	if (mOptions.opt_lloyd)   CGAL::lloyd_optimize_mesh_3(c3t3, domain, time_limit=30);
	if (mOptions.opt_perturb) CGAL::perturb_mesh_3(c3t3, domain, time_limit=15);
	if (mOptions.opt_exude)   CGAL::exude_mesh_3(c3t3, sliver_bound=10, time_limit=10);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);

	/*
	std::cout << "Looking for the 0d features after opt" << std::endl;
	int found_addpoints = 0;
	bool found_origin   = false;
	for (auto a: indexed.vertices) {
		for (auto p: addPoints) {
			if (p.x() == a[0] && p.y() == a[1] && p.z() == a[2]){
				if (p.x() == 0. && p.y() == 0. and p.z() == 0.) {
					found_origin=true;
				} else {
					++found_addpoints;	
				}
			}
		}
	}
	std::cout << "Found " << found_addpoints << "/" << mOptions.n_orbitpoints << " orbitpoints" << std::endl;
	std::cout << (found_origin?"origin found":"ERROR: origin LOST") << std::endl;
	*/

}

template<class TKernel2, class TKernel>
void
meshSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions)
{
    IndexedTetMesh indexed;
    meshSphere<TKernel2>(indexed, mOptions);
    indexed.convert(tri);
}

