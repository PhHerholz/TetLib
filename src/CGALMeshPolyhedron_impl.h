
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>

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


template<class TKernel>
typename TKernel::FT
sphere_function (const typename TKernel::Point& p)
{ return CGAL::squared_distance(p, typename TKernel::Point(CGAL::ORIGIN))-1; }

template<class TKernel>
void 
meshSphere(IndexedTetMesh& indexed, const double cellSize)
{
    using namespace CGAL;
    using namespace CGAL::parameters;

	//typedef typename TKernel Kernel;
    typedef typename TKernel::FT FT;
    typedef typename TKernel::Point_3 Point;
	typedef typename CGAL::Labeled_mesh_domain_3<TKernel> Mesh_domain;
    typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default, CGAL::Sequential_tag>::type Tr;
	typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
    typedef typename CGAL::Mesh_complex_3_in_triangulation_3<Tr> TMesh;
    typedef typename CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
    
	//Mesh domain (sphere)
	Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(sphere_function,
                                             typename TKernel::Sphere_3(CGAL::ORIGIN, 2.));
	// Mesh criteria
	Mesh_criteria criteria(facet_angle=30, facet_size=0.1, facet_distance=0.025, cell_radius_edge_ratio=2, cell_size=0.1);
  
	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    indexed = ::internal::extractIndexed<TKernel>(c3t3);
}

template<class TKernel>
void
meshSphere(CGALTriangulation<TKernel>& tri, const double cellSize)
{
    IndexedTetMesh indexed;
    meshSphere<TKernel>(indexed, cellSize);
    indexed.convert(tri);
}

