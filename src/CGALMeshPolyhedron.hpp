#pragma once

#include "CGALTriangulation.hpp"


// only for c3t3
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

template<class TPolyhedron>
void
meshPolyhedron(const TPolyhedron& p, IndexedTetMesh& indexed, const double cellSize = 0.05);

template<class TPolyhedron, class TKernel>
void
meshPolyhedron(const TPolyhedron& p, CGALTriangulation<TKernel>& tri, const double cellSize = 0.05);


// ###################
// # SPHERE MESHING: #
// ###################

struct meshingOptions;

template<class TKernel>
typename TKernel::FT
sphere_function (const typename TKernel::Point& p);

// SPHERE
template<class TKernel>
void
meshSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="");
template<class TKernel2, class TKernel>
void 
meshSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="");

// DOUBLESPHERE:
template<class TKernel>
void 
meshDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="", std::string decreglaplacianoutpath="");
template<class TKernel2, class TKernel>
void
meshDoubleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="", std::string decreglaplacianoutpath="");

// SINGLESPHERE (contains origin and a middleShell):
template<class TKernel>
void 
meshSingleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="");
template<class TKernel2, class TKernel>
void
meshSingleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="");

// DEC LAPLACIAN with from power diagram of a c3t3
template<class TKernel>
void 
calcDECLaplacianRegularFromC3t3(
		//CGAL::Mesh_complex_3_in_triangulation_3<CGAL::Mesh_triangulation_3<CGAL::Labeled_mesh_domain_3<TKernel>>> c3t3,
		CGAL::Mesh_complex_3_in_triangulation_3<CGAL::Mesh_3_regular_triangulation_3_wrapper<CGAL::Robust_weighted_circumcenter_filtered_traits_3<CGAL::Epick>, CGAL::Triangulation_data_structure_3<CGAL::Mesh_vertex_base_3<CGAL::Robust_weighted_circumcenter_filtered_traits_3<CGAL::Epick>, CGAL::Labeled_mesh_domain_3<CGAL::Epick, int, std::pair<int, int> >, CGAL::Regular_triangulation_vertex_base_3<CGAL::Robust_weighted_circumcenter_filtered_traits_3<CGAL::Epick>, CGAL::Triangulation_ds_vertex_base_3<void> > >, CGAL::Compact_mesh_cell_base_3<CGAL::Robust_weighted_circumcenter_filtered_traits_3<CGAL::Epick>, CGAL::Labeled_mesh_domain_3<CGAL::Epick, int, std::pair<int, int> >, void>, CGAL::Sequential_tag> >, int, int> c3t3,
		Eigen::SparseMatrix<double>& L, Eigen::SparseMatrix<double>* M);

#include "CGALMeshPolyhedron_impl.h"
