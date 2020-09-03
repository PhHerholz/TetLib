#pragma once

#include "CGALTriangulation.hpp"

template<class TPolyhedron>
void
meshPolyhedron(const TPolyhedron& p, IndexedTetMesh& indexed, const double cellSize = 0.05);

template<class TPolyhedron, class TKernel>
void
meshPolyhedron(const TPolyhedron& p, CGALTriangulation<TKernel>& tri, const double cellSize = 0.05);

struct meshingOptions;

template<class TKernel>
typename TKernel::FT
sphere_function (const typename TKernel::Point& p);

template<class TKernel>
void
meshSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="");

template<class TKernel2, class TKernel>
void 
meshSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="");


// DOUBLESPHERE:
template<class TKernel>
void 
meshDoubleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="");

template<class TKernel2, class TKernel>
void
meshDoubleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="");

// SINGLESPHERE (contains origin and a middleShell):
template<class TKernel>
void 
meshSingleSphere(IndexedTetMesh& indexed, meshingOptions mOptions, std::string regweightsoutpath="");
template<class TKernel2, class TKernel>
void
meshSingleSphere(CGALTriangulation<TKernel>& tri, meshingOptions mOptions, std::string regweightsoutpath="");

#include "CGALMeshPolyhedron_impl.h"
